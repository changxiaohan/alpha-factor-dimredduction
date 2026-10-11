# =============================================================================
# _factor_ext2.R -- factor pool extension (Step 1): 12 -> 25 factors
# -----------------------------------------------------------------------------
# WHY: _diag_dim.R showed that the 12-factor pool has an effective size of only
#      ~3.2 independent factors, and an ORACLE-weighted 12-factor composite
#      (ICIR 0.548) barely beats a single factor (volatility, 0.522). Dimension
#      reduction has nothing to work with until the pool carries more
#      information. So: widen the pool first, reduce later.
#
# OUTPUT: _factor_panel2.parquet  (date, code, <all 25 factors>, already
#         winsorized + cross-sectionally z-scored by date)
#         The OLD _factor_panel.parquet is NOT touched, so results stay
#         comparable with the previous round.
#
# WHAT IS NEW (industry-standard constructs; NOT taught in any course):
#   price/volume : reversal21, ivol60, beta60, skew60, maxret21, amihud60,
#                  turnchg, hlrange60, vpcorr60
#   financial    : sp (revenue/mcap), ep_np (net profit/mcap), bps_yoy
#   KEPT AS-IS   : momentum, volatility, size, liquidity  (old definitions)
#
# FIELD SEMANTICS -- verified by _probe_units.R on 2025 data, do NOT guess:
#   volume            = traded shares
#   turnover          = turnover RATE as a FRACTION (0.022 = 2.2%), not a volume
#   outstanding_share = free-float shares  (= volume / turnover)
#   amount / volume   = VWAP price (CNY/share)     log-corr with close = 0.99
#   amount / turnover = free-float market cap (CNY) == outstanding_share * price
#   => `size` = log(amount/turnover) = log free-float mkt cap -- which is what
#      the OLD panel called `size` too. It was RIGHT all along; only its inline
#      comment there ("不复权价") mislabelled it.
#   => the OLD value_bp / value_ep divided per-share figures (bps, eps) by
#      MARKET CAP instead of by PRICE -> off by a factor of share count.
#      This script uses px = amount/volume for those two ratios.
#
# HARD RULES carried over from the previous round:
#   - never use frollapply (it returns fewer elements than n when n < window)
#   - every rolling stat uses align="right" + shift()  => no lookahead
#   - financial availability uses the LEGAL disclosure deadline (PIT)
#   - console output is ASCII only (Windows GBK console)
#
# Usage (run from the project directory):
#   $env:DFZQ_PANEL_DAYS="400"     # smoke, ~1 min
#   $env:DFZQ_PANEL_DAYS="0"       # full   (0 = all days; do NOT rely on
#                                  #         deleting the line -- the env var
#                                  #         persists in the shell session)
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_factor_ext2.R" `
#       | Tee-Object "_factor_ext2.log"
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")
el <- function() as.numeric(difftime(Sys.time(), T0, units = "secs"))

T0 <- Sys.time()
W  <- 60L                                    # main rolling window (trading days)
NP <- as.integer(Sys.getenv("DFZQ_PANEL_DAYS", "0")); if (is.na(NP)) NP <- 0L

## ================= 0. LOAD =================
hr("0. LOAD")
daily <- as.data.table(arrow::read_parquet("data/daily.parquet"))
daily[, date := as.Date(date)]
daily[, code := norm_code(code)]
all_days   <- sort(unique(daily$date))
panel_days <- if (NP > 0L) tail(all_days, NP) else all_days
P("[0] daily: %d rows | %d codes | panel %d days (%s ~ %s)%s",
  nrow(daily), uniqueN(daily$code), length(panel_days),
  min(panel_days), max(panel_days), if (NP > 0L) "  [SMOKE]" else "  [FULL]")

## ================= 1. PRICE / VOLUME FACTORS =================
hr("1. PRICE / VOLUME FACTORS")
dt <- copy(daily)[, .(code, date, open, high, low, close, volume, amount,
                      turnover, outstanding_share)]
setorder(dt, code, date)
dt[, ret := close / shift(close) - 1, by = code]
dt[, mkt_ret := mean(ret, na.rm = TRUE), by = date]      # equal-weight market proxy
P("[1] market proxy = cross-sectional mean of daily returns")

# --- counts first: every windowed factor requires a FULL window -----------------
dt[, nret := frollsum(as.integer(is.finite(ret)), W, align = "right"), by = code]

# --- rolling moments of returns (population form; n/(n-1) applied to vol only) --
dt[, rbar := frollmean(ret, W, align = "right", na.rm = TRUE), by = code]
dt[, r2   := frollmean(ret^2, W, align = "right", na.rm = TRUE), by = code]
dt[, r3   := frollmean(ret^3, W, align = "right", na.rm = TRUE), by = code]
dt[, rvar := r2 - rbar^2, by = code]
dt[, skew60 := fifelse(rvar > 0, (r3 - 3 * rbar * r2 + 2 * rbar^3) / rvar^1.5, NA_real_),
   by = code]
dt[, volatility := fifelse(nret > 1L,
                           sqrt(pmax(rvar, 0) * nret / (nret - 1L) * 252), NA_real_)]

# --- market moments: beta and idiosyncratic vol --------------------------------
dt[, mbar := frollmean(mkt_ret, W, align = "right", na.rm = TRUE), by = code]
dt[, m2   := frollmean(mkt_ret^2, W, align = "right", na.rm = TRUE), by = code]
dt[, rmix := frollmean(ret * mkt_ret, W, align = "right", na.rm = TRUE), by = code]
dt[, mvar := m2 - mbar^2, by = code]
dt[, beta60 := fifelse(mvar > 0, (rmix - rbar * mbar) / mvar, NA_real_), by = code]
dt[, ivol60 := fifelse(mvar > 0,
                       sqrt(pmax(rvar - beta60^2 * mvar, 0)) * sqrt(252), NA_real_),
   by = code]

# --- other price / volume constructs -------------------------------------------
# mcap = free-float market cap; px = VWAP price (see FIELD SEMANTICS above)
dt[, mcap := fifelse(turnover > 0, amount / turnover, NA_real_)]
dt[, px   := fifelse(volume   > 0, amount / volume,   NA_real_)]
dt[, reversal  := -(close / shift(close, 21L) - 1), by = code]
dt[, momentum  := shift(close, 21L) / shift(close, 271L) - 1, by = code]
dt[, maxret21  := frollmax(ret, 21L, align = "right", na.rm = TRUE), by = code]
dt[, amihud60  := frollmean(abs(ret) / pmax(amount, 1), W, align = "right", na.rm = TRUE),
   by = code]
dt[, hlrange60 := frollmean((high - low) / close, W, align = "right", na.rm = TRUE),
   by = code]
dt[, size      := log(pmax(mcap, 1e-6))]
dt[, liquidity := frollmean(turnover, W, align = "right", na.rm = TRUE), by = code]
dt[, turnchg   := frollmean(turnover, 5L, align = "right", na.rm = TRUE) /
                  frollmean(turnover, W,  align = "right", na.rm = TRUE) - 1, by = code]

# rolling corr(volume, ret) -- built from moments, no frollapply
dt[, vbar := frollmean(volume, W, align = "right", na.rm = TRUE), by = code]
dt[, v2   := frollmean(volume^2, W, align = "right", na.rm = TRUE), by = code]
dt[, vrm  := frollmean(volume * ret, W, align = "right", na.rm = TRUE), by = code]
dt[, vpcorr60 := fifelse(rvar > 0 & (v2 - vbar^2) > 0,
                         (vrm - vbar * rbar) / sqrt((v2 - vbar^2) * rvar), NA_real_),
   by = code]

FACS_PNEW <- c("reversal", "ivol60", "beta60", "skew60", "maxret21",
               "amihud60", "turnchg", "hlrange60", "vpcorr60")
FACS_POLD <- c("momentum", "volatility", "size", "liquidity")
# windowed factors need a full 60-day window (momentum needs 271 -> its own NA)
WINDED <- c("volatility", "ivol60", "beta60", "skew60", "amihud60",
            "hlrange60", "vpcorr60", "liquidity", "turnchg")
for (f in WINDED) dt[is.na(nret) | nret < W, (f) := NA_real_]

P("[1] price factors ready (%.0f s)", el())

pf   <- dt[date %in% panel_days, c("date", "code", FACS_POLD, FACS_PNEW), with = FALSE]
pxdt <- dt[date %in% panel_days, .(date, code, px)]     # price, for valuation ratios
mcdt <- dt[date %in% panel_days, .(date, code, mcap)]   # free-float mkt cap, for sp/ep_np
P("[1] price panel: %d rows x %d cols", nrow(pf), ncol(pf))

## ================= 2. FINANCIAL FACTORS (PIT) =================
hr("2. FINANCIAL FACTORS (point-in-time)")
fin <- as.data.table(arrow::read_parquet("data/financial.parquet"))
fin[, code := norm_code(code)]
IND <- c(roe = "净资产收益率(ROE)", gpm = "毛利率", npr = "销售净利率",
         bps = "每股净资产", eps = "基本每股收益", npg = "归属母公司净利润增长率",
         rg  = "营业总收入增长率", lev = "资产负债率", eq = "股东权益合计(净资产)",
         rev = "营业总收入", np = "归母净利润")
per_cols <- names(fin)[grepl("^\\d{8}$", names(fin))]
per_cols <- sort(per_cols[per_cols >= "20131231"])
P("[2] financial: %d rows | %d period cols | indicators matched: %d",
  nrow(fin), length(per_cols), sum(unname(IND) %in% unique(fin$指标)))

fl <- melt(fin[指标 %in% unname(IND)], id.vars = c("code", "指标"),
           measure.vars = per_cols, variable.name = "period", value.name = "value",
           variable.factor = FALSE)
setnames(fl, "指标", "ind")
fl[, key := names(IND)[match(ind, IND)]]
fl <- fl[!is.na(key)]
fl[, value := suppressWarnings(as.numeric(value))]
fl <- fl[is.finite(value)]
fl[, report_date := as.Date(period, format = "%Y%m%d")]
fl[, mm := as.integer(substr(period, 5L, 6L))]
# PIT: LEGAL disclosure deadline (annual +120d, H1 +62d, Q1 +30d, Q3 +31d)
fl[, avail_date := report_date + fcase(mm == 12L, 120L, mm == 6L, 62L,
                                       mm == 3L, 30L, default = 31L)]
fl <- unique(fl, by = c("code", "report_date", "key"))
fw <- dcast(fl, code + report_date + avail_date ~ key, value.var = "value")
P("[2] wide financial: %d rows x %d cols (%.0f s)", nrow(fw), ncol(fw), el())

# --- YoY growth of bps: shift by 4 quarters, but ONLY accept a true 1-year gap --
setorder(fw, code, report_date)
fw[, bps_yoy  := bps / shift(bps, 4L)  - 1, by = code]
fw[, gap4     := as.numeric(report_date - shift(report_date, 4L)), by = code]
fw[!(gap4 >= 355 & gap4 <= 375), bps_yoy := NA_real_]
P("[2] bps_yoy accepted: %d / %d rows | day-gap quartiles: %s",
  sum(is.finite(fw$bps_yoy)), nrow(fw),
  paste(as.integer(quantile(fw$gap4, c(.25, .5, .75), na.rm = TRUE)), collapse = "/"))

# --- roll forward: latest DISCLOSED statement as of each trading day -----------
fw[, c("gap4") := NULL]
setkey(fw, code, avail_date)
grid <- CJ(code = unique(fw$code), date = panel_days)
grid[, avail_date := date]
setkey(grid, code, avail_date)
matched <- fw[grid, roll = TRUE]
matched <- matched[!is.na(report_date)]
P("[2] matched: %d (date,code) rows | coverage vs panel = %.1f%%",
  nrow(matched), 100 * nrow(matched) / nrow(grid))

## unit sanity: rev/np/eq are in CNY and must be comparable to market cap (CNY)
chk <- merge(merge(matched, mcdt, by = c("date", "code")), pxdt,
             by = c("date", "code"))
P("[2] unit check (median over matched rows):")
P("     rev = %.3g   np = %.3g   eq = %.3g   mktcap = %.3g   px = %.3g",
  median(chk$rev, na.rm = TRUE), median(chk$np, na.rm = TRUE),
  median(chk$eq, na.rm = TRUE), median(chk$mcap, na.rm = TRUE),
  median(chk$px, na.rm = TRUE))
P("     rev/mktcap = %.3g (should be O(0.1-1)) | bps/px = %.3g",
  median(chk$rev / chk$mcap, na.rm = TRUE), median(chk$bps / chk$px, na.rm = TRUE))

## ================= 3. ASSEMBLE =================
hr("3. ASSEMBLE")
## price for valuation ratios = VWAP (amount/volume), NOT amount/turnover
m  <- merge(matched, pxdt, by = c("date", "code"))
m  <- merge(m, mcdt, by = c("date", "code"))
m[, `:=`(
  value_bp    = bps / px,                  # book-to-price  (per-share / PRICE)
  value_ep    = eps / px,                  # earnings yield  (per-share / PRICE)
  quality_roe = roe,
  quality_gpm = gpm,
  quality_npr = npr,
  growth_np   = npg,
  growth_rev  = rg,
  leverage    = -lev,
  sp          = rev / mcap,                # revenue / free-float mkt cap
  ep_np       = np  / mcap                 # net profit / free-float mkt cap
)]
FACS_FNEW <- c("sp", "ep_np", "bps_yoy")
FACS_FOLD <- c("value_bp", "value_ep", "quality_roe", "quality_gpm",
               "quality_npr", "growth_np", "growth_rev", "leverage")
FACS_ALL  <- c(FACS_POLD, FACS_PNEW, FACS_FOLD, FACS_FNEW)
P("[3] total factors: %d  (price %d + financial %d)",
  length(FACS_ALL), length(c(FACS_POLD, FACS_PNEW)), length(c(FACS_FOLD, FACS_FNEW)))

allf <- merge(m[, c("date", "code", FACS_FOLD, FACS_FNEW), with = FALSE],
              pf[, c("date", "code", FACS_POLD, FACS_PNEW), with = FALSE],
              by = c("date", "code"))
dup <- nrow(allf) - uniqueN(allf, by = c("date", "code"))
if (dup > 0L) P("[3] WARN: %d duplicate (date,code) rows removed", dup)
allf <- unique(allf, by = c("date", "code"))
P("[3] factor table: %d rows x %d cols | %d factors",
  nrow(allf), ncol(allf), length(FACS_ALL))

## winsorize(1%/99%) + cross-sectional z-score, by date
for (f in FACS_ALL) {
  allf[, (f) := {
    v <- get(f); q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
    v <- pmin(pmax(v, q[1]), q[2])
    s <- sd(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) (v - mean(v, na.rm = TRUE)) / s else NA_real_
  }, by = date]
}
setorder(allf, date, code)
arrow::write_parquet(allf, "_factor_panel2.parquet")
P("[3] -> _factor_panel2.parquet (%.1f MB)", file.size("_factor_panel2.parquet") / 1e6)

## ================= 4. DIAGNOSTICS =================
hr("4a. COVERAGE BY YEAR (%)")
allf[, yr := year(date)]
cov <- rbindlist(lapply(FACS_ALL, function(f)
  allf[, .(factor = f, cov = mean(is.finite(get(f)))), by = yr]))
print(dcast(cov, yr ~ factor, value.var = "cov")[
  , lapply(.SD, function(x) round(100 * x, 1)), by = yr])

hr("4b. EXPOSURE CORRELATION (pooled, complete rows) -- |rho| desc top 15")
M <- as.matrix(allf[, ..FACS_ALL])
keep <- rowSums(!is.finite(M)) == 0L
P("  rows used = %d / %d", sum(keep), nrow(M))
cm <- cor(M[keep, , drop = FALSE])
off <- cm[upper.tri(cm)]
P("  mean|rho| %.3f | max %.3f | pairs |rho|>0.5: %d / %d",
  mean(abs(off)), max(abs(off)), sum(abs(off) > 0.5), length(off))
pd <- data.table(a = rep(FACS_ALL, length(FACS_ALL))[upper.tri(cm)],
                 b = rep(FACS_ALL, each = length(FACS_ALL))[upper.tri(cm)], rho = off)
print(head(pd[order(-abs(rho))][, .(a, b, rho = round(rho, 3))], 15))

hr("4c. NEW FACTORS vs OLD ONES -- max |corr| against any old factor")
for (f in c(FACS_PNEW, FACS_FNEW)) {
  r <- abs(cm[f, setdiff(FACS_ALL, f)])
  P("  %-12s max|rho| = %.3f (with %s) | mean|rho| = %.3f",
    f, max(r), names(which.max(r)), mean(r))
}

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("NEXT: re-run _diag_dim.R against _factor_panel2.parquet (point it at the new")
P("      panel) and see whether the effective size (PR) and the block")
P("      structure become clearer than the current 3.22 / 12.")
