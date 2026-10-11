# =============================================================================
# _probe_units.R -- what ARE these fields, really?
# -----------------------------------------------------------------------------
# Question raised by _factor_ext2.R: median(outstanding_share * amount/turnover)
# came out at 2.0e18 CNY, which is absurd for an A-share market cap (should be
# ~1e10). Before trusting `mktcap` -- and before trusting the OLD size /
# value_bp / value_ep, which use the same implied price -- find out what the
# raw fields actually are.
#
# Usage (run from the project directory):
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_probe_units.R" `
#       | Tee-Object "_probe_units.log"
# Output: console only, no files written. ASCII output only.
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")

d <- as.data.table(arrow::read_parquet(
  "data/daily.parquet",
  col_select = c("date", "code", "close", "volume", "amount",
                 "turnover", "outstanding_share")))
d[, date := as.Date(date)][, code := norm_code(code)]

hr("SAMPLE ROWS -- 2025-06-03, a few well-known names")
s <- d[date == as.Date("2025-06-03") & code %in% c("600519", "000001", "600036",
                                                   "002594", "300750", "601398")]
print(s)

hr("POOLED MEDIANS (all rows, 2025)")
y <- d[year(date) == 2025]
print(sapply(y[, .(close, volume, amount, turnover, outstanding_share)],
             median, na.rm = TRUE))

hr("IMPLIED PRICE CANDIDATES (median over 2025, per code, then pooled)")
y[, `:=`(px_at = amount / turnover, px_av = amount / volume,
         px_vt = amount / (volume * turnover))]
print(sapply(y[, .(px_at, px_av, px_vt)], median, na.rm = TRUE))

hr("CORRELATION OF EACH CANDIDATE WITH qfq close (log scale, per-day cross-sec)")
cand <- c("px_at", "px_av", "px_vt")
for (k in cand) {
  cc <- y[is.finite(get(k)) & is.finite(close) & get(k) > 0 & close > 0,
          .(r = suppressWarnings(cor(log(get(k)), log(close)))), by = date]
  P("  log(%s) vs log(close) : mean cross-sec corr = %.4f", k, mean(cc$r))
}

hr("IMPLIED MARKET CAP CANDIDATES (median over 2025)")
P("  outstanding_share * close            = %.4g", median(y$outstanding_share * y$close, na.rm = TRUE))
P("  outstanding_share * amount/turnover  = %.4g", median(y$outstanding_share * y$px_at, na.rm = TRUE))
P("  outstanding_share * amount/volume    = %.4g", median(y$outstanding_share * y$px_av, na.rm = TRUE))
P("  outstanding_share (itself)           = %.4g", median(y$outstanding_share, na.rm = TRUE))
P("  volume (itself)                      = %.4g", median(y$volume, na.rm = TRUE))
P("  amount (itself)                      = %.4g", median(y$amount, na.rm = TRUE))
P("  turnover (itself)                    = %.4g", median(y$turnover, na.rm = TRUE))
P("")
P("  For reference, a median A-share: price ~10 CNY, shares ~5e8,")
P("  free-float mktcap ~1e10 CNY, daily amount ~1e8 CNY, turnover ~2%%.")
cat("\nDONE\n")
