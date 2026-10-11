# =============================================================================
# _diag_dim.R -- Signal-line Step 0 diagnostics: is there room for dimension
#                reduction (factor clustering) in this data?
# -----------------------------------------------------------------------------
# Questions this script answers:
#   Q1 Redundancy : does the factor pool contain exploitable block structure?
#                   (no structure => clustering is a no-op, no room to reduce)
#   Q2 Stability  : is that structure stable over time?
#                   (unstable     => clustered weights are not trustworthy)
#   Q3 Info loss  : does squeezing J factors into G blocks lose signal?
#                   (composite-score IC vs best single factor)
#
# Acceptance-criteria source: SIMPLE_ALPHA_STRATEGY_OVERVIEW Copy.pdf
#   Sec 8.1 signal-layer metrics: Rank IC / ICIR / positive-IC ratio /
#           group monotonicity / coverage / factor correlation
#   Sec 5   "highly correlated, same-source factors must not be double
#           counted"  <- exactly what dimension reduction must deliver
#   Sec 8.3 train / validate / frozen-OOS time split
#
# Lessons carried over from the previous round (优化点 B):
#   - Do NOT conclude from a Marchenko-Pastur test: its verdict is driven
#     entirely by the N/T ratio (last round: 14 / 1 / 0 "residual factors"
#     on the same data depending on the setting). Use descriptive measures
#     instead: top-eigenvalue share / participation ratio / cum. variance.
#   - Always run a placebo: shuffle IC series in time; the correlation
#     matrix must collapse back to the identity.
#
# Usage (run from the project directory):
#   $env:DFZQ_PANEL_DAYS="400"        # smoke test, 1-2 min
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_diag_dim.R" `
#       | Tee-Object "_diag_dim.log"
# Console output is ASCII-only (Windows GBK console).
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")
el <- function() as.numeric(difftime(Sys.time(), T0, units = "secs"))
pm <- function(M, d = 2) print(round(M, d))

nw_se <- function(x, lag = 5L) {
  x <- x[is.finite(x)]; n <- length(x); if (n < 3L) return(NA_real_)
  e <- x - mean(x); s <- sum(e^2) / n
  for (l in seq_len(min(lag, n - 1L)))
    s <- s + 2 * (1 - l / (lag + 1)) * sum(e[(l + 1L):n] * e[1:(n - l)]) / n
  sqrt(s / n)
}

T0 <- Sys.time()
H  <- 21L
NP <- as.integer(Sys.getenv("DFZQ_PANEL_DAYS", "0")); if (is.na(NP)) NP <- 0L
SMOKE <- NP > 0L

# =============================================================================
# 0. DATA INVENTORY -- what does the "big environment" actually give us?
# =============================================================================
hr("0. DATA INVENTORY")
daily <- as.data.table(arrow::read_parquet(
  "data/daily.parquet", col_select = c("date", "code", "close", "amount", "turnover")))
daily[, date := as.Date(date)]
daily[, code := norm_code(code)]
all_days   <- sort(unique(daily$date))
panel_days <- if (NP > 0L) tail(all_days, NP) else all_days
P("[0] daily      : %d rows | %d codes | %s ~ %s",
  nrow(daily), uniqueN(daily$code), min(all_days), max(all_days))
P("[0] NPANEL     : %d days (%s ~ %s)%s", length(panel_days),
  min(panel_days), max(panel_days), if (SMOKE) "  [SMOKE]" else "  [FULL]")

dly <- daily[date %in% panel_days]
yr_cov <- dly[, .(n_codes = uniqueN(code)), by = .(yr = year(date))]
P("[0] stocks per year (panel):")
print(yr_cov)

PANEL <- Sys.getenv("DFZQ_PANEL", "_factor_panel.parquet")
fp <- as.data.table(arrow::read_parquet(PANEL))
fp[, date := as.Date(date)]
fp <- fp[date %in% panel_days]
FACS <- setdiff(names(fp), c("date", "code", "yr"))
J <- length(FACS)
P("[0] panel file: %s", PANEL)
P("[0] factor panel: %d rows | %d factors | %s ~ %s",
  nrow(fp), J, min(fp$date), max(fp$date))
P("[0] factors: %s", paste(FACS, collapse = ", "))

# =============================================================================
# 1. FORWARD RETURNS + DAILY CROSS-SECTIONAL IC  (Sec 8.1: Rank IC)
# =============================================================================
hr("1. DAILY CROSS-SECTIONAL RANK IC")
daily[, fwd := shift(close, -H) / close - 1, by = code]
fwd <- daily[is.finite(fwd), .(date, code, fwd)]
P("[1] forward %d-day returns: %d rows (%.0f s)", H, nrow(fwd), el())

m <- merge(fp, fwd, by = c("date", "code"))
P("[1] merged: %d rows (%.0f s)", nrow(m), el())

ic <- rbindlist(lapply(FACS, function(f) {
  d <- m[is.finite(get(f)), .(date, z = get(f), fwd)]
  d[is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = date][, factor := f][]
}))
ic <- ic[is.finite(ic)]
P("[1] IC series: %d rows | %d days (%.0f s)", nrow(ic), uniqueN(ic$date), el())

sm <- ic[, {
  se <- nw_se(ic, 5L)
  .(n = .N, ic_mean = mean(ic), ic_sd = sd(ic),
    icir = mean(ic) / sd(ic), t_hac = mean(ic) / se,
    pos = mean(ic > 0))
}, by = factor]
sm[, abs_t := abs(t_hac)]
setorder(sm, -abs_t)
P("[1] IC summary  (|t|>2 => detectable predictive power; sign matters):")
print(sm[, .(factor, n, ic_mean = round(ic_mean, 4), icir = round(icir, 3),
             t_hac = round(t_hac, 2), pos = round(pos, 3))])
P("  factors with POSITIVE ic_mean: %d / %d", sum(sm$ic_mean > 0), J)
P("  factors with |t| > 2         : %d / %d", sum(sm$abs_t > 2), J)

# =============================================================================
# 2. Q1 -- REDUNDANCY: IC correlation matrix, spectrum, effective size
# =============================================================================
hr("2. Q1 REDUNDANCY -- IC CORRELATION STRUCTURE")
iw <- dcast(ic, date ~ factor, value.var = "ic")
CM <- cor(as.matrix(iw[, -1]), use = "pairwise.complete.obs")
P("[2] IC correlation matrix (J = %d, T = %d days):", J, nrow(iw))
dimnames(CM) <- list(FACS, FACS)
pm(CM, 2)

off <- CM[upper.tri(CM)]
P("\n[2] off-diagonal |rho| distribution:")
P("    mean %.3f | median %.3f | max %.3f | p90 %.3f",
  mean(abs(off)), median(abs(off)), max(abs(off)),
  quantile(abs(off), 0.9, names = FALSE))
P("    pairs with |rho|>0.3 : %d / %d (%.0f%%)",
  sum(abs(off) > 0.3), length(off), 100 * mean(abs(off) > 0.3))
P("    pairs with |rho|>0.5 : %d / %d (%.0f%%)",
  sum(abs(off) > 0.5), length(off), 100 * mean(abs(off) > 0.5))
P("    pairs with |rho|>0.7 : %d / %d (%.0f%%)",
  sum(abs(off) > 0.7), length(off), 100 * mean(abs(off) > 0.7))
P("    most-correlated pairs (|rho| desc, top 10):")
pairs_dt <- data.table(
  a = rep(FACS, J)[upper.tri(CM)],
  b = rep(FACS, each = J)[upper.tri(CM)],
  rho = off)
print(head(pairs_dt[order(-abs(rho))][, .(a, b, rho = round(rho, 3))], 10))

ev <- eigen(CM, symmetric = TRUE)$values
Jf <- sum(ev > 1e-8)
share1   <- ev[1] / sum(ev)                 # top-eigenvalue share
share10  <- sum(ev[seq_len(min(10, Jf))]) / sum(ev)
pr       <- sum(ev)^2 / sum(ev^2)           # participation ratio (effective #)
cumev    <- cumsum(ev) / sum(ev)
npc90    <- which(cumev >= 0.90)[1]
npc95    <- which(cumev >= 0.95)[1]
P("\n[2] spectrum of the IC correlation matrix (sum(lambda) = J = %d):", J)
print(data.table(k = seq_along(ev), lambda = round(ev, 3),
                 share = round(ev / sum(ev), 3), cum = round(cumev, 3)))
P("\n[2] effective-size measures (descriptive, NOT an MP test):")
P("    lambda1 share        = %.3f   (1/J = %.3f if NO structure)", share1, 1 / J)
P("    top-10 share         = %.3f", share10)
P("    participation ratio  = %.2f   (J = %d; near J => no structure)", pr, J)
P("    Kaiser (# lambda > 1)= %d / %d", sum(ev > 1), J)
P("    PCs for 90%% / 95%%     = %d / %d", npc90, npc95)
P("    min lambda           = %.3e  (negative => non-PSD pairwise matrix)", min(ev))

# loading of the top eigenvector: which factors form the "common mode"
ld <- data.table(factor = FACS, load1 = eigen(CM, symmetric = TRUE)$vectors[, 1])
ld <- ld[order(-abs(load1))]
P("\n[2] loadings on the 1st eigenvector (|.| desc) -- the dominant common mode:")
print(head(ld[, .(factor, load1 = round(load1, 3))], 12))

# =============================================================================
# 3. Q2 -- STABILITY: split-half and year-by-year
# =============================================================================
hr("3. Q2 STABILITY -- split-half + by-year")
dd   <- sort(unique(iw$date)); mid <- dd[floor(length(dd) / 2)]
h1   <- iw[date <= mid]; h2 <- iw[date >  mid]
CM1  <- cor(as.matrix(h1[, -1]), use = "pairwise.complete.obs")
CM2  <- cor(as.matrix(h2[, -1]), use = "pairwise.complete.obs")
tr   <- upper.tri(CM)
P("[3] split at %s : %d days / %d days", mid, nrow(h1), nrow(h2))
P("    corr(structure H1, structure H2) on %d off-diagonal pairs = %.3f",
  sum(tr), cor(CM1[tr], CM2[tr]))
s1 <- ic[date <= mid, .(s1 = mean(ic)), by = factor]
s2 <- ic[date >  mid, .(s2 = mean(ic)), by = factor]
sg <- merge(s1, s2, by = "factor")
P("    corr(IC mean H1, IC mean H2) across %d factors = %.3f", nrow(sg),
  cor(sg$s1, sg$s2))
P("    factors whose IC SIGN flips between halves: %d / %d",
  sum(sign(sg$s1) != sign(sg$s2)), nrow(sg))
print(sg[, .(factor, H1 = round(s1, 4), H2 = round(s2, 4),
             flip = sign(s1) != sign(s2))][order(H1)])

icy <- ic[, .(ic = mean(ic)), by = .(factor, yr = year(date))]
wicy <- dcast(icy, factor ~ yr, value.var = "ic")
P("\n[3] IC mean by year (sign flips across years => unstable signal):")
ycols <- setdiff(names(wicy), "factor")
print(cbind(data.table(factor = wicy$factor),
            wicy[, lapply(.SD, function(x) round(x, 3)), .SDcols = ycols]))

# =============================================================================
# 4. Q1b -- PLACEBO: shuffle IC in time; structure must disappear
# =============================================================================
hr("4. PLACEBO -- time-shuffled IC must show NO structure")
M0 <- as.matrix(iw[, -1])
set.seed(20261010L)
Mp <- M0
for (j in seq_len(ncol(Mp))) Mp[, j] <- Mp[sample.int(nrow(Mp)), j]
CMp <- cor(Mp, use = "pairwise.complete.obs")
offp <- CMp[upper.tri(CMp)]
evp  <- eigen(CMp, symmetric = TRUE)$values
P("[4] observed : mean|rho| %.3f | lambda1 share %.3f | PR %.2f",
  mean(abs(off)), share1, pr)
P("[4] placebo  : mean|rho| %.3f | lambda1 share %.3f | PR %.2f",
  mean(abs(offp)), evp[1] / sum(evp), sum(evp)^2 / sum(evp^2))
P("    READ: if placebo PR is also far below J, the 'structure' is an artifact")
P("          of finite-sample noise (J small, T large is NOT innocent).")

# =============================================================================
# 4b. EFFECTIVE SAMPLE -- overlapping forward windows inflate everything
# -----------------------------------------------------------------------------
# IC_f(t) uses the return over [t, t+H], so consecutive ICs share H-1 of H days.
# The IC series is therefore a heavily smoothed object: the effective number of
# INDEPENDENT observations is ~T/H, not T. Two smoothed series can correlate
# strongly even when the underlying daily quantities are unrelated. The
# time-shuffle placebo does NOT catch this (shuffling destroys the smoothing
# too). The only honest check is to re-estimate on non-overlapping blocks.
# =============================================================================
hr(sprintf("4b. EFFECTIVE SAMPLE -- overlapping %d-day forward returns", H))
acf1 <- mean(sapply(FACS, function(f) { x <- iw[[f]]; x <- x[is.finite(x)]
  suppressWarnings(cor(x[-1], x[-length(x)])) }))
acfH <- mean(sapply(FACS, function(f) { x <- iw[[f]]; x <- x[is.finite(x)]
  n <- length(x) - H
  suppressWarnings(cor(x[seq_len(n)], x[seq_len(n) + H])) }))
P("[4b] mean lag-1  autocorr of IC series = %.3f", acf1)
P("[4b] mean lag-%d autocorr of IC series = %.3f (high => smooth series)", H, acfH)
P("[4b] effective independent blocks ~ T/H = %d / %d = %.1f", nrow(iw), H, nrow(iw) / H)

sub_idx <- seq(1L, nrow(iw), by = H)
CMs <- cor(as.matrix(iw[sub_idx, -1]), use = "pairwise.complete.obs")
offs <- CMs[upper.tri(CMs)]
evs  <- eigen(CMs, symmetric = TRUE)$values
P("[4b] NON-OVERLAPPING subsample: %d days", length(sub_idx))
P("     mean|rho| %.3f | pairs|rho|>0.5 : %.0f%% | lambda1 share %.3f | PR %.2f",
  mean(abs(offs)), 100 * mean(abs(offs) > 0.5),
  evs[1] / sum(evs), sum(evs)^2 / sum(evs^2))
P("     observed (overlapping, sec 2): mean|rho| %.3f | PR %.2f", mean(abs(off)), pr)
P("     READ: if this collapses toward the placebo level (~0.035), the")
P("           structure in sec 2 was an artifact of overlapping windows.")

# =============================================================================
# 5. EXPOSURE CORRELATION (control) -- the plan's alternative similarity
# =============================================================================
hr("5. CONTROL -- exposure correlation (not IC correlation)")
samp <- fp[date %in% dd[seq(1, length(dd), by = max(1L, length(dd) %/% 40L))]]
Xs <- as.matrix(samp[, ..FACS])
keep <- rowSums(!is.finite(Xs)) == 0L
CX <- cor(Xs[keep, , drop = FALSE])
offx <- CX[upper.tri(CX)]
P("[5] rows used: %d / %d", sum(keep), nrow(Xs))
P("    mean|rho| %.3f | max %.3f | pairs |rho|>0.5 : %d / %d",
  mean(abs(offx)), max(abs(offx)), sum(abs(offx) > 0.5), length(offx))
P("    corr(IC-corr, exposure-corr) on shared pairs = %.3f", {
  ok <- is.finite(CM[tr]) & is.finite(CX[tr]); cor(CM[tr][ok], CX[tr][ok])
})
px_dt <- data.table(a = rep(FACS, J)[upper.tri(CX)],
                    b = rep(FACS, each = J)[upper.tri(CX)], rho = offx)
P("    top exposure pairs (|rho| desc, top 10):")
print(head(px_dt[order(-abs(rho))][, .(a, b, rho = round(rho, 3))], 10))

# --- 5b. SANITY: are the top IC-correlated pairs REALLY the same factor? -----
# Two factors can correlate at the IC level only if they rank stocks similarly
# each day. Check the pair-level daily cross-sectional correlation directly.
# Near 1.0 between economically unrelated factors => construction bug upstream
# (in _factor_ext.R), NOT genuine redundancy.
P("")
P("[5b] SANITY CHECK on the top IC-correlated pairs -- mean daily")
P("     cross-sectional correlation of the raw z-scored exposures:")
chk <- head(pairs_dt[order(-abs(rho))], 6)
sub <- fp[date %in% dd]
for (i in seq_len(nrow(chk))) {
  a <- chk$a[i]; b <- chk$b[i]
  v <- sub[is.finite(get(a)) & is.finite(get(b)),
           .(r = suppressWarnings(cor(get(a), get(b)))), by = date]
  P("     %-12s vs %-12s : IC-corr %+.3f | exposure-corr %+.3f (n=%d days)",
    a, b, chk$rho[i], mean(v$r), nrow(v))
}
P("     READ: exposure-corr ~= IC-corr => the two columns carry the same")
P("           cross-sectional information (~duplicate). If the pair is")
P("           economically unrelated, suspect the factor construction.")

# =============================================================================
# 6. Q3 -- INFORMATION LOSS: does combining many factors beat one?
# =============================================================================
hr("6. Q3 INFORMATION LOSS -- composite IC vs best single factor")
sgn <- setNames(sign(sm$ic_mean), sm$factor)
sgn[!is.finite(sgn) | sgn == 0] <- 1
for (f in FACS) m[, (paste0("s_", f)) := get(f) * sgn[[f]]]
scols <- paste0("s_", FACS)
m[, sc_ew := rowMeans(.SD, na.rm = TRUE), .SDcols = scols]

compose <- function(cols, nm) {
  if (length(cols) == 1L) m[, sc := get(cols)]
  else m[, sc := rowMeans(.SD, na.rm = TRUE), .SDcols = cols]
  d <- m[, .(date, z = sc, fwd)][is.finite(z) & is.finite(fwd)]
  s <- d[, .(ic = if (.N < 30L) NA_real_ else
             suppressWarnings(cor(z, fwd, method = "spearman"))), by = date][is.finite(ic)]
  se <- nw_se(s$ic, 5L)
  data.table(variant = nm, nfac = length(cols), days = nrow(s),
             ic_mean = mean(s$ic), icir = mean(s$ic) / sd(s$ic),
             t_hac = mean(s$ic) / se, pos = mean(s$ic > 0))
}
best1 <- sm$factor[1]
top3  <- paste0("s_", sm$factor[seq_len(min(3, J))])   # sign-adjusted columns
res <- rbindlist(list(
  compose(scols, "EW_all_sign_adj"),
  compose(top3, "EW_top3_by_|t|_sign_adj"),
  compose(paste0("s_", best1), paste0("best1:", best1))
))
P("[6] composite-score IC (sign-adjusted equal weight; full-sample signs = oracle,")
P("    so read this as an UPPER BOUND, not an OOS result):")
print(res[, .(variant, nfac, days, ic_mean = round(ic_mean, 4),
              icir = round(icir, 3), t_hac = round(t_hac, 2), pos = round(pos, 3))])

# =============================================================================
# 7. VERDICT
# =============================================================================
hr("7. VERDICT -- is there room for dimension reduction?")
P("  Q1 redundancy : mean|rho| %.3f | pairs|rho|>0.5 = %.0f%% | PR %.2f / J = %d",
  mean(abs(off)), 100 * mean(abs(off) > 0.5), pr, J)
P("  Q2 stability  : corr(H1,H2) structure = %.3f | IC-mean corr = %.3f",
  cor(CM1[tr], CM2[tr]), cor(sg$s1, sg$s2))
P("  Q3 info loss  : EW_all ICIR = %.3f vs best single ICIR = %.3f",
  res$icir[1], res$icir[3])
P("")
P("  Rule of thumb (descriptive, no hard threshold):")
P("    - PR close to J  AND few |rho|>0.5 pairs  => little redundancy;")
P("      clustering would be cosmetic. Widen the pool BEFORE reducing.")
P("    - clear |rho|>0.5 clusters that are stable across halves => real room.")
P("    - EW_all ICIR <= best single ICIR => the extra factors add little,")
P("      so dimension reduction cannot LOSE much either (and cannot gain much).")
P("")
P("  NOTE: signs above are estimated on the full sample -> oracle. Stage 1/2")
P("        must re-estimate every weight on the TRAINING window only (Sec 8.3).")

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("NEXT: if redundancy is thin, extend the factor pool (widely used factors)")
P("      and re-run this script; then decide G* via BIC + F test in Stage 1.")
