# =============================================================================
# _stage4_equiv.R -- is the compression loss NEGLIGIBLE? (two tests, not one)
# -----------------------------------------------------------------------------
# "Not much different" is an EQUIVALENCE claim, not a "failed to reject" claim.
# This script therefore reports TWO things, because they answer different
# questions and will disagree:
#
# (1) NESTED-MODEL F TEST -- the ETX5460 Project 2 / HAR convention
#       F = ((RSS_r - RSS_u)/q) / (RSS_u/(N - k_u))
#     Compression imposes "coefficients are equal within a block", i.e.
#       gamma_j = beta_g / |block g|,  q = J - G restrictions.
#     (HAR in AR(10): q = 10 - 3 = 7.  Here: q = 24 - G.)
#     The model is a panel regression of forward return on the factor
#     exposures with DAY FIXED EFFECTS (so slopes are common across days,
#     exactly parallel to HAR being one time-series regression).
#     *** WARNING *** with N ~ 4.2e6 the residual df is ~4.2e6, so this test
#     has enormous power: it will reject "the restriction holds exactly" for
#     ANY nonzero loss. It answers "is the restriction literally true", NOT
#     "is the loss ignorable".
#
# (2) EQUIVALENCE (TOST) + 95% CONFIDENCE INTERVAL -- answers "can we SAY it
#     is not much different". For the daily cross-sectional IC difference and
#     for the daily rank correlation between the compressed and uncompressed
#     scores, we build a HAC (Newey-West) 95% CI and check whether it lies
#     entirely inside a pre-declared negligible band (-delta, +delta).
#     A HAC CI is required: the IC series is autocorrelated (overlapping
#     21-day forward returns, lag-1 autocorr ~ 0.96).
#
# delta is reported at 5% / 10% / 20% of the reference IC so the choice of
# "ignorable" is explicit and not smuggled in.
#
# Usage:  $env:DFZQ_PANEL="_factor_panel2.parquet"
#         & "...\Rscript.exe" "_stage4_equiv.R" | Tee-Object "_stage4_equiv.log"
# ASCII output only. OOS stays frozen.
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")
el <- function() as.numeric(difftime(Sys.time(), T0, units = "secs"))

nw_se <- function(x, lag = 5L) {
  x <- x[is.finite(x)]; n <- length(x); if (n < 3L) return(NA_real_)
  e <- x - mean(x); s <- sum(e^2) / n
  for (l in seq_len(min(lag, n - 1L)))
    s <- s + 2 * (1 - l / (lag + 1)) * sum(e[(l + 1L):n] * e[1:(n - l)]) / n
  sqrt(s / n)
}

T0 <- Sys.time(); H <- 21L
PANEL <- Sys.getenv("DFZQ_PANEL", "_factor_panel2.parquet")
START <- local({ s <- Sys.getenv("DFZQ_START", "")
                 if (nzchar(s)) as.Date(s) else as.Date("2017-01-01") })
TRAIN_END <- as.Date("2022-01-01"); VALID_END <- as.Date("2024-01-01")

## ================= 0. LOAD =================
hr("0. LOAD")
fp <- as.data.table(arrow::read_parquet(PANEL)); fp[, date := as.Date(date)]
FACS <- setdiff(names(fp), c("date", "code", "yr")); J <- length(FACS)
daily <- as.data.table(arrow::read_parquet(
  "data/daily.parquet", col_select = c("date", "code", "close")))
daily[, date := as.Date(date)][, code := norm_code(code)]
setorder(daily, code, date)
daily[, fwd := shift(close, -H) / close - 1, by = code]
m <- merge(fp[date >= START], daily[is.finite(fwd), .(date, code, fwd)], by = c("date", "code"))
rm(fp, daily); gc()
m[, win := fifelse(date < TRAIN_END, "TRAIN", fifelse(date < VALID_END, "VALID", "OOS"))]
setorder(m, date, code)
P("[0] %d rows | %d factors | TRAIN/VALID/OOS = %d/%d/%d days (OOS frozen)",
  nrow(m), J, uniqueN(m[win == "TRAIN", date]), uniqueN(m[win == "VALID", date]),
  uniqueN(m[win == "OOS", date]))

## ================= 1. SIGNS FROM TRAIN =================
hr("1. SIGNS FROM TRAIN")
ic <- rbindlist(lapply(FACS, function(f) {
  d <- m[is.finite(get(f)), .(win, date, z = get(f), fwd)]
  d[is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = .(win, date)][, factor := f][]
}))
ic <- ic[is.finite(ic)]
icm <- ic[win == "TRAIN", .(m = mean(ic)), by = factor]; setkey(icm, factor)
sgn <- vapply(FACS, function(f) { s <- sign(icm[f, m]); if (!is.finite(s) || s == 0) 1 else s }, numeric(1))
names(sgn) <- FACS

zcol <- function(f) { x <- m[[f]]; x[!is.finite(x)] <- 0; x }
setcol <- function(nm, v) { if (is.null(m[[nm]])) m[, (nm) := NA_real_]; set(m, j = nm, value = as.numeric(v)) }
score_eff <- function(e) { v <- rep(0, nrow(m))
  for (j in seq_along(FACS)) if (e[j] != 0) v <- v + e[j] * zcol(FACS[j]); v }
eff_block <- function(bm, w) { e <- numeric(J); names(e) <- FACS; lev <- sort(unique(bm))
  for (g in seq_along(lev)) { fs <- FACS[bm == lev[g]]; e[fs] <- w[g] * sgn[fs] / length(fs) }; e }

## ================= 2. BLOCKS =================
hr("2. BLOCKS (exposure-correlation Ward, TRAIN only)")
tr_days <- sort(unique(m[win == "TRAIN", date]))
sd_days <- tr_days[seq(1L, length(tr_days), by = 3L)]
CM <- cor(as.matrix(m[date %in% sd_days, ..FACS]), use = "pairwise.complete.obs")
dimnames(CM) <- list(FACS, FACS)
hc <- hclust(as.dist(1 - CM), method = "ward.D2")

## ================= 3. NESTED-MODEL F TEST =================
hr("3. NESTED-MODEL F TEST (HAR / ETX5460 P2 convention)")
P("  model: fwd ~ factor exposures, day fixed effects (slopes common across days)")
P("  unrestricted: J = %d slopes | restricted (G blocks): G slopes | q = J - G", J)

# per-day sufficient statistics on DAY-DEMEANED data (day FE is removed by
# within-day demeaning; an intercept is then unnecessary)
setkey(m, date)
st <- list(); Ntot <- 0L; dropped <- 0L
for (i in seq_along(tr_days)) {
  sub <- m[.(tr_days[i])]
  X <- as.matrix(sub[, ..FACS]); y <- sub$fwd
  ok <- is.finite(y); X <- X[ok, , drop = FALSE]; y <- y[ok]
  X[!is.finite(X)] <- 0
  y <- y - mean(y)                       # day FE
  X <- sweep(X, 2L, colMeans(X))         # day FE on regressors
  # A day whose design is rank-deficient (e.g. a factor with zero coverage that
  # day, so its column is constant) cannot support the unrestricted model.
  # Skip it and report the count -- do NOT silently ridge it away.
  if (qr(X)$rank < ncol(X)) { dropped <- dropped + 1L; next }
  st[[length(st) + 1L]] <- list(G = crossprod(X), b = crossprod(X, y),
                                yy = sum(y^2), n = length(y))
  Ntot <- Ntot + length(y)
}
n_days <- length(st)
P("  days used = %d | days SKIPPED (rank-deficient design) = %d", n_days, dropped)
P("  N = %d pooled stock-days | days = %d | residual df (unrestricted) = %d",
  Ntot, n_days, Ntot - n_days - J)

rss_un <- sum(vapply(st, function(s) s$yy - sum(s$b * solve(s$G, s$b)), numeric(1)))
rss_of <- function(A) sum(vapply(st, function(s) {
  Ag <- crossprod(A, s$G %*% A); bb <- crossprod(A, s$b)
  sol <- tryCatch(solve(Ag, bb), error = function(e) NULL)
  if (is.null(sol) || any(!is.finite(sol))) return(NA_real_)
  s$yy - sum(bb * sol) }, numeric(1)))
ss_tot <- sum(vapply(st, function(s) s$yy, numeric(1)))
P("  RSS_unrestricted (J = %d) = %.6g | within-day total SS = %.6g", J, rss_un, ss_tot)
P("  R2 (unrestricted) = %.4f%%  <-- so a 3%% RSS rise is a 3%% relative loss",
  100 * (1 - rss_un / ss_tot))
P("                                   of an already tiny explained share")

ftab <- rbindlist(lapply(c(15L, 12L, 10L, 8L, 6L, 4L, 3L, 2L), function(G) {
  bm <- cutree(hc, k = G); lev <- sort(unique(bm))
  A <- matrix(0, J, length(lev), dimnames = list(FACS, NULL))
  for (g in seq_along(lev)) A[bm == lev[g], g] <- 1 / sum(bm == lev[g])
  rss_r <- rss_of(A)
  q <- J - G; df2 <- Ntot - n_days - J
  Fst <- ((rss_r - rss_un) / q) / (rss_un / df2)
  data.table(G = G, q = q, RSS_r = rss_r, F = Fst, p = pf(Fst, q, df2, lower.tail = FALSE),
             rel_loss = (rss_r - rss_un) / rss_un)
}))
print(ftab[, .(G, q, RSS_r = signif(RSS_r, 7), F = signif(F, 4),
               p = format.pval(p, digits = 3, eps = 1e-300),
               rel_loss_pct = round(100 * rel_loss, 4))])
P("")
P("  READ: residual df ~ %d => this test has enormous power. It asks whether",
  Ntot - n_days - J)
P("        the block restriction holds EXACTLY, not whether the loss matters.")
P("        A large F with rel_loss_pct ~ 0.01%% still means the loss is tiny.")
P("        For contrast, the HAR test had T-k = 1022 and F = 0.643, p = 0.720.")

## ================= 4. EQUIVALENCE (TOST) =================
hr("4. EQUIVALENCE TEST -- can we SAY the output is 'not much different'?")
ref_sc <- score_eff(rep(1 / J, J) * sgn)
setcol("sc_ref", ref_sc)
setcol("sc_tmp", ref_sc)
ref_ic_va <- m[win == "VALID" & is.finite(sc_tmp) & is.finite(fwd),
               .(ic = suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = date][order(date)]
ref_mean <- mean(ref_ic_va$ic)
P("  reference (uncompressed) VALID mean IC = %.4f | days = %d", ref_mean, nrow(ref_ic_va))

etab <- rbindlist(lapply(c(15L, 12L, 10L, 8L, 6L, 4L, 3L, 2L), function(G) {
  bm <- cutree(hc, k = G); lev <- sort(unique(bm))
  s_g <- score_eff(eff_block(bm, rep(1 / length(lev), length(lev))))
  setcol("sc_tmp", s_g)
  # (a) daily rank correlation between the compressed and uncompressed scores
  d1 <- m[win == "VALID" & is.finite(sc_tmp) & is.finite(sc_ref),
          .(r = suppressWarnings(cor(sc_tmp, sc_ref, method = "spearman"))), by = date]
  d1 <- d1[is.finite(r)]
  # (b) daily IC of the compressed score, and the paired difference
  iv <- m[win == "VALID" & is.finite(sc_tmp) & is.finite(fwd),
          .(ic = suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = date][order(date)]
  dd <- merge(iv, ref_ic_va, by = "date", suffixes = c(".g", ".ref"))
  df <- dd$ic.g - dd$ic.ref
  se_df <- nw_se(df, 5L); se_rc <- nw_se(d1$r, 5L)
  data.table(G = G,
             corr_mean = mean(d1$r), corr_lo = mean(d1$r) - 1.96 * se_rc, corr_hi = mean(d1$r) + 1.96 * se_rc,
             dIC_mean = mean(df), dIC_lo = mean(df) - 1.96 * se_df, dIC_hi = mean(df) + 1.96 * se_df,
             dIC_pct = 100 * mean(df) / ref_mean)
}))
print(etab[, .(G, corr_mean = round(corr_mean, 3),
               corr_95CI = sprintf("[%.3f, %.3f]", corr_lo, corr_hi),
               dIC_mean = round(dIC_mean, 5),
               dIC_95CI = sprintf("[%.5f, %.5f]", dIC_lo, dIC_hi),
               dIC_pct_of_ref = round(dIC_pct, 2))])

hr("4b. TOST DECISION at delta = 5% / 10% / 20% of the reference IC")
P("  (equivalence is declared only if the WHOLE 95%% CI lies inside (-delta, +delta))")
P("  delta(5%%) = %.5f | delta(10%%) = %.5f | delta(20%%) = %.5f",
  0.05 * abs(ref_mean), 0.10 * abs(ref_mean), 0.20 * abs(ref_mean))
for (k in 1:nrow(etab)) {
  P("   G=%2d : dIC = %+.5f  CI [%+.5f, %+.5f]  ->  EQUIV 5%%:%s 10%%:%s 20%%:%s | NI(10%%) non-inferior:%s",
    etab$G[k], etab$dIC_mean[k], etab$dIC_lo[k], etab$dIC_hi[k],
    ifelse(etab$dIC_lo[k] > -0.05 * abs(ref_mean) & etab$dIC_hi[k] < 0.05 * abs(ref_mean), "EQUIV", "no"),
    ifelse(etab$dIC_lo[k] > -0.10 * abs(ref_mean) & etab$dIC_hi[k] < 0.10 * abs(ref_mean), "EQUIV", "no"),
    ifelse(etab$dIC_lo[k] > -0.20 * abs(ref_mean) & etab$dIC_hi[k] < 0.20 * abs(ref_mean), "EQUIV", "no"),
    ifelse(etab$dIC_lo[k] > -0.10 * abs(ref_mean), "YES", "no"))
}
P("")
P("  NOTE: two-sided EQUIVENCE fails when compression is significantly BETTER.")
P("        For a 'no loss' claim the ONE-SIDED non-inferiority test (NI) is the")
P("        right tool: it only asks whether the LOWER bound is worse than -delta.")

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("LOGIC: F test = 'is the restriction exactly true' (expect REJECT).")
P("       TOST   = 'can we SAY the loss is ignorable' (this is the claim asked for).")
