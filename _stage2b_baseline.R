# =============================================================================
# _stage2b_baseline.R -- Step 2b: HONEST baseline table for the signal line
# -----------------------------------------------------------------------------
# THE QUESTION THIS ANSWERS:
#   Stage 1 showed that compressing 24 factors into ~10 blocks lifts the
#   VALIDATION ICIR from 0.498 (equal-weight all 24) to ~0.70. But the single
#   strongest factor on the TRAINING window already has |ICIR| = 0.798.
#   So: does the whole "cluster -> weight -> shrink" pipeline actually BEAT
#   "just use the single best factor"? Until that is answered, dimension
#   reduction has not earned its keep.
#
# DISCIPLINE (Sec 8.3 + the previous round's lessons):
#   - every weight is estimated on TRAIN ONLY and then FROZEN onto VALID
#   - the FROZEN OOS window is NOT touched (env DFZQ_TOUCH_OOS=1 to override,
#     which burns it -- do that at most once, at the very end)
#   - no "walking the validation set": G is a DESIGN choice justified by
#     split-half stability + interpretability, not by the validation argmax
#   - every variant is compared on the SAME rows and the SAME metric (Sec 9.4)
#
# VARIANTS
#   B0  single best factor, chosen on TRAIN by |ICIR|
#   B1  all 24 factors, equal weight
#   B2  G blocks (exposure-correlation Ward), equal weight across blocks
#   B3  B2 + Stage 2 GLS weights  w ~ Sigma_IC^-1 mu_IC
#   B4  B3 + Stage 3 shrinkage    w ~ (Sigma_IC + lambda I)^-1 mu_IC
#   B5  no clustering at all: 24 factors weighted by their TRAIN ICIR
#       (B5 isolates the contribution of the CLUSTERING step, i.e. B2/B3 vs B5)
#
# Usage (run from the project directory):
#   $env:DFZQ_PANEL="_factor_panel2.parquet"
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_stage2b_baseline.R" `
#       | Tee-Object "_stage2b_baseline.log"
# ASCII output only.
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

T0    <- Sys.time()
H     <- 21L
PANEL <- Sys.getenv("DFZQ_PANEL", "_factor_panel2.parquet")
START <- local({ s <- Sys.getenv("DFZQ_START", "")
                 if (nzchar(s)) as.Date(s) else as.Date("2017-01-01") })
TRAIN_END <- as.Date("2022-01-01")
VALID_END <- as.Date("2024-01-01")
TOUCH_OOS <- nzchar(Sys.getenv("DFZQ_TOUCH_OOS", ""))
G_DESIGN  <- 10L                     # design choice (see header)

## ================= 0. LOAD =================
hr("0. LOAD")
fp <- as.data.table(arrow::read_parquet(PANEL))
fp[, date := as.Date(date)]
FACS <- setdiff(names(fp), c("date", "code", "yr")); J <- length(FACS)

daily <- as.data.table(arrow::read_parquet(
  "data/daily.parquet", col_select = c("date", "code", "close")))
daily[, date := as.Date(date)][, code := norm_code(code)]
setorder(daily, code, date)
daily[, fwd := shift(close, -H) / close - 1, by = code]

m <- merge(fp[date >= START], daily[is.finite(fwd), .(date, code, fwd)],
           by = c("date", "code"))
rm(fp, daily); gc()
m[, win := fifelse(date < TRAIN_END, "TRAIN",
                   fifelse(date < VALID_END, "VALID", "OOS"))]
setorder(m, date, code)
P("[0] panel %s | %d rows | %d factors | %s ~ %s", PANEL, nrow(m), J,
  min(m$date), max(m$date))
P("[0] days per window: TRAIN=%d  VALID=%d  OOS=%d%s",
  uniqueN(m[win == "TRAIN", date]), uniqueN(m[win == "VALID", date]),
  uniqueN(m[win == "OOS", date]),
  if (TOUCH_OOS) "   [OOS WILL BE TOUCHED -- it is no longer frozen]" else "   [OOS frozen, not evaluated]")

## ================= 1. TRAIN IC -> SIGNS AND ICIR =================
hr("1. TRAIN IC (weights and signs come from here, TRAIN ONLY)")
ic <- rbindlist(lapply(FACS, function(f) {
  d <- m[is.finite(get(f)), .(win, date, z = get(f), fwd)]
  d[is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = .(win, date)][, factor := f][]
}))
ic <- ic[is.finite(ic)]
sm <- ic[, .(n = .N, ic_mean = mean(ic), icir = mean(ic) / sd(ic),
             t_hac = mean(ic) / nw_se(ic, 5L)), by = .(win, factor)]

icm_tr <- ic[win == "TRAIN", .(m = mean(ic)), by = factor]; setkey(icm_tr, factor)
icir_tr <- ic[win == "TRAIN", .(icir = mean(ic) / sd(ic)), by = factor]; setkey(icir_tr, factor)
sgn <- vapply(FACS, function(f) {
  s <- sign(icm_tr[f, m]); if (!is.finite(s) || s == 0) 1 else s }, numeric(1))
names(sgn) <- FACS
P("[1] TRAIN |ICIR| ranking (top 8):")
tmp <- data.table(factor = FACS, sgn = sgn,
                  icir = icir_tr[FACS, icir], ic_mean = icm_tr[FACS, m])
print(head(tmp[order(-abs(icir))][, .(factor, sgn, ic_mean = round(ic_mean, 4),
                                      icir = round(icir, 3))], 8))

## ================= 2. SCORE BUILDERS =================
hr("2. BUILDING SCORES (missing factor value -> 0, i.e. neutral)")
zcol <- function(f) { x <- m[[f]]; x[!is.finite(x)] <- 0; x }

# score from EFFECTIVE weights e (already includes each factor's train sign)
score_eff <- function(e) {
  v <- rep(0, nrow(m))
  for (j in seq_along(FACS)) if (e[j] != 0) v <- v + e[j] * zcol(FACS[j])
  v
}
# effective weights for "equal weight inside each block, w across blocks"
eff_block <- function(bm, w) {
  e <- numeric(J); names(e) <- FACS
  lev <- sort(unique(bm))
  for (g in seq_along(lev)) {
    fs <- FACS[bm == lev[g]]
    e[fs] <- w[g] * sgn[fs] / length(fs)
  }
  e
}
eval_score <- function(scv) {
  # NB: use set(), NOT `m[, sc := sc]`. data.table resolves the RHS of `:=`
  # with columns taking precedence, so once a column `sc` exists the assignment
  # becomes a no-op and every later call silently re-scores the FIRST vector.
  if (is.null(m[["sc_tmp"]])) m[, sc_tmp := NA_real_]
  set(m, j = "sc_tmp", value = as.numeric(scv))
  d <- m[is.finite(sc_tmp) & is.finite(fwd),
         .(ic = if (.N < 30L) NA_real_ else
             suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = .(win, date)]
  d <- d[is.finite(ic)]
  d[, .(ic_mean = mean(ic), icir = mean(ic) / sd(ic),
        t_hac = mean(ic) / nw_se(ic, 5L), pos = mean(ic > 0), days = .N), by = win]
}

## ================= 3. CLUSTERING (EXPOSURE CORR, TRAIN ONLY) =================
hr(sprintf("3. BLOCKS (exposure-correlation Ward, TRAIN only, G = %d by design)", G_DESIGN))
tr_days <- sort(unique(m[win == "TRAIN", date]))
sd_days <- tr_days[seq(1L, length(tr_days), by = 3L)]
CM <- cor(as.matrix(m[date %in% sd_days, ..FACS]), use = "pairwise.complete.obs")
dimnames(CM) <- list(FACS, FACS)
hc <- hclust(as.dist(1 - CM), method = "ward.D2")
bm <- cutree(hc, k = G_DESIGN)
for (k in sort(unique(bm))) P("  block %d: %s", k, paste(names(bm)[bm == k], collapse = ", "))

## ================= 4. BLOCK IC SERIES ON TRAIN -> mu, Sigma =================
hr("4. BLOCK IC ON TRAIN -> mu_IC and Sigma_IC (inputs to GLS)")
blk <- lapply(sort(unique(bm)), function(k) score_eff(eff_block(bm, { w <- numeric(G_DESIGN); w[k] <- 1; w })))
bi <- rbindlist(lapply(seq_along(blk), function(g) {
  m[, z := blk[[g]]]
  m[win == "TRAIN" & is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = date][is.finite(ic)][, block := g][]
}))
bw <- dcast(bi, date ~ block, value.var = "ic")
BICM <- as.matrix(bw[, -1]); colnames(BICM) <- paste0("b", colnames(BICM))
mu <- colMeans(BICM, na.rm = TRUE)
Sig <- cov(BICM, use = "pairwise.complete.obs")
P("[4] block IC means on TRAIN: %s", paste(sprintf("b%d=%.4f", seq_along(mu), mu), collapse = " "))
P("[4] cond(Sigma) = %.3g | T = %d days | G = %d", kappa(Sig), nrow(BICM), ncol(BICM))

## ================= 5. WEIGHTS =================
hr("5. WEIGHTS (all estimated on TRAIN, then frozen)")
w_eq <- rep(1 / G_DESIGN, G_DESIGN)
w_gls <- tryCatch(as.numeric(solve(Sig, mu)), error = function(e) NULL)
if (is.null(w_gls) || any(!is.finite(w_gls))) { P("[5] Sigma singular -> GLS fallback to equal"); w_gls <- w_eq }
w_gls <- w_gls / sum(abs(w_gls))
P("[5] equal-weight      : %s", paste(sprintf("%+.3f", w_eq), collapse = " "))
P("[5] GLS  w ~ S^-1 mu  : %s", paste(sprintf("%+.3f", w_gls), collapse = " "))
P("[5]   sum|w| = %.3f (normalized); negative/large weights are exactly what", sum(abs(w_gls)))
P("[5]   Stage 3 shrinkage is meant to tame.")

## --- Stage 3 shrinkage: pick lambda by K-fold CV INSIDE the training window ---
# NOTE: do NOT pick lambda by TRAIN ICIR. Shrinkage moves the weights AWAY from
# the in-sample optimum, so TRAIN ICIR is monotone DECREASING in lambda and that
# rule would always return lambda = 0 -- which is exactly what happened in the
# first full run (B3 and B4 came out numerically identical, i.e. Stage 3 was
# never actually tested). Proper CV: fit mu/Sigma on K-1 folds, score the
# held-out fold.
cv_folds  <- 5L
cv_dates  <- bw$date
fold_vec  <- rep(seq_len(cv_folds), length.out = length(cv_dates))
cv_icir <- function(lam) {
  v <- numeric(cv_folds)
  for (k in seq_len(cv_folds)) {
    Bk <- BICM[fold_vec != k, , drop = FALSE]
    mk <- colMeans(Bk, na.rm = TRUE)
    Sk <- cov(Bk, use = "pairwise.complete.obs")
    wk <- tryCatch(as.numeric(solve(Sk + lam * diag(G_DESIGN), mk)), error = function(e) NULL)
    if (is.null(wk) || any(!is.finite(wk))) { v[k] <- NA_real_; next }
    wk <- wk / sum(abs(wk))
    set(m, j = "sc_tmp", value = score_eff(eff_block(bm, wk)))
    d <- m[date %in% cv_dates[fold_vec == k] & is.finite(sc_tmp) & is.finite(fwd),
           .(ic = if (.N < 30L) NA_real_ else
               suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = date]
    d <- d[is.finite(ic)]
    v[k] <- if (nrow(d) > 5L) mean(d$ic) / sd(d$ic) else NA_real_
  }
  mean(v, na.rm = TRUE)
}
lam_grid <- c(0, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3) * mean(diag(Sig))
P("[5] Stage 3 lambda chosen by %d-fold CV INSIDE TRAIN (NOT by TRAIN ICIR):", cv_folds)
cvres <- data.table(lam = lam_grid, cv_icir = NA_real_)
for (k in seq_along(lam_grid)) {
  cvres$cv_icir[k] <- cv_icir(lam_grid[k])
  P("     lam=%10.6g -> CV ICIR %s", lam_grid[k],
    if (is.finite(cvres$cv_icir[k])) sprintf("%.3f", cvres$cv_icir[k]) else "NA")
}
kk  <- which.max(cvres$cv_icir)
LAM <- lam_grid[kk]
w_shr <- tryCatch(as.numeric(solve(Sig + LAM * diag(G_DESIGN), mu)), error = function(e) NULL)
if (is.null(w_shr) || any(!is.finite(w_shr))) w_shr <- w_gls
w_shr <- w_shr / sum(abs(w_shr))
P("[5] chosen lam = %.6g (CV ICIR %.3f) | w = %s", LAM, cvres$cv_icir[kk],
  paste(sprintf("%+.3f", w_shr), collapse = " "))

## ================= 6. THE TABLE =================
hr("6. BASELINE TABLE -- TRAIN (in-sample ref) and VALID (out-of-sample)")
b0 <- tmp[order(-abs(icir))][1]
variants <- list(
  list(id = "B0_single_best_factor", sc = score_eff({ e <- numeric(J); names(e) <- FACS
                                                       e[b0$factor] <- sgn[b0$factor]; e })),
  list(id = "B1_all24_equal",        sc = score_eff({ e <- rep(1 / J, J); names(e) <- FACS; e * sgn })),
  list(id = "B2_blocks_equal",       sc = score_eff(eff_block(bm, w_eq))),
  list(id = "B3_blocks_GLS",         sc = score_eff(eff_block(bm, w_gls))),
  list(id = "B4_blocks_GLS_shrunk",  sc = score_eff(eff_block(bm, w_shr))),
  # NB: icir ALREADY carries the direction -- do NOT multiply by sgn again
  # (that double-flip is what made the first run report B5 with ICIR -1.0).
  list(id = "B5_noCluster_ICIR_wt",  sc = score_eff({ e <- icir_tr[FACS, icir] / sum(abs(icir_tr[FACS, icir])); e }))
)
rows <- rbindlist(lapply(variants, function(v) {
  t <- eval_score(v$sc)
  keep <- if (TOUCH_OOS) c("TRAIN", "VALID", "OOS") else c("TRAIN", "VALID")
  t[win %in% keep][, variant := v$id][]
}))
setcolorder(rows, c("variant", "win", "ic_mean", "icir", "t_hac", "pos", "days"))
print(rows[, .(variant, win, ic_mean = round(ic_mean, 4), icir = round(icir, 3),
               t_hac = round(t_hac, 2), pos = round(pos, 3), days)])

hr("6b. HEADLINE -- VALIDATION ICIR")
print(rows[win == "VALID"][order(-icir), .(variant, icir = round(icir, 3), t_hac = round(t_hac, 2))])
best_v <- rows[win == "VALID"][which.max(icir)]
P("")
P("[6b] best on VALIDATION : %s (ICIR %.3f)", best_v$variant, best_v$icir)
P("[6b] B0 = the honest hurdle: a pipeline only earns its keep if it beats")
P("     B0 on VALIDATION, and the gap must survive the HAC t values above.")
P("[6b] NOTE: B0's factor was chosen on TRAIN, so this is a fair comparison.")

## ================= 6c. PAIRED TEST: is the gap over B0 real? =================
hr("6c. PAIRED TEST -- daily IC differences on VALIDATION")
ic_daily <- function(sc) {
  if (is.null(m[["sc_tmp"]])) m[, sc_tmp := NA_real_]
  set(m, j = "sc_tmp", value = as.numeric(sc))
  m[is.finite(sc_tmp) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))),
    by = .(win, date)][is.finite(ic)]
}
S <- lapply(variants, function(v) ic_daily(v$sc))
names(S) <- vapply(variants, function(v) v$id, character(1))
cmp <- list(c("B4_blocks_GLS_shrunk", "B0_single_best_factor"),
            c("B3_blocks_GLS",        "B0_single_best_factor"),
            c("B2_blocks_equal",      "B0_single_best_factor"),
            c("B5_noCluster_ICIR_wt", "B0_single_best_factor"),
            c("B4_blocks_GLS_shrunk", "B2_blocks_equal"))
pt <- rbindlist(lapply(cmp, function(pr) {
  a <- S[[pr[1]]][win == "VALID"][order(date)]
  b <- S[[pr[2]]][win == "VALID"][order(date)]
  d <- merge(a[, .(date, ia = ic)], b[, .(date, ib = ic)], by = "date")
  df <- d$ia - d$ib
  data.table(A = pr[1], B = pr[2], n = nrow(d), mean_diff = mean(df),
             t_hac = mean(df) / nw_se(df, 5L), win_rate = mean(df > 0))
}))
print(pt[, .(A, B, n, mean_diff = round(mean_diff, 5),
             t_hac = round(t_hac, 2), win_rate = round(win_rate, 3))])
P("[6c] |t_hac| > %.2f needed (Bonferroni over %d comparisons)",
   qnorm(1 - 0.05 / (2 * nrow(pt))), nrow(pt))
P("[6c] READ: a 5%% ICIR gap with |t| < 2 is NOT evidence the pipeline wins.")

## ================= 7. SENSITIVITY TO G (design choice, not tuned on VALID) ==
hr(sprintf("7. SENSITIVITY TO G -- does the conclusion depend on G = %d?", G_DESIGN))
grows <- rbindlist(lapply(c(6L, 8L, 10L, 12L, 15L), function(g) {
  bmg <- cutree(hc, k = g)
  bmg_list <- lapply(sort(unique(bmg)), function(k) score_eff(eff_block(bmg, { w <- numeric(g); w[k] <- 1; w })))
  big <- rbindlist(lapply(seq_along(bmg_list), function(i) {
    m[, z := bmg_list[[i]]]
    m[win == "TRAIN" & is.finite(z) & is.finite(fwd),
      .(ic = if (.N < 30L) NA_real_ else suppressWarnings(cor(z, fwd, method = "spearman"))),
      by = date][is.finite(ic)][, block := i][]
  }))
  bgw <- dcast(big, date ~ block, value.var = "ic")
  BM <- as.matrix(bgw[, -1]); mu_g <- colMeans(BM, na.rm = TRUE)
  Sg <- cov(BM, use = "pairwise.complete.obs")
  wg <- tryCatch(as.numeric(solve(Sg, mu_g)), error = function(e) NULL)
  if (is.null(wg) || any(!is.finite(wg))) wg <- rep(1 / g, g)
  wg <- wg / sum(abs(wg))
  ws <- tryCatch(as.numeric(solve(Sg + LAM * diag(g), mu_g)), error = function(e) NULL)
  if (is.null(ws) || any(!is.finite(ws))) ws <- wg
  ws <- ws / sum(abs(ws))
  tt <- eval_score(score_eff(eff_block(bmg, wg)))
  ts <- eval_score(score_eff(eff_block(bmg, ws)))
  data.table(G = g,
             gls_train = tt[win == "TRAIN", icir], gls_valid = tt[win == "VALID", icir],
             shr_train = ts[win == "TRAIN", icir], shr_valid = ts[win == "VALID", icir])
}))
print(grows[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("VERDICT LOGIC: if B3/B4 do not beat B0 on VALIDATION, dimension reduction")
P("has not earned its keep and should be stopped (same test as 优化点 B).")
