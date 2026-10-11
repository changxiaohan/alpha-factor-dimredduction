# =============================================================================
# _stage3_compress.R -- the COMPRESSION CURVE
# -----------------------------------------------------------------------------
# THE QUESTION (correctly posed this time):
#   Dimension reduction is NOT meant to beat the uncompressed signal. It is
#   meant to shrink J factors into G blocks while (a) losing almost nothing in
#   performance and (b) massively cutting the downstream computational burden.
#   So the benchmark is the UNCOMPRESSED version, not a single factor.
#
# WHAT IS MEASURED, for each candidate G:
#   1. score agreement : mean/min daily cross-sectional RANK correlation between
#                        the G-block score and the uncompressed 24-factor score
#                        -> this is literally "how different is the output"
#   2. performance     : Rank IC mean / ICIR on TRAIN and on VALIDATION
#   3. paired test     : HAC t of the daily IC DIFFERENCE vs the 24-factor score
#                        (a matching ICIR with an insignificant difference is
#                         exactly the "no loss" claim we want to support)
#   4. downstream cost : #parameters = G(G+1)/2 and cond(Sigma_IC), measured
#                        against the effective sample T_eff ~ T/21 (the forward
#                        returns overlap by 21 days, so the IC series has only
#                        ~T/21 independent observations)
#
# DELIBERATELY INCLUDED for the adversarial review that follows:
#   - TRAIN and VALID are reported side by side (is the story window-dependent?)
#   - the MINIMUM daily score correlation, not just the mean (a reviewer will
#     ask whether some days differ wildly)
#   - within-block coherence: mean |exposure corr| inside vs across blocks
#   - a PLACEBO partition (random blocks at the same G) as a control
#
# OOS stays frozen. ASCII output only.
#
# Usage (run from the project directory):
#   $env:DFZQ_PANEL="_factor_panel2.parquet"
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_stage3_compress.R" `
#       | Tee-Object "_stage3_compress.log"
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
P("[0] %d rows | %d factors | TRAIN/VALID/OOS days = %d/%d/%d  (OOS frozen)",
  nrow(m), J, uniqueN(m[win == "TRAIN", date]), uniqueN(m[win == "VALID", date]),
  uniqueN(m[win == "OOS", date]))
T_eff <- uniqueN(m[win == "TRAIN", date]) / H
P("[0] effective training sample T_eff = T/H = %.1f independent observations", T_eff)

## ================= 1. TRAIN SIGNS =================
hr("1. SIGNS FROM TRAIN ONLY")
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
P("[1] signs aligned on TRAIN for all %d factors", J)

## ================= 2. SCORE MACHINERY =================
zcol <- function(f) { x <- m[[f]]; x[!is.finite(x)] <- 0; x }
score_eff <- function(e) {
  v <- rep(0, nrow(m))
  for (j in seq_along(FACS)) if (e[j] != 0) v <- v + e[j] * zcol(FACS[j])
  v
}
eff_block <- function(bm, w) {           # w: numeric over sort(unique(bm))
  e <- numeric(J); names(e) <- FACS; lev <- sort(unique(bm))
  for (g in seq_along(lev)) { fs <- FACS[bm == lev[g]]; e[fs] <- w[g] * sgn[fs] / length(fs) }
  e
}
setcol <- function(nm, v) { if (is.null(m[[nm]])) m[, (nm) := NA_real_]; set(m, j = nm, value = as.numeric(v)) }

ic_by_day <- function(scv, w) {
  setcol("sc_tmp", scv)
  d <- m[win == w & is.finite(sc_tmp) & is.finite(fwd),
         .(ic = if (.N < 30L) NA_real_ else
             suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = date]
  d[is.finite(ic)][order(date)]
}
KEEPWIN <- if (nzchar(Sys.getenv("DFZQ_TOUCH_OOS", ""))) c("TRAIN", "VALID", "OOS") else c("TRAIN", "VALID")
perf <- function(scv) {
  setcol("sc_tmp", scv)
  d <- m[is.finite(sc_tmp) & is.finite(fwd),
         .(ic = if (.N < 30L) NA_real_ else
             suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))),
         by = .(win, date)]
  d <- d[is.finite(ic)]
  d[, .(ic_mean = mean(ic), icir = mean(ic) / sd(ic),
        t_hac = mean(ic) / nw_se(ic, 5L), days = .N), by = win][win %in% KEEPWIN]
}
agree <- function(a, b, w) {              # daily cross-sectional rank corr of two scores
  setcol("sa_", a); setcol("sb_", b)
  d <- m[win == w & is.finite(sa_) & is.finite(sb_),
         .(r = suppressWarnings(cor(sa_, sb_, method = "spearman"))), by = date]
  d$r[is.finite(d$r)]
}

## ================= 3. BLOCKS =================
hr("3. EXPOSURE-CORRELATION WARD CLUSTERING (TRAIN ONLY)")
tr_days <- sort(unique(m[win == "TRAIN", date]))
sd_days <- tr_days[seq(1L, length(tr_days), by = 3L)]
CM <- cor(as.matrix(m[date %in% sd_days, ..FACS]), use = "pairwise.complete.obs")
dimnames(CM) <- list(FACS, FACS)
hc <- hclust(as.dist(1 - CM), method = "ward.D2")
P("[3] similarity built from %d days | mean|rho| = %.3f", length(sd_days),
  mean(abs(CM[upper.tri(CM)])))

## ================= 4. REFERENCE = UNCOMPRESSED 24-FACTOR SCORE =================
hr("4. REFERENCE: uncompressed equal-weight 24-factor score (= the G=J case)")
e_all <- rep(1 / J, J); names(e_all) <- FACS
s_ref <- score_eff(e_all * sgn)
ref <- perf(s_ref)
print(ref[, .(win, ic_mean = round(ic_mean, 4), icir = round(icir, 3),
              t_hac = round(t_hac, 2), days)])
ref_ic_va <- ic_by_day(s_ref, "VALID")
P("[4] this row is the benchmark: compression must match IT, not a single factor.")

## ================= 5. THE COMPRESSION CURVE =================
hr("5. COMPRESSION CURVE (block-equal-weight score vs the uncompressed score)")
grid <- c(24L, 15L, 12L, 10L, 8L, 6L, 5L, 4L, 3L, 2L, 1L)
blk_ic <- function(bm) {                   # G x T matrix of block IC series on TRAIN
  lev <- sort(unique(bm))
  out <- lapply(seq_along(lev), function(g) {
    w <- numeric(length(lev)); w[g] <- 1
    setcol("sc_tmp", score_eff(eff_block(bm, w)))
    m[win == "TRAIN" & is.finite(sc_tmp) & is.finite(fwd),
      .(ic = if (.N < 30L) NA_real_ else
          suppressWarnings(cor(sc_tmp, fwd, method = "spearman"))), by = date]
  })
  dt <- rbindlist(lapply(seq_along(out), function(i) out[[i]][, .(date, ic)][, block := i][]))
  dcast(dt, date ~ block, value.var = "ic")
}
res <- rbindlist(lapply(grid, function(G) {
  bm <- if (G >= J) setNames(seq_len(J), FACS) else cutree(hc, k = G)
  lev <- sort(unique(bm))
  s_g <- score_eff(eff_block(bm, rep(1 / length(lev), length(lev))))
  pv <- perf(s_g)
  # score agreement with the uncompressed score
  r_tr <- agree(s_g, s_ref, "TRAIN"); r_va <- agree(s_g, s_ref, "VALID")
  # paired test vs uncompressed, on VALID daily IC
  iv <- ic_by_day(s_g, "VALID")
  d <- merge(iv[, .(date, ia = ic)], ref_ic_va[, .(date, ib = ic)], by = "date")
  df <- d$ia - d$ib
  tv <- if (sd(df) > 0) mean(df) / nw_se(df, 5L) else NA_real_
  # downstream cost
  BM <- as.matrix(blk_ic(bm)[, -1]); Sg <- cov(BM, use = "pairwise.complete.obs")
  data.table(
    G = length(lev),
    corr_TRAIN_mean = mean(r_tr), corr_TRAIN_min = min(r_tr),
    corr_VALID_mean = mean(r_va), corr_VALID_min = min(r_va),
    ic_TRAIN = pv[win == "TRAIN", ic_mean], icir_TRAIN = pv[win == "TRAIN", icir],
    ic_VALID = pv[win == "VALID", ic_mean], icir_VALID = pv[win == "VALID", icir],
    dIC_VALID = mean(df), t_vs_uncompressed = tv,
    cond_Sigma = kappa(Sg), params = length(lev) * (length(lev) + 1) / 2,
    params_per_obs = (length(lev) * (length(lev) + 1) / 2) / T_eff)
}))
print(res[, .(G, corr_VALID_mean = round(corr_VALID_mean, 3),
              corr_VALID_min = round(corr_VALID_min, 3),
              ic_VALID = round(ic_VALID, 4), icir_VALID = round(icir_VALID, 3),
              dIC_VALID = round(dIC_VALID, 5), t_vs_uncompressed = round(t_vs_uncompressed, 2),
              cond_Sigma = signif(cond_Sigma, 3), params, params_per_obs = round(params_per_obs, 2))])
P("")
P("  REFERENCE ROW (uncompressed, G = J = %d): ic_VALID %.4f | icir_VALID %.3f | params %d",
  J, ref[win == "VALID", ic_mean], ref[win == "VALID", icir], J * (J + 1) / 2)
P("  READ: a good G has corr_VALID_mean high, t_vs_uncompressed |t| < 2 (no loss),")
P("        and params_per_obs well below 1 (downstream estimation is well-posed).")

## ================= 6. BLOCK COHERENCE + PLACEBO =================
hr("6. SANITY / PLACEBO -- is the partition doing anything the data supports?")
coh <- function(bm) {
  lev <- sort(unique(bm)); within <- c(); across <- c()
  for (g in seq_along(lev)) {
    fs <- FACS[bm == lev[g]]
    if (length(fs) > 1) within <- c(within, abs(CM[fs, fs][upper.tri(CM[fs, fs])]))
    across <- c(across, abs(CM[fs, setdiff(FACS, fs)]))
  }
  c(within = mean(within), across = mean(across))
}
for (G in c(6L, 10L, 12L)) {
  bm <- cutree(hc, k = G); co <- coh(bm)
  P("  exposure Ward G=%2d : mean|rho| within blocks %.3f vs across %.3f", G, co[1], co[2])
}
set.seed(20261011L)
for (G in c(6L, 10L)) {
  # Placebo must keep the SAME block sizes as the Ward partition, otherwise a
  # balanced partition makes block-equal-weight identical to plain equal weight
  # (e_f = sgn_f / J whenever every block has size J/G) and the control tests
  # nothing at all -- which is what happened in the first run at G = 6.
  perm <- sample(J); CMp <- CM[perm, perm]; dimnames(CMp) <- list(FACS, FACS)
  rp <- cutree(hclust(as.dist(1 - CMp), method = "ward.D2"), k = G)
  s_p <- score_eff(eff_block(rp, rep(1 / G, G)))
  pp <- perf(s_p); iv <- ic_by_day(s_p, "VALID")
  d <- merge(iv[, .(date, ia = ic)], ref_ic_va[, .(date, ib = ic)], by = "date")
  df <- d$ia - d$ib
  P("  RANDOM partition G=%2d: ic_VALID %.4f | icir_VALID %.3f | dIC %.5f (t=%.2f)",
    G, pp[win == "VALID", ic_mean], pp[win == "VALID", icir],
    mean(df), mean(df) / nw_se(df, 5L))
}
P("  READ: if a RANDOM partition matches the Ward partition, the clustering is")
P("        not adding information -- only the compression itself is.")

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
