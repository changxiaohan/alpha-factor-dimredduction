# =============================================================================
# _stage5_challenge.R -- adversarial review of the "compression is harmless" claim
#                        (v2: delta = 10%, plus CAUSE DIAGNOSTICS)
# -----------------------------------------------------------------------------
# delta default 10% (env DFZQ_DELTA_PCT). Applied to two scales:
#   (a) predictive: |dIC| <= delta x reference IC   (HAC 95% CI must fit inside)
#   (b) output    : compressed/uncompressed rank correlation must stay >= 1-delta
#
# THE CENTRAL MECHANISM THIS SCRIPT EXPOSES
#   block-equal-weight + equal-across-blocks gives each factor the weight
#        e_f = sgn_f / (G * |block(f)|)
#   so "compression" is simultaneously COMPRESSION and IMPLICIT REWEIGHTING.
#   Two controls separate them:
#     NEUTRAL compression : w_g prop. to |block_g|  =>  e_f = sgn_f/J exactly,
#                           i.e. the score is IDENTICAL to plain equal weight.
#                           Pure compression therefore costs nothing in score,
#                           its only cost is the model restriction (the F test).
#     EQUAL-BLOCK         : w_g = 1/G  =>  the reweighting effect.
#   Plus: correlate each factor's weight multiplier against its TRAIN ICIR. If
#   positive, compression systematically up-weights the strong factors -- that
#   is the mechanism behind the jagged dIC curve.
#
# CHALLENGES
#   A  which reference (equal-weight vs best-uncompressed ICIR-weighted)
#   A2 neutral vs equal-block compression  (the compression / reweighting split)
#   A3 weight-multiplier vs TRAIN ICIR correlation  (cause of the jagged dIC)
#   A4 monotonicity of dIC in G (a genuine loss should be monotone)
#   B  portfolio-level agreement + the NATURAL DISAGREEMENT BAND between
#      uncompressed weightings
#   C  partition stability across 4 disjoint training sub-periods
#   D  placebo distribution + does within-block coherence relate to performance
#   F  drop the lottery/volatility cluster and redo
#   E  window swap: run the same script with $env:DFZQ_SWAP="1"
#
# Usage:  $env:DFZQ_PANEL="_factor_panel2.parquet"
#         & "...\Rscript.exe" "_stage5_challenge.R" | Tee-Object "_stage5_challenge.log"
# ASCII only. OOS stays frozen.
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
ci95 <- function(x) c(lo = mean(x) - 1.96 * nw_se(x, 5L), hi = mean(x) + 1.96 * nw_se(x, 5L))

T0 <- Sys.time(); H <- 21L
PANEL <- Sys.getenv("DFZQ_PANEL", "_factor_panel2.parquet")
START <- local({ s <- Sys.getenv("DFZQ_START", ""); if (nzchar(s)) as.Date(s) else as.Date("2017-01-01") })
B1 <- as.Date("2022-01-01"); B2 <- as.Date("2024-01-01")
DELTA <- as.numeric(Sys.getenv("DFZQ_DELTA_PCT", "10")) / 100
SWAP  <- nzchar(Sys.getenv("DFZQ_SWAP", ""))

## ================= 0. LOAD =================
hr("0. LOAD")
fp <- as.data.table(arrow::read_parquet(PANEL)); fp[, date := as.Date(date)]
FACS <- setdiff(names(fp), c("date", "code", "yr")); J <- length(FACS)
daily <- as.data.table(arrow::read_parquet("data/daily.parquet",
                                           col_select = c("date", "code", "close")))
daily[, date := as.Date(date)][, code := norm_code(code)]
setorder(daily, code, date)
daily[, fwd := shift(close, -H) / close - 1, by = code]
m <- merge(fp[date >= START], daily[is.finite(fwd), .(date, code, fwd)], by = c("date", "code"))
rm(fp, daily); gc()
if (SWAP) {
  m[, win := fifelse(date >= B1 & date < B2, "TRAIN", fifelse(date < B1, "VALID", "OOS"))]
} else {
  m[, win := fifelse(date < B1, "TRAIN", fifelse(date < B2, "VALID", "OOS"))]
}
setorder(m, date, code); setkey(m, date)
P("[0] %d rows | %d factors | delta = %.0f%% | WINDOW MODE = %s",
  nrow(m), J, 100 * DELTA, if (SWAP) "SWAPPED (train on the LATER window)" else "normal")
P("[0] TRAIN %s ~ %s (%d days) | VALID %s ~ %s (%d days) | OOS %d days frozen",
  min(m[win == "TRAIN", date]), max(m[win == "TRAIN", date]), uniqueN(m[win == "TRAIN", date]),
  min(m[win == "VALID", date]), max(m[win == "VALID", date]), uniqueN(m[win == "VALID", date]),
  uniqueN(m[win == "OOS", date]))

## ---- helpers ---------------------------------------------------------------
zcol <- function(f) { x <- m[[f]]; x[!is.finite(x)] <- 0; x }
setcol <- function(nm, v) { if (is.null(m[[nm]])) m[, (nm) := NA_real_]; set(m, j = nm, value = as.numeric(v)) }
score_eff <- function(e) { v <- rep(0, nrow(m))
  for (j in seq_along(FACS)) if (e[j] != 0) v <- v + e[j] * zcol(FACS[j]); v }
ic_series <- function(scv, w = "VALID") { setcol("s_", scv)
  d <- m[win == w & is.finite(s_) & is.finite(fwd),
         .(ic = if (.N < 30L) NA_real_ else suppressWarnings(cor(s_, fwd, method = "spearman"))),
         by = date]
  d[is.finite(ic)][order(date)] }
rankcorr <- function(a, b, w = "VALID") { setcol("sa_", a); setcol("sb_", b)
  d <- m[win == w & is.finite(sa_) & is.finite(sb_),
         .(r = suppressWarnings(cor(sa_, sb_, method = "spearman"))), by = date]
  d$r[is.finite(d$r)] }
topsets <- function(scv, k = 100L) { setcol("s_", scv)
  m[win == "VALID" & is.finite(s_), ][order(date, -s_), head(.SD, k), by = date][, .(date, code)] }
ov <- function(a, b) mean(paste(a$date, a$code) %in% paste(b$date, b$code))
turnover_of <- function(t) { setorder(t, date, code); dts <- sort(unique(t$date))
  prev <- t[date == dts[1], code]; tot <- 0
  for (i in 2:length(dts)) { cur <- t[date == dts[i], code]
    tot <- tot + 1 - sum(cur %in% prev) / length(cur); prev <- cur }
  tot / (length(dts) - 1) }

ic <- rbindlist(lapply(FACS, function(f) {
  d <- m[is.finite(get(f)), .(win, date, z = get(f), fwd)]
  d[is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = .(win, date)][, factor := f][]
}))
ic <- ic[is.finite(ic)]
icm <- ic[win == "TRAIN", .(m = mean(ic)), by = factor]; setkey(icm, factor)
icir_tr <- ic[win == "TRAIN", .(icir = mean(ic) / sd(ic)), by = factor]; setkey(icir_tr, factor)
sgn <- vapply(FACS, function(f) { v <- sign(icm[f, m]); if (!is.finite(v) || v == 0) 1 else v }, numeric(1))
names(sgn) <- FACS

o_le <- rep(1 / J, J); names(o_le) <- FACS
o_opt <- icir_tr[FACS, icir]; o_opt <- o_opt / sum(abs(o_opt))
ref_le  <- score_eff(o_le * sgn); ref_opt <- score_eff(o_opt)
ic_le <- ic_series(ref_le); ic_opt <- ic_series(ref_opt)
P("[0] VALID mean IC: equal-weight 24 = %.4f | ICIR-weighted 24 = %.4f (the 'best uncompressed')",
  mean(ic_le$ic), mean(ic_opt$ic))

tr_days <- sort(unique(m[win == "TRAIN", date]))
sd_days <- tr_days[seq(1L, length(tr_days), by = 3L)]
CM <- cor(as.matrix(m[date %in% sd_days, ..FACS]), use = "pairwise.complete.obs")
dimnames(CM) <- list(FACS, FACS)
hc <- hclust(as.dist(1 - CM), method = "ward.D2")

# weight multiplier of each factor vs plain equal weight, for a block map + block w
mult_of <- function(bm, w) { lev <- sort(unique(bm)); mu <- numeric(J); names(mu) <- FACS
  for (g in seq_along(lev)) { mem <- FACS[bm == lev[g]]; mu[mem] <- w[g] / length(mem) }
  mu * J }
score_from <- function(bm, w) { mu <- mult_of(bm, w); score_eff((mu / J) * sgn) }
w_eqblk  <- function(bm) { lev <- sort(unique(bm)); rep(1 / length(lev), length(lev)) }   # 1/G
w_neutral<- function(bm) { lev <- sort(unique(bm)); n <- sapply(lev, function(g) sum(bm == g)); n / sum(n) }
cohgap <- function(bm) { lev <- sort(unique(bm)); wi <- c(); ac <- c()
  for (g in seq_along(lev)) { fs <- FACS[bm == lev[g]]
    if (length(fs) > 1) wi <- c(wi, abs(CM[fs, fs][upper.tri(CM[fs, fs])]))
    ac <- c(ac, abs(CM[fs, setdiff(FACS, fs)])) }
  mean(wi) - mean(ac) }

GRID <- c(12L, 10L, 8L, 6L, 4L, 3L, 2L)

## ================= 3. CHALLENGE A + the compression/reweighting split ======
hr("3. CHALLENGE A -- reference choice, and the compression/reweighting split")
P("  ref_eq  VALID IC %.4f | ref_opt (best uncompressed) %.4f", mean(ic_le$ic), mean(ic_opt$ic))
P("  delta band: %.5f (vs ref_eq) | %.5f (vs ref_opt)",
  DELTA * mean(ic_le$ic), DELTA * mean(ic_opt$ic))

A <- rbindlist(lapply(GRID, function(G) {
  bm <- cutree(hc, k = G); lev <- sort(unique(bm))
  s_eqblk <- score_from(bm, w_eqblk(bm))
  s_neut  <- score_from(bm, w_neutral(bm))
  iv  <- ic_series(s_eqblk); ivn <- ic_series(s_neut)
  d1 <- merge(iv, ic_le,  by = "date", suffixes = c(".g", ".e"))
  d2 <- merge(iv, ic_opt, by = "date", suffixes = c(".g", ".o"))
  dn <- merge(ivn, ic_le, by = "date", suffixes = c(".g", ".e"))
  f1 <- d1$ic.g - d1$ic.e; f2 <- d2$ic.g - d2$ic.o; fn <- dn$ic.g - dn$ic.e
  mu <- mult_of(bm, w_eqblk(bm))
  data.table(G = G,
             dIC_vs_eq = mean(f1), lo_eq = ci95(f1)["lo"], hi_eq = ci95(f1)["hi"],
             NI_eq = ci95(f1)["lo"] > -DELTA * mean(ic_le$ic),
             dIC_vs_opt = mean(f2), lo_opt = ci95(f2)["lo"], hi_opt = ci95(f2)["hi"],
             NI_opt = ci95(f2)["lo"] > -DELTA * mean(ic_opt$ic),
             dIC_NEUTRAL = mean(fn), maxabs_neur = max(abs(fn)),
             mult_spread = sd(log(mu)),
             mult_vs_icir = suppressWarnings(cor(mu, icir_tr[FACS, icir])))
}))
print(A[, .(G, dIC_vs_eq = round(dIC_vs_eq, 5),
            CI_eq = sprintf("[%+.5f,%+.5f]", lo_eq, hi_eq), NI_eq,
            dIC_vs_opt = round(dIC_vs_opt, 5), NI_opt,
            mult_spread = round(mult_spread, 3), mult_vs_icir = round(mult_vs_icir, 3))])
P("")
P("  NEUTRAL compression (block weight prop. to block size):")
print(A[, .(G, dIC_NEUTRAL = signif(dIC_NEUTRAL, 2), max_abs_daily = signif(maxabs_neur, 2))])
P("  READ: NEUTRAL must be ~0 by construction => pure compression costs NOTHING in")
P("        score; everything in dIC_vs_eq is REWEIGHTING. mult_vs_icir > 0 means")
P("        compression systematically up-weights the strong factors (the cause of")
P("        the jagged dIC curve).")

## ---- A4 monotonicity -------------------------------------------------------
hr("3b. A4 -- is the loss MONOTONE in G? (a real information loss should be)")
ok <- A[is.finite(dIC_vs_eq)]
P("  Spearman(G, dIC vs ref_eq) = %.3f   (monotone loss would be strongly POSITIVE)",
  suppressWarnings(cor(ok$G, ok$dIC_vs_eq, method = "spearman")))
P("  Spearman(G, dIC vs ref_opt) = %.3f",
  suppressWarnings(cor(ok$G, ok$dIC_vs_opt, method = "spearman")))

## ================= 4. CHALLENGE B ==========================================
hr("4. CHALLENGE B -- portfolio agreement, and the NATURAL DISAGREEMENT BAND")
o_ic <- ic[win == "TRAIN", .(w = mean(ic)), by = factor]
w_ic <- setNames(o_ic$w, o_ic$factor)[FACS]; w_ic <- w_ic / sum(abs(w_ic))
ref_ic_w <- score_eff(w_ic)
t_le <- topsets(ref_le); t_opt <- topsets(ref_opt); t_ic <- topsets(ref_ic_w)
P("  UNCOMPRESSED variants (all 24 columns) -- pairwise top-100 overlap:")
P("    equal vs ICIR-wtd = %.3f | equal vs IC-wtd = %.3f | ICIR vs IC = %.3f",
  ov(t_le, t_opt), ov(t_le, t_ic), ov(t_opt, t_ic))
P("    daily turnover: equal %.3f | ICIR %.3f | IC %.3f",
  turnover_of(t_le), turnover_of(t_opt), turnover_of(t_ic))
P("  ^ that is the 'natural disagreement band': different weightings of the SAME")
P("    24 factors already disagree by this much.")
B <- rbindlist(lapply(GRID, function(G) {
  bm <- cutree(hc, k = G); t_g <- topsets(score_from(bm, w_eqblk(bm)))
  data.table(G = G, ov_le = ov(t_g, t_le), ov_opt = ov(t_g, t_opt), ov_ic = ov(t_g, t_ic),
             turnover = turnover_of(t_g)) }))
print(B[, .(G, ov_le = round(ov_le, 3), ov_opt = round(ov_opt, 3),
            ov_ic = round(ov_ic, 3), turnover = round(turnover, 3))])

## ================= 5. CHALLENGE C ==========================================
hr("5. CHALLENGE C -- partition stability across 4 disjoint training sub-periods")
ch <- split(tr_days, cut(seq_along(tr_days), 4L, labels = FALSE))
parts <- lapply(ch, function(dd) { ds <- dd[seq(1L, length(dd), by = 3L)]
  Cx <- cor(as.matrix(m[date %in% ds, ..FACS]), use = "pairwise.complete.obs")
  dimnames(Cx) <- list(FACS, FACS)
  cutree(hclust(as.dist(1 - Cx), method = "ward.D2"), k = 10L) })
agree2 <- function(a, b) { pr <- combn(J, 2L); mean((a[pr[1, ]] == a[pr[2, ]]) == (b[pr[1, ]] == b[pr[2, ]])) }
P("  pairwise agreement at G=10 across the 4 sub-periods: %s",
  paste(sprintf("%.3f", unlist(sapply(1:3, function(i) sapply((i + 1):4, function(k) agree2(parts[[i]], parts[[k]]))))),
        collapse = ", "))

## ================= 6. CHALLENGE D ==========================================
hr("6. CHALLENGE D -- placebo distribution AND does coherence drive performance?")
set.seed(20261011L)
for (G in c(10L, 8L)) {
  bm <- cutree(hc, k = G); w_ward <- mean(ic_series(score_from(bm, w_eqblk(bm)))$ic)
  gap_ward <- cohgap(bm)
  r <- numeric(30); gp <- numeric(30)
  for (s in 1:30) {
    perm <- sample(J); Cp <- CM[perm, perm]; dimnames(Cp) <- list(FACS, FACS)
    rp <- cutree(hclust(as.dist(1 - Cp), method = "ward.D2"), k = G)
    r[s] <- mean(ic_series(score_from(rp, w_eqblk(rp)))$ic); gp[s] <- cohgap(rp)
  }
  P("  G=%2d : Ward %.4f | random %.4f +- %.4f | Ward percentile %.0f%%",
    G, w_ward, mean(r), sd(r), 100 * mean(r < w_ward))
  P("         coherence gap: Ward %.3f | random %.3f +- %.3f",
    gap_ward, mean(gp), sd(gp))
  P("         corr(random partitions: coherence gap vs VALID IC) = %.3f",
    suppressWarnings(cor(gp, r)))
}
P("  READ: if that last correlation is ~0, the clustering CRITERION (high within-")
P("        block correlation) has no relation to performance -- Ward's edge, when")
P("        it exists, comes only from the block-SIZE distribution, not the grouping.")

## ================= 7. CHALLENGE F ==========================================
hr("7. CHALLENGE F -- drop the lottery/volatility cluster and redo")
VOL <- c("maxret21", "ivol60", "volatility", "liquidity", "hlrange60", "beta60", "skew60")
F2 <- setdiff(FACS, VOL)
P("  dropped %d | remaining J' = %d : %s", length(VOL), length(F2), paste(F2, collapse = ", "))
sgn2 <- setNames(sgn[F2], F2)
cm2 <- cor(as.matrix(m[date %in% sd_days, ..F2]), use = "pairwise.complete.obs")
dimnames(cm2) <- list(F2, F2)
off2 <- abs(cm2[upper.tri(cm2)]); ev2 <- eigen(cm2, symmetric = TRUE)$values
P("  remaining pool: mean|rho| %.3f | |rho|>0.5 %.0f%% | PR %.2f / %d",
  mean(off2), 100 * mean(off2 > 0.5), sum(ev2)^2 / sum(ev2^2), length(F2))
hc2 <- hclust(as.dist(1 - cm2), method = "ward.D2")
o2 <- icir_tr[F2, icir]; o2 <- o2 / sum(abs(o2))
setcol("sc_r", score_eff({ e <- numeric(J); names(e) <- FACS; e[F2] <- o2; e }))
ref2 <- m[win == "VALID" & is.finite(sc_r) & is.finite(fwd),
          .(ic = suppressWarnings(cor(sc_r, fwd, method = "spearman"))), by = date][is.finite(ic)]
P("  reduced-pool uncompressed IC %.4f (full pool %.4f)", mean(ref2$ic), mean(ic_opt$ic))
for (G in c(6L, 4L, 3L)) {
  bm2 <- cutree(hc2, k = G); lev2 <- sort(unique(bm2))
  e <- numeric(J); names(e) <- FACS
  for (g in seq_along(lev2)) { mem <- F2[bm2 == lev2[g]]; e[mem] <- sgn2[mem] / (length(lev2) * length(mem)) }
  iv <- ic_series(score_eff(e)); d <- merge(iv, ref2, by = "date", suffixes = c(".g", ".r"))
  df <- d$ic.g - d$ic.r
  P("  reduced G=%2d : IC %.4f | dIC %+.5f CI [%+.5f,%+.5f] | NI(%d%%): %s",
    G, mean(iv$ic), mean(df), ci95(df)["lo"], ci95(df)["hi"], 100 * DELTA,
    ifelse(ci95(df)["lo"] > -DELTA * abs(mean(ref2$ic)), "YES", "no"))
}

## ================= 8. CONSTRUCTIVE TEST: replace the clustering criterion ===
hr("8. CONSTRUCTIVE TEST -- ICIR-BALANCED blocks instead of correlation blocks")
P("  Section 6 shows Ward's criterion (max within-block correlation) is unrelated")
P("  or negatively related to performance. Alternative criterion: group factors so")
P("  that every block carries a similar SHARE of TRAIN |ICIR| (functional blocks),")
P("  greedy-balanced, instead of grouping by similarity.")
for (G in c(6L, 4L, 3L, 2L)) {
  ord <- FACS[order(-abs(icir_tr[FACS, icir]))]
  sum_ic <- numeric(G); grp <- integer(J); names(grp) <- FACS
  for (f in ord) { g <- which.min(sum_ic); grp[f] <- g; sum_ic[g] <- sum_ic[g] + abs(icir_tr[f, icir]) }
  s_w <- score_from(grp, w_eqblk(grp)); iv <- ic_series(s_w)
  d1 <- merge(iv, ic_le,  by = "date", suffixes = c(".g", ".e")); f1 <- d1$ic.g - d1$ic.e
  d2 <- merge(iv, ic_opt, by = "date", suffixes = c(".g", ".o")); f2 <- d2$ic.g - d2$ic.o
  P("  ICIR-bal G=%2d : IC %.4f | vs eq  %+.5f CI[%+.5f,%+.5f] NI:%s | vs opt %+.5f NI:%s",
    G, mean(iv$ic), mean(f1), ci95(f1)["lo"], ci95(f1)["hi"],
    ifelse(ci95(f1)["lo"] > -DELTA * mean(ic_le$ic), "YES", "no"),
    mean(f2), ifelse(ci95(f2)["lo"] > -DELTA * mean(ic_opt$ic), "YES", "no"))
  P("        blocks (%s): %s", paste(sprintf("%.2f", sum_ic), collapse = "/"),
    paste(sapply(sort(unique(grp)), function(g) paste(FACS[grp == g], collapse = "+")), collapse = " | "))
}

hr("DONE")
P("Total elapsed %.0f s (%.1f min) | delta = %.0f%% | mode = %s",
  el(), el() / 60, 100 * DELTA, if (SWAP) "SWAPPED" else "normal")
