# =============================================================================
# _stage1_cluster.R -- Signal-line Stage 1: cluster the 24 factors into G blocks
# -----------------------------------------------------------------------------
# QUESTION: how many blocks should the 24 factors be compressed into, and on
#           WHICH similarity measure? The plan originally fixed "IC correlation"
#           with exposure correlation as a control; _diag_dim.R showed the two
#           are nearly unrelated (corr = 0.19), so this script runs BOTH and
#           compares them side by side.
#
# METHOD (per the plan, Stage 1):
#   1. similarity matrix on the TRAINING window only (Sec 8.3: never touch
#      validation / frozen OOS when estimating)
#   2. Ward hierarchical clustering on (1 - corr)
#   3. walk the dendrogram from G = J down to G = 1, score every G with
#        - BIC = n*ln(RSS/n) + k*ln(n),  k = G   (the plan's formula)
#          RSS from the daily cross-sectional regression of forward return on
#          the BLOCK exposures (equal-weight within block). Computed from
#          per-day sufficient statistics X'X / X'y, transformed by the
#          aggregation matrix A, so no per-G re-scan of the data is needed.
#        - composite ICIR on TRAIN and on VALIDATION (equal-weight blocks,
#          within-block equal weight, block signs from TRAIN only)
#   4. split-half stability of the partition
#
# WHY BOTH CRITERIA ARE REPORTED -- do not trust BIC alone here:
#   with n ~ 4e6 pooled stock-days, the penalty k*ln(n) is only ~15, while
#   every extra block keeps cutting RSS. BIC is therefore expected to
#   degenerate to G = J. This is the same disease as the previous round's
#   MP test: the verdict is driven by N, not by structure. The validation-set
#   ICIR is the tie-breaker (Sec 8.3 explicitly allows comparing a small set
#   of candidate configurations on the validation window).
#
# Usage (run from the project directory):
#   $env:DFZQ_PANEL="_factor_panel2.parquet"
#   & "C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe" "_stage1_cluster.R" `
#       | Tee-Object "_stage1_cluster.log"
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
                 if (nzchar(s)) as.Date(s) else as.Date("2017-01-01") })  # 2016 = warm-up
TRAIN_END <- as.Date("2022-01-01")
VALID_END <- as.Date("2024-01-01")

## ================= 0. LOAD =================
hr("0. LOAD")
fp <- as.data.table(arrow::read_parquet(PANEL))
fp[, date := as.Date(date)]
FACS <- setdiff(names(fp), c("date", "code", "yr"))
J <- length(FACS)

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
P("[0] panel %s | %d rows | %d factors | %s ~ %s", PANEL, nrow(m), J,
  min(m$date), max(m$date))
P("[0] days per window: %s", paste(sapply(c("TRAIN", "VALID", "OOS"),
    function(w) sprintf("%s=%d", w, uniqueN(m[win == w, date]))), collapse = "  "))
P("[0] J = %d factors: %s", J, paste(FACS, collapse = ", "))

## ================= 1. PER-FACTOR IC BY WINDOW =================
hr("1. PER-FACTOR IC (memo: needed for block signs and the IC-based similarity)")
setkey(m, date)
ic <- rbindlist(lapply(FACS, function(f) {
  d <- m[is.finite(get(f)), .(win, date, z = get(f), fwd)]
  d[is.finite(z) & is.finite(fwd),
    .(ic = if (.N < 30L) NA_real_ else
        suppressWarnings(cor(z, fwd, method = "spearman"))),
    by = .(win, date)][, factor := f][]
}))
ic <- ic[is.finite(ic)]
P("[1] IC series: %d rows (%.0f s)", nrow(ic), el())

sm <- ic[, .(n = .N, ic_mean = mean(ic), icir = mean(ic) / sd(ic),
             t_hac = mean(ic) / nw_se(ic, 5L)), by = .(win, factor)]
sm[, a := abs(t_hac)]; setorder(sm, win, -a); sm[, a := NULL]
P("[1] TRAIN IC summary (|t|>2 => usable, sign matters):")
print(sm[win == "TRAIN", .(factor, n, ic_mean = round(ic_mean, 4),
                           icir = round(icir, 3), t_hac = round(t_hac, 2))])

icm_tr <- ic[win == "TRAIN", .(m = mean(ic)), by = factor]
setkey(icm_tr, factor)

## ================= 2. TWO SIMILARITY MATRICES (TRAIN ONLY) =================
hr("2. SIMILARITY MATRICES (TRAINING WINDOW ONLY)")

# --- (A) IC correlation ---------------------------------------------------------
ict <- dcast(ic[win == "TRAIN"], date ~ factor, value.var = "ic")
CM_ic <- cor(as.matrix(ict[, -1]), use = "pairwise.complete.obs")
dimnames(CM_ic) <- list(FACS, FACS)
P("[2A] IC correlation: T = %d days | mean|rho| = %.3f | max = %.3f",
  nrow(ict), mean(abs(CM_ic[upper.tri(CM_ic)])), max(abs(CM_ic[upper.tri(CM_ic)])))

# --- (B) exposure correlation (every 3rd train day, pairwise complete) ----------
tr_days <- sort(unique(m[win == "TRAIN", date]))
sd_days <- tr_days[seq(1L, length(tr_days), by = 3L)]
Xs <- as.matrix(m[date %in% sd_days, ..FACS])
CM_ex <- cor(Xs, use = "pairwise.complete.obs")
dimnames(CM_ex) <- list(FACS, FACS)
P("[2B] exposure correlation: %d rows x %d days | mean|rho| = %.3f | max = %.3f",
  nrow(Xs), length(sd_days),
  mean(abs(CM_ex[upper.tri(CM_ex)])), max(abs(CM_ex[upper.tri(CM_ex)])))
P("[2 ] corr(the two similarity measures) on shared pairs = %.3f",
  { u <- upper.tri(CM_ic); cor(CM_ic[u], CM_ex[u]) })
rm(Xs); gc()

## ================= 3. PER-DAY SUFFICIENT STATISTICS (TRAIN) =================
hr("3. PER-DAY SUFFICIENT STATISTICS (TRAIN) -- for the block regression")
n_tr <- length(tr_days)
Gs <- vector("list", n_tr); bs <- vector("list", n_tr)
yy <- numeric(n_tr); nn <- integer(n_tr); imputed <- 0; total <- 0
for (i in seq_len(n_tr)) {
  sub <- m[.(tr_days[i])]
  X <- as.matrix(sub[, ..FACS]); y <- sub$fwd
  ok <- is.finite(y); X <- X[ok, , drop = FALSE]; y <- y[ok]
  total <- total + length(y); imputed <- imputed + sum(!is.finite(X))
  X[!is.finite(X)] <- 0                    # missing z-score -> neutral (0)
  Gs[[i]] <- crossprod(X); bs[[i]] <- crossprod(X, y)
  yy[i] <- sum(y^2); nn[i] <- length(y)
}
N <- sum(nn)
P("[3] TRAIN: %d days | n = %d stock-days | imputed rate = %.2f%%",
  n_tr, N, 100 * imputed / total)

BAD <- 0L
rss_of <- function(A) {
  tot <- 0
  for (i in seq_len(n_tr)) {
    Ag <- crossprod(A, Gs[[i]] %*% A); bb <- crossprod(A, bs[[i]])
    sol <- tryCatch(solve(Ag, bb), error = function(e) NULL)
    if (is.null(sol) || any(!is.finite(sol))) { BAD <<- BAD + 1L; next }
    tot <- tot + yy[i] - sum(bb * sol)
  }
  tot
}

## ================= 4. WALK THE DENDROGRAM; SCORE EVERY G =================
hr("4. BLOCK COUNT SELECTION -- BIC vs validation ICIR, for BOTH bases")

m_tv <- m[win != "OOS"]
setkey(m_tv, date)

# composite IC (Spearman) by window, given a score vector aligned to m_tv
icir_of <- function(v) {
  m_tv[, sc := v]
  d <- m_tv[is.finite(sc) & is.finite(fwd),
            .(ic = if (.N < 30L) NA_real_ else
                suppressWarnings(cor(sc, fwd, method = "spearman"))),
            by = .(win, date)]
  d <- d[is.finite(ic)]
  out <- d[, .(ic_mean = mean(ic), icir = mean(ic) / sd(ic),
               t_hac = mean(ic) / nw_se(ic, 5L), n = .N), by = win]
  out
}

run_base <- function(CM, label) {
  hc <- hclust(as.dist(1 - CM), method = "ward.D2")
  # bottom-up: leaves = single factors, internal node J+i merges hc$merge[i, ]
  ns <- vector("list", 2L * J - 1L)     # node -> score vector (NA-free)
  mem <- vector("list", 2L * J - 1L)    # node -> member factor names
  for (j in seq_len(J)) {
    # Align each factor by its OWN train IC sign FIRST. A single block-level
    # sign is WRONG whenever a block mixes positive- and negative-IC factors
    # (e.g. G = 1 mixes all 24) -- that bug made the first full run's curve
    # meaningless. After leaf alignment, block scores are already oriented.
    s0 <- sign(icm_tr[FACS[j], m]); if (!is.finite(s0) || s0 == 0) s0 <- 1
    v <- copy(m_tv[[FACS[j]]]); v[!is.finite(v)] <- 0
    ns[[j]] <- s0 * v; mem[[j]] <- FACS[j]
  }
  rows <- list()
  for (i in seq_len(J - 1L)) {
    a <- hc$merge[i, 1L]; b <- hc$merge[i, 2L]
    a <- if (a < 0L) -a else J + a
    b <- if (b < 0L) -b else J + b
    na <- length(mem[[a]]); nb <- length(mem[[b]])
    ns[[J + i]] <- (na * ns[[a]] + nb * ns[[b]]) / (na + nb)
    mem[[J + i]] <- c(mem[[a]], mem[[b]])
  }
  cur <- seq_len(J)                     # roots (node ids), start at G = J
  G <- J
  while (TRUE) {
    # block map: factor -> block id, from the current roots
    bm <- integer(J); names(bm) <- FACS
    for (g in seq_along(cur)) bm[mem[[cur[g]]]] <- g
    A <- matrix(0, J, length(cur), dimnames = list(FACS, NULL))
    for (g in seq_along(cur)) A[bm == g, g] <- 1 / sum(bm == g)
    r <- rss_of(A)
    bic <- N * log(r / N) + G * log(N)
    # composite: equal weight across blocks, sign from TRAIN IC means
    tot <- rep(0, nrow(m_tv)); sg_txt <- character(length(cur))
    for (g in seq_along(cur)) {
      mm <- mem[[cur[g]]]
      tot <- tot + ns[[cur[g]]]          # leaves are already sign-aligned
      sg_txt[g] <- paste(substr(mm, 1, 4), collapse = ",")
    }
    iv <- icir_of(tot / length(cur))
    rows[[length(rows) + 1L]] <- data.table(
      base = label, G = G, rss = r, bic = bic,
      icir_train = iv[win == "TRAIN", icir], icir_valid = iv[win == "VALID", icir],
      t_train = iv[win == "TRAIN", t_hac], t_valid = iv[win == "VALID", t_hac],
      blocks = paste(sg_txt, collapse = " | "))
    if (G == 1L) break
    # apply merge (J - G + 1) to drop from G to G-1
    i <- J - G + 1L
    a <- hc$merge[i, 1L]; b <- hc$merge[i, 2L]
    a <- if (a < 0L) -a else J + a
    b <- if (b < 0L) -b else J + b
    cur <- c(cur[!(cur %in% c(a, b))], J + i)
    G <- G - 1L
  }
  list(hc = hc, tab = rbindlist(rows))
}

res_ex <- run_base(CM_ex, "EXPOSURE")
res_ic <- run_base(CM_ic, "IC")
res <- rbind(res_ex$tab, res_ic$tab)
P("[4] scanning done (%.0f s) | singular days skipped: %d", el(), BAD)

for (lb in c("EXPOSURE", "IC")) {
  P("")
  P("---- base = %s ----", lb)
  print(res[base == lb, .(G, rss = signif(rss, 5), bic = signif(bic, 5),
                          icir_train = round(icir_train, 3),
                          icir_valid = round(icir_valid, 3),
                          t_valid = round(t_valid, 2))])
  sub   <- res[base == lb][is.finite(bic)]
  g_bic <- sub[which.min(bic)]$G
  g_val <- res[base == lb][which.max(icir_valid)]$G
  g_tr  <- res[base == lb][which.max(icir_train)]$G
  P("  argmin BIC          : G = %d   <-- the plan's rule", g_bic)
  P("  best VALID ICIR     : G = %d   <-- the tie-breaker", g_val)
  P("  best TRAIN ICIR     : G = %d", g_tr)
  P("  G = 1 ICIR (train/valid) : %.3f / %.3f",
    res[base == lb & G == 1, icir_train], res[base == lb & G == 1, icir_valid])
  b_tr <- sm[win == "TRAIN"][which.max(abs(icir))]
  b_va <- sm[win == "VALID"][which.max(abs(icir))]
  P("  best single-factor  : TRAIN |ICIR| = %.3f (%s) | VALID |ICIR| = %.3f (%s)",
    abs(b_tr$icir), b_tr$factor, abs(b_va$icir), b_va$factor)
}

## ================= 5. PARTITIONS SIDE BY SIDE =================
hr("5. BLOCK PARTITIONS AT THE VALIDATION-OPTIMAL G")
g_ex <- res[base == "EXPOSURE"][which.max(icir_valid)]$G
g_ic <- res[base == "IC"][which.max(icir_valid)]$G
P("[5] chosen G: exposure = %d | IC = %d", g_ex, g_ic)
for (lb in c("EXPOSURE", "IC")) {
  g <- if (lb == "EXPOSURE") g_ex else g_ic
  hc <- if (lb == "EXPOSURE") res_ex$hc else res_ic$hc
  bm <- cutree(hc, k = g)
  P("")
  P("---- %s base, G = %d ----", lb, g)
  for (k in sort(unique(bm))) P("  block %d: %s", k, paste(names(bm)[bm == k], collapse = ", "))
}

## ================= 6. SPLIT-HALF STABILITY =================
hr("6. SPLIT-HALF STABILITY OF THE PARTITION")
half <- tr_days[floor(length(tr_days) / 2)]
agree <- function(bm1, bm2) {
  pr <- combn(J, 2L)
  a1 <- bm1[pr[1, ]] == bm1[pr[2, ]]; a2 <- bm2[pr[1, ]] == bm2[pr[2, ]]
  mean(a1 == a2)
}
ict1 <- dcast(ic[win == "TRAIN" & date <= half], date ~ factor, value.var = "ic")
ict2 <- dcast(ic[win == "TRAIN" & date >  half], date ~ factor, value.var = "ic")
CM1 <- cor(as.matrix(ict1[, -1]), use = "pairwise.complete.obs"); dimnames(CM1) <- list(FACS, FACS)
CM2 <- cor(as.matrix(ict2[, -1]), use = "pairwise.complete.obs"); dimnames(CM2) <- list(FACS, FACS)
P("[6] IC base    : pair-agreement between halves at same G")
for (g in c(2L, 3L, 4L, 5L, 8L)) P("     G=%2d -> %.3f", g, agree(cutree(hclust(as.dist(1 - CM1), "ward.D2"), g),
                                                                   cutree(hclust(as.dist(1 - CM2), "ward.D2"), g)))
d1 <- sd_days[sd_days <= half]; d2 <- sd_days[sd_days > half]
CMe1 <- cor(as.matrix(m[date %in% d1, ..FACS]), use = "pairwise.complete.obs"); dimnames(CMe1) <- list(FACS, FACS)
CMe2 <- cor(as.matrix(m[date %in% d2, ..FACS]), use = "pairwise.complete.obs"); dimnames(CMe2) <- list(FACS, FACS)
P("[6] exposure base: pair-agreement between halves at same G")
for (g in c(2L, 3L, 4L, 5L, 8L)) P("     G=%2d -> %.3f", g, agree(cutree(hclust(as.dist(1 - CMe1), "ward.D2"), g),
                                                                   cutree(hclust(as.dist(1 - CMe2), "ward.D2"), g)))
P("[6] corr(H1, H2) of the similarity matrices themselves: IC = %.3f | exposure = %.3f",
  cor(CM1[upper.tri(CM1)], CM2[upper.tri(CM2)]),
  cor(CMe1[upper.tri(CMe1)], CMe2[upper.tri(CMe2)]))

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("READ: BIC is expected to degenerate (n is huge) -- trust the VALIDATION ICIR")
P("      and the split-half agreement. Compare which base gives blocks that are")
P("      both STABLE and ECONOMICALLY interpretable.")
