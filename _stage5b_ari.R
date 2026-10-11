# =============================================================================
# _stage5b_ari.R -- RI and ARI for the block partition across 4 training sub-periods
# -----------------------------------------------------------------------------
# The earlier report quoted the pair-agreement statistic, which is the RAND INDEX
# (RI). RI is NOT chance-corrected, so with J = 24 factors cut into G = 10 blocks
# it reads optimistically: a lot of pairs land "together" by chance. The honest
# number is the ADJUSTED Rand Index (ARI). This script reports both, so the
# report can quote the right one.
#
#   RI  = (A + D) / C(n,2),  A = #pairs together in both, D = #pairs apart in both
#   ARI = (A - E[A]) / (max(A) - E[A]),  E[A] = S_i * S_j / C(n,2)
#   where S_i = sum_i C(a_i,2) over blocks of partition 1, S_j likewise for 2.
#   ARI: 1 = identical, 0 = chance level, negative = worse than chance.
#
# Usage: & "...\Rscript.exe" "_stage5b_ari.R" | Tee-Object "_stage5b_ari.log"
# Output: console only, no files written. ASCII output only.
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")

fp <- as.data.table(arrow::read_parquet(Sys.getenv("DFZQ_PANEL", "_factor_panel2.parquet")))
fp[, date := as.Date(date)]
FACS <- setdiff(names(fp), c("date", "code", "yr")); J <- length(FACS)
m <- fp[date >= as.Date("2017-01-01") & date < as.Date("2022-01-01")]
tr_days <- sort(unique(m$date))
P("train days = %d | J = %d factors", length(tr_days), J)

honest <- function(a, b) {                 # contingency-table form
  n <- length(a); tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  tot <- c2(n)
  A  <- sum(c2(as.vector(tab)))            # pairs together in BOTH
  Si <- sum(c2(rowSums(tab))); Sj <- sum(c2(colSums(tab)))
  D  <- tot - Si - Sj + A                  # pairs apart in BOTH
  RI <- (A + D) / tot
  EA <- Si * Sj / tot
  ARI <- (A - EA) / (0.5 * (Si + Sj) - EA)
  c(RI = RI, ARI = ARI)
}
brute <- function(a, b) {                  # check RI by direct pair enumeration
  pr <- combn(J, 2L); mean((a[pr[1, ]] == a[pr[2, ]]) == (b[pr[1, ]] == b[pr[2, ]]))
}

hr("PARTITION AGREEMENT ACROSS 4 TRAINING SUB-PERIODS (G = 10)")
for (G in c(10L, 6L)) {
  brk <- split(tr_days, cut(seq_along(tr_days), 4L, labels = FALSE))
  parts <- lapply(brk, function(dd) {
    ds <- dd[seq(1L, length(dd), by = 3L)]
    Cx <- cor(as.matrix(m[date %in% ds, ..FACS]), use = "pairwise.complete.obs")
    dimnames(Cx) <- list(FACS, FACS)
    cutree(hclust(as.dist(1 - Cx), method = "ward.D2"), k = G)
  })
  P("")
  P("G = %d   (block sizes per sub-period: %s)", G,
    paste(sapply(parts, function(p) paste(sort(as.integer(table(p))), collapse = "/")), collapse = "  "))
  P("  pair |   RI   |   ARI  | RI(brute-force) | RI-RI = 0 ?")
  for (i in 1:3) for (k in (i + 1):4) {
    s <- honest(parts[[i]], parts[[k]]); bf <- brute(parts[[i]], parts[[k]])
    P("  %d-%d  | %.3f  | %.3f |     %.3f       | %s",
      i, k, s["RI"], s["ARI"], bf, ifelse(abs(bf - s["RI"]) < 1e-9, "OK", "MISMATCH"))
  }
  P("  READ: ARI is the chance-corrected number; quote ARI, not RI.")
}
cat("\nDONE\n")
