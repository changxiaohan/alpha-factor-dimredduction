# =============================================================================
# _crosssec_wls.R  ——  风险线：横截面因子收益估计（优化点 B）+ 风险矩阵 Σ
# -----------------------------------------------------------------------------
# 对应《横截面风险估计_优化点B_思路流程》：
#   Stage0 数据对齐 → Stage1 特异方差 σ² → Stage2【唯一改动】ridge → WLS
#   → Stage3 下游 F/D/Σ → Stage4 三层评价（L1 f 预测 / L2 Σ 风险精度 / L3 组合层）
#
# 用法：
#   完整流程 :  Rscript _crosssec_wls.R
#   单元测试 :  先 Sys.setenv(DFZQ_SOURCE_ONLY="1") 再 source()
#
# ⚠️ 口径（同《数据获取方案》）：
#   * close 是前复权价 → 仅用于收益，不用于算市值
#   * 【2026-10-11 更正 / 见 _probe_units.R】字段语义实测：
#       amount / volume   = 价格（元/股）
#       amount / turnover = 流通市值（元）
#       volume / turnover = 流通股本
#     旧注释只写"市值用 amount/turnover"，而 _factor_ext.R 把它当成"价格"去算
#     B/P、E/P，量纲差一个股本因子 → 价值因子 IC 被压没（t 仅 1.78）。
#     改用真价格后 value_bp 的 t 升到 8.14。凡涉及价格/市值/换手，先跑 _probe_units.R。
#   * 暴露矩阵 B = [行业哑变量(0/1) | 风格暴露(已横截面标准化)]
#   * 横截面回归是【恰好识别】(K 矩/K 参数) → GMM 与权重 W 无关；
#     WLS 不是"换 W"，而是"换最优工具 z_i = B_i/σ_i²"
#
# ⚠️ 2026-10-07 修复记录（详见《风险线_横截面WLS_代码审查》）：
#   [H-1] build_style_exposures 用 frollapply：组长度 < 窗口时它返回 window-1 个
#         元素（实测 n=54 -> 返回 59），赋值长度不匹配 → 全量数据必崩
#         （27 只股票历史 < 60 天）。改用 frollmean 恒等式 E[r²]-E[r]²（兼收性能）。
#   [H-2] 行业哑变量缺失：industry_at() 已写好但从未接入 B → 已接入。
#   [A]   代码前缀不一致：daily 带 sh/sz/bj 前缀、industry 不带 → 匹配率 0.00%。
#         新增 norm_code() 统一口径 → 5162/5162 = 100%。
#   [B]   industry_at() 用 start<=d<=end，但 end 是"下一次变更日"不是"有效期至"
#         （只有 9.6% 的股票 last_end >= 面板末日）→ 匹配率掉到 6.9%~56.3%。
#         改用 LOFF（取 start<=d 中最晚一条）→ 99.9%~100%。
#   [M-2] σ² 两套实现并存 → 统一为 EWMA(0.94)，删除 est_specific_var。
#   [M-3] cs_wls 原先没有正则项，而 cs_ridge 有 λ → 两者不是"只差加权"。
#         现两处都加【尺度相对】稳定项 stab。
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

# =============================================================================
# 0. 全局参数（取值依据写在旁边，改之前务必读一遍）
# =============================================================================

# ---- EWMA 衰减因子 λ --------------------------------------------------------
# EWMA 递归： s_t = (1-λ)·x_{t-1} + λ·s_{t-1}；半衰期 h = -ln2/lnλ
#
# λ_SPECIFIC = 0.94（半衰期 ≈ 11.3 天）—— 用于个股特异方差 σ²
#   【为什么是 0.94】
#   (1) 模拟标定：在 GARCH(1,1)（α=0.08, β=0.90，日频股票典型量级）上，用 QLIKE
#       损失做的 OOS 最优 λ* ≈ 0.9275（半衰期 9.2 天）。0.94 只比最优差 0.0008；
#       而"半衰期 60 天"（λ=0.9885）比最优差 0.0550，约 70 倍。
#       ⇒ 传统拍脑袋的 60 天半衰期明显过慢，必须改。
#   (2) 为什么不直接上 MLE：QLIKE 曲面在最优附近极平（λ=0.97 也只差 0.017），
#       "大致对"就够；而 in-sample MLE 在三种情形下会系统性偏向【过快】衰减——
#       无波动聚集时（真值 λ→1 却选出 0.91）、厚尾 t(ν=3) 时（选 0.96 vs 最优 0.995）、
#       样本短时（T=250 直接撞网格下界 0.80）。滚动重估虽原理干净，但实测
#       OOS 损失与固定 0.94 打平。⇒ 用稳健的行业标准值更划算。
#   (3) 0.94 即 RiskMetrics 的经典日频取值，可对外直接引用。
LAMBDA_SPECIFIC <- 0.94

# λ_FACTOR = 0.97（半衰期 ≈ 22.8 天）—— 用于因子协方差矩阵 F
#   【为什么比 σ² 慢】
#   F 是 K×K 的"共同风险"矩阵，其动态比个股特异波动慢一个量级；且 λ 越小
#   （衰减越快）估计越吵，K 越大越明显（接入行业后 K≈35）。σ² 与 F 的
#   最优衰减速度本来就不同，故分开设两个常数。
LAMBDA_FACTOR <- 0.97

# EWMA 预热期：0.94^250 ≈ 2e-7，初值影响可忽略。面板需比演示窗口多回溯这么久。
EWMA_BURN <- 250L

# ---- 正则稳定项（尺度相对）--------------------------------------------------
# 【为什么不用裸 λ】横截面回归 N≫K（约 5000 vs 35），B'B 条件数极好，
#   数值上并不需要 ridge。文档《横截面风险估计》§111/§253 只要一个"小 λ
#   做数值稳定"（目的是让"唯一差别 = 加权"更纯粹）。
#   但裸 λ（如原先 cs_ridge 的 1e-3）【没有尺度归一化】：B 是横截面 z-score，
#   tr(B'B)/K ≈ N ≈ 5000，于是 1e-3/5000 = 2e-7 —— 实际等于 0，且换个尺度
#   （比如 B 不 z-score）含义就完全变了。
#   ⇒ 改用尺度相对形式： pen = δ · tr(A)/K，其中 A = B'B（或 B'Ω⁻¹B）。
#      这样 δ 无量纲、与 N 和暴露尺度无关，ridge 与 WLS 口径也一致。
STAB_DELTA <- 1e-6

# =============================================================================
# 0b. 工具函数
# =============================================================================

find_col <- function(df, ...) {
  keys <- c(...)
  for (cn in names(df)) if (any(vapply(keys, grepl, logical(1), x = cn))) return(cn)
  NA_character_
}

cs_zscore <- function(x) {
  x <- as.numeric(x); s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' 统一股票代码口径：去掉 sh/sz/bj 前缀
#' 【为什么需要】daily.parquet 的 code 带前缀（sh600000），industry.parquet 不带
#'   （600000）。不统一则行业 PIT 连接匹配率为 0.00%（实测）。去前缀后 100%。
norm_code <- function(x) sub("^(sh|sz|bj)", "", as.character(x))

#' 滞后一期的 EWMA（防前视）：s_t = (1-λ)·x_{t-1} + λ·s_{t-1}
#' 用 stats::filter 的递归模式，O(n)，比 frollapply 快一个量级。
ewma_lag1 <- function(x, lambda) {
  n <- length(x)
  if (n < 2L) return(rep(NA_real_, n))
  u <- c(0, (1 - lambda) * x[-n])
  s <- as.numeric(stats::filter(u, lambda, method = "recursive"))
  s[1] <- NA_real_        # t=1 无历史 → NA（filter 的递归从 0 起步，需显式置 NA）
  s
}

# =====================================================================
# A. 核心估计量（纯函数，重点单测）
# =====================================================================

#' Stage2-baseline：ridge 横截面回归  f = (B'B + pen·I)^{-1} B'r
#'   pen = lambda + stab · tr(B'B)/K   （尺度相对稳定项，见文件头 STAB_DELTA 说明）
cs_ridge <- function(r, B, lambda = 0, stab = 0) {
  B <- as.matrix(B); r <- as.numeric(r)
  p <- ncol(B)
  A <- crossprod(B)
  if (stab > 0) A <- A + (stab * sum(diag(A)) / p) * diag(p)
  if (lambda != 0) A <- A + lambda * diag(p)
  f <- tryCatch(as.numeric(solve(A, crossprod(B, r))),
                error = function(e) rep(NA_real_, p))
  names(f) <- colnames(B)
  f
}

#' Stage2-★优化★：WLS（异方差有效）  f = (B'Ω⁻¹B)^{-1} B'Ω⁻¹r, Ω=diag(σ²)
#'   stab 与 cs_ridge 同口径 —— 保证两版"唯一差别 = 加权"，而不是"ridge 有正则、
#'   WLS 没有"。pen 相对 tr(B'Ω⁻¹B)/K 归一，与 σ² 的量级无关。
cs_wls <- function(r, B, sigma2, stab = 0) {
  B <- as.matrix(B); r <- as.numeric(r); sigma2 <- as.numeric(sigma2)
  if (length(sigma2) != nrow(B)) stop("sigma2 长度须等于股票数")
  if (any(!is.finite(sigma2)) || any(sigma2 <= 0)) stop("sigma2 须为正且有限")
  w <- 1 / sigma2                       # Ω⁻¹ 对角
  Bw <- B * w                           # 每行乘 w_i（R 按列回收 -> 行缩放）
  A <- crossprod(B, Bw)                 # B'Ω⁻¹B
  p <- ncol(B)
  if (stab > 0) A <- A + (stab * sum(diag(A)) / p) * diag(p)
  b <- crossprod(B, r * w)              # B'Ω⁻¹r
  f <- tryCatch(as.numeric(solve(A, b)), error = function(e) rep(NA_real_, p))
  names(f) <- colnames(B)
  f
}

#' EWMA 协方差（对 T×K 矩阵，行 = 时间，列 = 变量）
#' 默认 λ = LAMBDA_FACTOR（因子协方差比特异方差衰减更慢，见文件头）
#' 注：Σw=1 得到的是加权二阶矩（未做 1/(1-Σw²) 的无偏修正）；对风险模型可接受。
ewma_cov <- function(X, lambda = LAMBDA_FACTOR) {
  X <- as.matrix(X)
  ok <- rowSums(!is.finite(X)) == 0L
  X <- X[ok, , drop = FALSE]
  Tn <- nrow(X); K <- ncol(X)
  if (Tn < K + 2L) return(matrix(NA_real_, K, K, dimnames = list(colnames(X), colnames(X))))
  wgt <- lambda^((Tn - 1L):0L); wgt <- wgt / sum(wgt)   # 最近一期权重最大
  mu <- colSums(X * wgt)                                # 加权均值
  Xc <- sweep(X, 2, mu, "-")
  crossprod(Xc * sqrt(wgt))                             # Σ_t w_t x_t x_t'
}

#' Stage3：主动风险矩阵  Σ = B F B' + D   （⚠️ N×N，仅供小 N / 单测使用）
#' 生产环境算组合风险请用 pred_risk()（因子形式 O(NK)，不构造 N×N）
risk_matrix <- function(B, F, D) {
  B <- as.matrix(B)
  S <- B %*% as.matrix(F) %*% t(B)
  diag(S) <- diag(S) + as.numeric(D)
  S
}

#' 组合主动风险  sqrt(a'Σa)，a = w - w_bench   （⚠️ 同 risk_matrix，N×N）
active_risk <- function(w, w_bench, Sigma) {
  a <- as.numeric(w) - as.numeric(w_bench)
  sqrt(as.numeric(crossprod(a, as.matrix(Sigma) %*% a)))
}

#' ★ 预测组合风险：sqrt( w'Σw )，用因子形式避免构造 N×N
#'   w'Σw = (B'w)'F(B'w) + Σ_i w_i² D_i
pred_risk <- function(w, B, F, D) {
  x <- as.numeric(crossprod(as.matrix(B), as.numeric(w)))   # K 维
  as.numeric(sqrt(crossprod(x, as.matrix(F) %*% x) + sum(as.numeric(w)^2 * as.numeric(D))))
}

#' ★ 个股预测风险（Σ 的对角），因子形式 O(NK)：sqrt((BF)_i·B_i + D_i)
diag_risk <- function(B, F, D) {
  BF <- as.matrix(B) %*% as.matrix(F)
  sqrt(rowSums(BF * as.matrix(B)) + as.numeric(D))
}

# =====================================================================
# B. 暴露矩阵
# =====================================================================

#' PIT 行业查询：给定 (code,date) 返回当时行业 id
#' ind 需含列：股票代码/计入日期/变更日期/行业id（或已改名 code/start/end/industry）
#'
#' 【口径】LOFF（last observation carried forward）：取 start <= date 中【最晚】的一条。
#'   不能要求 start<=date<=end —— 实测 end 是"下一次变更日"而非"有效期至"，
#'   只有 9.6% 的股票 last_end >= 面板末日，硬卡 end 会把匹配率从 100% 压到 6.9%~56.3%。
#'   另外 30.8% 的股票存在区间重叠，非等值连接会返回重复行 → LOFF 天然规避。
#' 返回的 code 是【调用方传入的原样 code】，便于与 daily 面板合并。
industry_at <- function(ind, code, date) {
  ind2 <- as.data.table(ind)
  if (!"code" %in% names(ind2))
    setnames(ind2, c(find_col(ind2, "股票代码"), find_col(ind2, "计入日期"),
                     find_col(ind2, "变更日期"), find_col(ind2, "行业id")),
             c("code", "start", "end", "industry"), skip_absent = TRUE)
  ind2 <- ind2[, .(code = norm_code(code), start = as.Date(start),
                   industry = as.character(industry))]
  ind2 <- ind2[!is.na(code) & !is.na(start) & !is.na(industry)]
  ind2 <- unique(ind2, by = c("code", "start"))   # 同 (code,start) 只留一条
  setorder(ind2, code, start)
  # code_raw 保留调用方传入的【原始】代码（可能带 sh/sz/bj 前缀），
  # 因为下游要与 daily 面板（带前缀）合并；norm_code 只用于连接。
  q <- data.table(code = norm_code(code), date = as.Date(date),
                  code_raw = as.character(code))
  # ⚠️ 必须用【等值】连接 + roll=TRUE 实现 LOFF：本版 data.table 的 roll
  #    不支持非等值连接（on = .(code, start <= date) 会报 "roll is not implemented
  #    for non-equi joins yet"）。on=(code, start=date) + roll=TRUE 语义完全相同：
  #    取 start <= date 中最大的一条（last observation carried forward）。
  ind2[q, on = .(code, start = date), roll = TRUE,
       .(code = i.code_raw, date = i.date, industry = x.industry), nomatch = NA]
}

#' 风格暴露（价格类）：size / momentum / volatility / liquidity
#' 返回长表 date, code, exposure, value（已横截面标准化）
#' ⚠️ 行业哑变量【不】走这里 —— 它必须保持 0/1，不能被 cs_zscore 标准化。
#'    行业行由 build_industry_exposure() 单独产出后再拼接。
build_style_exposures <- function(daily, lookback = 250L, skip = 21L, vol_win = 60L) {
  dt <- copy(daily)[, .(date, code, close, amount, turnover)]
  setorder(dt, code, date)
  dt[, ret := close / shift(close) - 1, by = code]
  dt[, size := fifelse(turnover > 0, log(pmax(amount / turnover, 1e-6)), NA_real_)]
  dt[, liquidity := frollmean(turnover, vol_win, align = "right", na.rm = TRUE), by = code]
  # [H-1 修复] 原为 frollapply(ret, vol_win, sd, ...)：组长度 < 窗口时 frollapply
  # 返回 window-1 个元素，与组长度不匹配 → 全量数据崩溃。改用 frollmean 恒等式
  # E[r²]-E[r]²（数学上与 sd² 等价，单测已验证相对差 < 5e-16），并且快一个量级。
  # [M-2 修复] Bessel 修正必须用【窗口内实际有效观测数 k】而不是名义窗口 vol_win：
  # na.rm=TRUE 时窗口部分缺失也会给值，用 vol_win 会低估 sqrt((k-1)/k) 倍。
  # 同时用 (n_ok == vol_win) 强制"完整窗口才出值"，与 frollapply 的行为一致。
  dt[, m1 := frollmean(ret, vol_win, align = "right", na.rm = TRUE), by = code]
  dt[, m2 := frollmean(ret^2, vol_win, align = "right", na.rm = TRUE), by = code]
  dt[, n_ok := frollsum(as.integer(!is.na(ret)), vol_win, align = "right"), by = code]
  dt[, volatility := fifelse(n_ok == vol_win,
                             sqrt(pmax(m2 - m1^2, 0) * n_ok / (n_ok - 1) * 252),
                             NA_real_), by = code]
  dt[, momentum := fifelse(.I > lookback + skip,
                           shift(close, skip) / shift(close, skip + lookback) - 1, NA_real_),
     by = code]
  long <- melt(dt[, .(date, code, size, momentum, volatility, liquidity)],
               id.vars = c("date", "code"), variable.name = "exposure", value.name = "value")
  long[, value := cs_zscore(value), by = .(date, exposure)]
  long[]
}

#' ★ 行业哑变量（长表）：每个 (date, code) 一行，exposure = "IND<行业id>"，value = 1
#' 与 build_style_exposures 的输出 schema 完全一致，rbind 后交给 dcast 即成 0/1 矩阵。
#' 【编码选择】用**原始 0/1 哑变量 + 无截距**（而非去均值版）：
#'   数学上"原始哑变量（无截距）"与"去均值哑变量 + 截距"张成同一列空间
#'   ⇒ 拟合值、残差、D、F、Σ 完全相同，只是 f 的读法不同
#'   （原始：行业绝对收益；去均值：行业相对市场超额）。选原始版是因为
#'   cs_ridge/cs_wls 没有截距项，加原始哑变量正好满秩。附带好处：风格暴露被
#'   哑变量部分掉之后，风格因子收益自动变成【行业内】的收益（标准读法）。
build_industry_exposure <- function(daily, ind, dates = NULL) {
  d <- unique(copy(daily)[, .(date, code)])
  if (!is.null(dates)) d <- d[date %in% dates]
  q <- industry_at(ind, d$code, d$date)
  out <- q[!is.na(industry), .(date, code, exposure = paste0("IND", industry), value = 1)]
  out[]
}

#' ★ 收益面板（一次性算好，供 make_cs 查表）
#' 【M-3】替代 make_cs 里"每天 copy(daily)[date<=d] + .SD[.N] 扫全历史"的 O(历史/天) 做法。
#'   语义完全等价：shift(close) 在按 (code, date) 排序的序列上取到的就是"上一交易日收盘价"，
#'   与 dt[date < d][, .SD[.N], by = code] 相同（前提：同股票同日期不重复）。
#'   复杂度：O(N) 一次 + O(1) 查表（原为 O(n_days × 历史长度)）。
build_ret_panel <- function(daily) {
  d <- copy(daily)[, .(date, code, close)]
  setorder(d, code, date)
  d[, c_p := shift(close), by = code]
  d <- d[is.finite(c_p) & is.finite(close)]
  d[, .(date, code, ret = close / c_p - 1)]
}

#' 拼装横截面：把长表暴露转成 (date -> list(r, B, codes))
#  （数据量大，按需逐日构建，供评价环节使用）
make_cs_list <- function(daily, style_long, ind_at = NULL, dates = NULL) {
  r <- copy(daily)[, .(date, code, close)]
  setorder(r, code, date)
  r[, ret := close / shift(close) - 1, by = code]
  r <- r[is.finite(ret), .(date, code, ret)]
  if (!is.null(dates)) r <- r[date %in% dates]
  r
}

# =====================================================================
# C. 评价
# =====================================================================

#' L1：用 t 日估计的因子收益预测 t+1 横截面收益，算 OOS RMSE
#' f_by_date: data.table(date, term, f)  —— 由 cs_ridge / cs_wls 输出整理
l1_pred_rmse <- function(cs, f_by_date) {
  dates <- sort(unique(cs$ret$date))
  errs <- c()
  for (i in seq_len(length(dates) - 1L)) {
    t0 <- dates[i]; t1 <- dates[i + 1L]
    f0 <- f_by_date[date == t0]
    if (!nrow(f0)) next
    B0 <- cs$B[[as.character(t0)]]; r1 <- cs$ret[date == t1]
    if (is.null(B0)) next
    if (!all(f0$term %in% colnames(B0))) next
    pred <- as.numeric(B0[, f0$term, drop = FALSE] %*% f0$f)
    act <- r1$ret[match(cs$codes[[as.character(t0)]], r1$code)]
    ok <- is.finite(pred) & is.finite(act)
    if (any(ok)) errs <- c(errs, (pred[ok] - act[ok])^2)
  }
  if (!length(errs)) return(NA_real_)
  sqrt(mean(errs))
}

#' ★ L2：Bias Statistic（Barra 风格）
#' z = 实际收益 / 预测风险。模型正确时 E[z]≈0 且 sd(z)≈1。
#'   sd(z) > 1 → 【低估】风险（缺共同因子时必然如此）
#'   sd(z) < 1 → 【高估】风险
#' t_sd：检验 sd(z)=1 的近似 z 值，se(sd) ≈ 1/sqrt(2(n-1))（正态近似）
#'   |t_sd| > 2 ⇒ 显著偏离 1
bias_statistic <- function(z) {
  z <- z[is.finite(z)]
  n <- length(z)
  if (n < 10L) return(data.table(n = n, mean_z = NA_real_, sd_z = NA_real_, t_sd = NA_real_))
  m <- mean(z); s <- sd(z)
  data.table(n = n, mean_z = m, sd_z = s, t_sd = (s - 1) * sqrt(2 * (n - 1L)))
}

#' 生成测试组合权重矩阵（N 行 = 股票，M 列 = 组合）
#' method: "random" 随机多空（L1 归一，总杠杆 1）/ "equal" 等权多头
make_test_portfolios <- function(N, M = 100L, method = "random", seed = 1L) {
  if (method == "equal") return(matrix(1 / N, N, 1L))
  set.seed(seed)
  W <- matrix(stats::runif(N * M) - 0.5, N, M)
  lev <- colSums(abs(W))
  sweep(W, 2L, ifelse(lev > 0, lev, 1), "/")
}

# =====================================================================
# D. 主流程（规范两遍法：ridge 取残差 -> σ² 面板 -> WLS -> Σ -> 三层评价）
# =====================================================================

#' 取某日横截面：list(code, B, r, n_drop)
#'
#' 【M-3 性能修复】原实现每天 copy(daily)[date <= d] + dt[date < d][, .SD[.N], by=code]
#'   扫全历史（O(历史长度)/天，实测 2.63 s/日 → 全量 25 分钟）。
#'   改为消费预计算的 ret_panel（build_ret_panel 一次性算好），每天只做 O(当日股票数) 的查表。
#'   前置条件：ret_panel / style / ind_long 都已按 date 建 key（run_wls_module 里做）。
#'
#' 【M-6 健壮性修复】原来 B[!is.finite(B)] <- 0 会把缺失暴露【静默】插补成横截面均值。
#'   现改为【整行删除】(listwise deletion)：任一暴露缺失就把该股票当日剔除，并回报剔除数量。
#'   剔除后若样本不足 min_n，或出现退化列（常数列 → B'B 奇异），返回 NULL 让上层跳过该日。
make_cs <- function(ret_panel, style, d, ind_long = NULL, min_n = 50L) {
  Bw <- dcast(style[.(d)], code ~ exposure, value.var = "value")
  n_drop <- 0L

  if (!is.null(ind_long)) {
    iw <- dcast(ind_long[.(d)], code ~ exposure, value.var = "value")
    ind_cols <- setdiff(names(iw), "code")
    # ⚠️ 关键区分：行业哑变量的 NA 是【结构性 0】（"该股票不属于这个行业"），
    #    不是缺失数据 —— dcast 对未出现的 (code, exposure) 组合一律填 NA。
    #    而风格暴露的 NA 是【真缺失】（如动量需要 271 天历史）。
    #    前者填 0（正确），后者整行剔除（M-6 的"舍弃"）。
    for (cc in ind_cols) set(iw, which(is.na(iw[[cc]])), cc, 0)
    Bw <- merge(Bw, iw, by = "code", all.x = TRUE)
    for (cc in ind_cols) set(Bw, which(is.na(Bw[[cc]])), cc, 0)
    # 行业完全未知（所有哑变量为 0）的股票剔除 —— 塞进"无行业"桶会污染行业因子收益
    has_ind <- rowSums(as.matrix(Bw[, ..ind_cols])) > 0
    n_drop <- n_drop + sum(!has_ind)
    Bw <- Bw[has_ind]
  }
  if (!nrow(Bw)) return(NULL)

  m <- merge(Bw, ret_panel[.(d)], by = "code")
  m <- m[is.finite(ret)]
  if (!nrow(m)) return(NULL)

  # 注意：cols 必须在【行业 merge 之后】取，否则只会拿到 4 个风格列
  cols <- setdiff(names(Bw), "code")
  B <- as.matrix(m[, ..cols])
  r <- m$ret; code <- m$code

  # ---- M-6：风格暴露缺失 → 整行剔除（不再静默填 0 = 横截面均值）----
  # 行业哑变量的 NA 已在上面按"结构性 0"处理掉，所以这里剩下的 NA 只有真缺失。
  ok <- rowSums(!is.finite(B)) == 0L
  n_drop <- n_drop + sum(!ok)
  if (any(!ok)) { B <- B[ok, , drop = FALSE]; r <- r[ok]; code <- code[ok] }
  if (length(r) < min_n) return(NULL)

  # ---- 退化列检测：常数列使 B'B 奇异 ----
  const_col <- vapply(seq_len(ncol(B)),
                      function(j) all(B[, j] == B[1L, j]), logical(1))
  if (any(const_col)) return(NULL)

  list(code = code, B = B, r = r, n_drop = n_drop)
}

run_wls_module <- function(data_dir = "data", lambda_ridge = 0,
                           stab = STAB_DELTA, lambda_spec = LAMBDA_SPECIFIC,
                           lambda_F = LAMBDA_FACTOR,
                           n_days = 250, max_stocks = NULL,
                           M_port = 100L, seed = 1L) {
  cat("== 横截面模块（风险线，优化点 B）==\n")
  daily <- as.data.table(arrow::read_parquet(file.path(data_dir, "daily.parquet")))
  daily[, date := as.Date(date)]
  if (!is.null(max_stocks)) {            # 抽样：便于冒烟测试 / 快速迭代
    kc <- head(sort(unique(daily$code)), max_stocks)
    daily <- daily[code %in% kc]
    cat(sprintf("  [抽样] 仅用前 %d 只股票\n", max_stocks))
  }
  cat(sprintf("  daily: %d 行 / %d 只 / %s ~ %s\n", nrow(daily), uniqueN(daily$code),
              min(daily$date), max(daily$date)))

  # ---- 暴露矩阵：风格 ----
  style <- build_style_exposures(daily)
  all_days <- sort(unique(style$date))
  # 面板 = 演示窗口 + EWMA 预热期（预热期内 EWMA 初值影响未衰减完，不参与评价）
  panel_days <- tail(all_days, min(length(all_days), n_days + EWMA_BURN))
  demo_days  <- tail(panel_days, min(n_days, length(panel_days)))
  # 只用面板日的暴露：build_style_exposures 对全期 2430 天产出约 3800 万行，
  # 而下游 make_cs 只按面板日切片 → 提前过滤可省约 1 GB 内存并加速逐日筛选。
  style <- style[date %in% panel_days]
  cat(sprintf("  面板 %s~%s (%d 天) | 演示 %d 天 | EWMA 预热 %d 天\n",
              min(panel_days), max(panel_days), length(panel_days),
              length(demo_days), EWMA_BURN))

  # ---- 暴露矩阵：行业（只对面板日构建，省内存 + 加速 make_cs 的逐日筛选）----
  ind_long <- NULL
  ip <- file.path(data_dir, "industry.parquet")
  if (file.exists(ip)) {
    ind <- as.data.table(arrow::read_parquet(ip))
    ind_long <- build_industry_exposure(daily, ind, dates = panel_days)
    # 覆盖率 = 当日匹配到行业的股票数 / 当日有行情的股票数（逐日平均）
    cov_dt <- merge(ind_long[, .(n_ind = uniqueN(code)), by = date],
                    daily[date %in% panel_days, .(n_all = uniqueN(code)), by = date],
                    by = "date")
    cat(sprintf("  行业暴露: %d 行 | 覆盖率 %.1f%%（逐日均值）| 行业数 %d | K = %d\n",
                nrow(ind_long), 100 * cov_dt[, mean(n_ind / n_all)],
                uniqueN(sub("^IND", "", ind_long$exposure)),
                4L + uniqueN(sub("^IND", "", ind_long$exposure))))
  } else {
    cat("  ⚠️ 找不到 industry.parquet → B 只含风格暴露（行业风险未建模）\n")
  }

  # ---- 收益面板 + 建 key（M-3：一次性填表，之后逐日 O(log n) 查表）----
  t_pre <- Sys.time()
  ret_panel <- build_ret_panel(daily)
  setkey(ret_panel, date); setkey(style, date)
  if (!is.null(ind_long)) setkey(ind_long, date)
  cat(sprintf("  收益面板: %d 行（预计算 %.1f s）\n", nrow(ret_panel),
              as.numeric(difftime(Sys.time(), t_pre, units = "secs"))))

  t_cs <- Sys.time()
  cs <- lapply(panel_days, function(d) make_cs(ret_panel, style, d, ind_long))
  empty <- vapply(cs, is.null, logical(1))
  n_drop_tot <- sum(vapply(cs, function(x) if (is.null(x)) 0L else x$n_drop, 0L))
  cs <- cs[!empty]; panel_days <- panel_days[!empty]
  cat(sprintf("  可用横截面: %d 天（跳过 %d 天）| B 维度: %d 列 | 缺失暴露剔除 %d 行 | 逐日切片 %.1f s\n",
              length(cs), sum(empty), ncol(cs[[1]]$B), n_drop_tot,
              as.numeric(difftime(Sys.time(), t_cs, units = "secs"))))

  # ---- Pass 1: ridge（baseline）----
  f_ridge <- lapply(cs, function(x) cs_ridge(x$r, x$B, lambda_ridge, stab))
  resid_dt <- rbindlist(lapply(seq_along(cs), function(i) {
    data.table(date = panel_days[i], code = cs[[i]]$code,
               resid = cs[[i]]$r - as.numeric(cs[[i]]$B %*% f_ridge[[i]]))
  }))

  # ---- Stage 1: σ² 面板（EWMA(λ_spec)，滞后一期，防前视）----
  # 原实现为 frollmean(shift(resid²,1), win) 的【等权窗】；现统一为 EWMA。
  setorder(resid_dt, code, date)
  resid_dt[, sig2 := ewma_lag1(resid^2, lambda_spec), by = code]
  # 前 EWMA_BURN 天初值影响未衰减完 → 用当日横截面中位数兜底
  resid_dt[is.na(sig2) | seq_len(.N) <= EWMA_BURN, sig2 := NA_real_, by = code]
  resid_dt[is.na(sig2), sig2 := median(resid^2, na.rm = TRUE), by = date]
  resid_dt[, sig2 := { fl <- quantile(sig2, 0.05, na.rm = TRUE); pmax(sig2, fl) }, by = date]

  # ---- Pass 2: WLS（★唯一改动）—— 只对演示窗口 ----
  demo_idx <- which(panel_days %in% demo_days)
  f_wls <- vector("list", length(cs))
  for (i in demo_idx) {
    s <- resid_dt[date == panel_days[i], .(code, sig2)]
    s2 <- s$sig2[match(cs[[i]]$code, s$code)]
    bad <- !is.finite(s2) | s2 <= 0
    if (all(bad)) { f_wls[[i]] <- rep(NA_real_, ncol(cs[[i]]$B)); next }
    if (any(bad)) s2[bad] <- median(s2[!bad])
    f_wls[[i]] <- cs_wls(cs[[i]]$r, cs[[i]]$B, s2, stab)
  }

  comp <- rbindlist(lapply(demo_idx, function(i)
    data.table(date = panel_days[i], term = names(f_ridge[[i]]),
               f_ridge = f_ridge[[i]], f_wls = f_wls[[i]])))
  cat("\n-- 各暴露的因子收益：ridge vs WLS（仅风格暴露，行业略）--\n")
  sty <- comp[!grepl("^IND", term)]
  print(sty[, .(sd_ridge = sd(f_ridge, na.rm = TRUE), sd_wls = sd(f_wls, na.rm = TRUE),
                cor = cor(f_ridge, f_wls, use = "complete.obs")), by = term])

  # ---- Stage 3: F / D / Σ（两版各算一套，方法完全相同 → 差别只来自 f）----
  mkF <- function(fl) {
    M <- do.call(rbind, lapply(demo_idx, function(i) fl[[i]]))
    colnames(M) <- names(fl[[demo_idx[1]]])
    ewma_cov(M, lambda_F)
  }
  F_ridge <- mkF(f_ridge); F_wls <- mkF(f_wls)

  D_of <- function(i) {
    s <- resid_dt[date == panel_days[i], .(code, sig2)]
    s$sig2[match(cs[[i]]$code, s$code)]
  }

  # ---- Stage 4-L1：f 预测下一日横截面收益的 OOS RMSE ----
  l1 <- rbindlist(lapply(seq_len(length(demo_idx) - 1L), function(k) {
    i <- demo_idx[k]; j <- demo_idx[k + 1L]
    if (j != i + 1L) return(NULL)
    nxt <- cs[[j]]; idx <- match(cs[[i]]$code, nxt$code); ok <- !is.na(idx)
    if (!any(ok)) return(NULL)
    act <- nxt$r[idx[ok]]
    pr <- as.numeric(cs[[i]]$B[ok, , drop = FALSE] %*% f_ridge[[i]])
    pw <- as.numeric(cs[[i]]$B[ok, , drop = FALSE] %*% f_wls[[i]])
    data.table(ridge = (pr - act)^2, wls = (pw - act)^2)
  }), fill = TRUE)
  if (nrow(l1)) {
    rmse_r <- sqrt(mean(l1$ridge, na.rm = TRUE)); rmse_w <- sqrt(mean(l1$wls, na.rm = TRUE))
    cat(sprintf("\n★ L1 (f 预测下一日横截面收益) OOS RMSE:\n    ridge = %.6g   WLS = %.6g   %s\n",
                rmse_r, rmse_w, ifelse(rmse_w < rmse_r, "WLS 更优 ✅", "WLS 未占优 ⚠️")))
  }

  # ---- Stage 4-L2：Σ 风险精度（Bias Statistic）ridge vs WLS ----
  # 目的不是"Σ 准不准"，而是"WLS 版 Σ 是否比 ridge 版更准"。
  cat("\n★ L2 (Bias Statistic：z = 实际收益/预测风险，理想 sd(z)=1)\n")
  l2 <- rbindlist(lapply(seq_len(length(demo_idx) - 1L), function(k) {
    i <- demo_idx[k]; j <- demo_idx[k + 1L]
    if (j != i + 1L) return(NULL)
    nxt <- cs[[j]]; idx <- match(cs[[i]]$code, nxt$code); ok <- !is.na(idx)
    if (sum(ok) < 10L) return(NULL)
    act <- nxt$r[idx[ok]]
    B <- cs[[i]]$B[ok, , drop = FALSE]
    Dr <- D_of(i)[ok]; Dw <- Dr
    sr <- diag_risk(B, F_ridge, Dr); sw <- diag_risk(B, F_wls, Dw)
    Wp <- make_test_portfolios(sum(ok), M_port, "random", seed + k)
    pzr <- as.numeric(crossprod(Wp, act)) / apply(Wp, 2, function(w)
      pred_risk(w, B, F_ridge, Dr))
    pzw <- as.numeric(crossprod(Wp, act)) / apply(Wp, 2, function(w)
      pred_risk(w, B, F_wls, Dw))
    data.table(kind = c(rep("个股", sum(ok)), rep("随机组合", ncol(Wp))),
               ridge = c(act / sr, pzr), wls = c(act / sw, pzw))
  }), fill = TRUE)
  if (nrow(l2)) {
    L2 <- l2[, rbind(bias_statistic(ridge), bias_statistic(wls)), by = kind]
    L2[, method := rep(c("ridge", "wls"), length.out = .N)]
    print(L2[, .(kind, method, n, mean_z = round(mean_z, 4), sd_z = round(sd_z, 4),
                 t_sd = round(t_sd, 2))])
    cat("  读法：sd_z 越接近 1 越好；明显 >1 = 低估风险（缺共同因子时常见）\n")
  } else {
    cat("  [L2] 样本不足\n")
  }

  invisible(list(style = style, f = comp, F_ridge = F_ridge, F_wls = F_wls,
                 sigma2 = resid_dt, l1 = l1, l2 = l2))
}

if (!nzchar(Sys.getenv("DFZQ_SOURCE_ONLY"))) {
  # 直接 `Rscript _crosssec_wls.R` 会用【全量 5162 只 × 500 天面板】，
  # 按 make_cs 约 2.6 s/日估算 ≈ 25~35 分钟，且 cs 列表约占 700 MB 内存。
  # 想快速试跑，用环境变量控制（无需改代码）：
  #   PowerShell : $env:DFZQ_MAX_STOCKS="300"; $env:DFZQ_N_DAYS="120"; Rscript _crosssec_wls.R
  #   Bash       : DFZQ_MAX_STOCKS=300 DFZQ_N_DAYS=120 Rscript _crosssec_wls.R
  #   两者都不设 = 全量（默认）
  .ms <- Sys.getenv("DFZQ_MAX_STOCKS", "")
  .nd <- Sys.getenv("DFZQ_N_DAYS", "")
  out <- run_wls_module(
    max_stocks = if (nzchar(.ms)) as.integer(.ms) else NULL,
    n_days     = if (nzchar(.nd)) as.integer(.nd) else 250L
  )
  cat("\n完成。\n")
}
