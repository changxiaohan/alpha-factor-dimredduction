# =============================================================================
# _factor_ext.R -- 因子扩充（信号线 P0）：把 2 个价格因子扩到 10 个
# -----------------------------------------------------------------------------
# 背景：此前所有实验用的信号是 z(momentum) - z(volatility)，只有 2 个因子。
#   本脚本用【财务数据】补上 价值 / 质量 / 成长 / 杠杆 四类因子，
#   并【修正 PIT 口径】（代码审查发现原 +90 天规则对年报前视 31 天）。
#
# 【PIT 口径修正】
#   原：报告期 + 90 天（对年报 = 次年 3-31，但法定披露截止是 4-30 -> 前视 31 天）
#   现：按报告期类型用【法定披露截止】作可用日：
#       0331 -> +30 天(4-30) | 0630 -> +62 天(8-31) | 0930 -> +31 天(10-31) | 1231 -> +120 天(次年4-30)
#   ⚠️ 取法定截止日是【保守】做法（实际披露可能更早），保证不产生前视。
#
# 【因子清单】
#   价格类（已有口径）：momentum, volatility, size, liquidity
#   新 增：value_bp  (每股净资产/不复权价)      —— 账面市值比
#          value_ep  (基本每股收益/不复权价)    —— 盈利收益率
#          quality_roe (净资产收益率)
#          quality_gpm (毛利率)
#          growth_np   (归母净利润增长率)
#          leverage    (资产负债率，取负号 -> 低杠杆更好)
#
# ⚠️ 估值因子用 amount/turnover 作【不复权价】——close 是前复权价，不能配 BPS/EPS 使用
#
# 用法（约 3~5 分钟）：
#   & "...\Rscript.exe" "_factor_ext.R" | Tee-Object "_factor_ext.log"
# 输出：_factor_panel.csv（date, code, factor, value_z）
# 输出全 ASCII。
# =============================================================================
Sys.setenv(DFZQ_SOURCE_ONLY = "1")
suppressPackageStartupMessages({ library(data.table); library(arrow) })
source("_crosssec_wls.R")

P  <- function(...) cat(sprintf(...), "\n")
hr <- function(t) cat("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78), "\n", sep = "")
el <- function() as.numeric(difftime(Sys.time(), T0, units = "secs"))

T0 <- Sys.time()
NP <- as.integer(Sys.getenv("DFZQ_PANEL_DAYS", "0")); if (is.na(NP)) NP <- 0L

## ================= 0. 载入 =================
daily <- as.data.table(arrow::read_parquet("data/daily.parquet")); daily[, date := as.Date(date)]
daily[, code := norm_code(code)]                      # 统一代码口径（去 sh/sz/bj 前缀）
all_days <- sort(unique(daily$date))
panel_days <- if (NP > 0L) tail(all_days, NP) else all_days
P("[0] daily: %d rows | %d codes | %d panel days (%s ~ %s)",
  nrow(daily), uniqueN(daily$code), length(panel_days),
  min(panel_days), max(panel_days))

## ================= 1. 财务长表 + PIT 可用日（修正 +90 bug）=================
fin <- as.data.table(arrow::read_parquet("data/financial.parquet"))
fin[, code := norm_code(code)]
IND <- c(roe   = "净资产收益率(ROE)",
         gpm   = "毛利率",
         npr   = "销售净利率",
         bps   = "每股净资产",
         eps   = "基本每股收益",
         npg   = "归属母公司净利润增长率",
         rg    = "营业总收入增长率",
         lev   = "资产负债率",
         eq    = "股东权益合计(净资产)",
         rev   = "营业总收入",
         np    = "归母净利润")
per_cols <- names(fin)[grepl("^\\d{8}$", names(fin))]
per_cols <- sort(per_cols[per_cols >= "20131231"])     # 面板 2016 起，多留一年缓冲
P("[1] financial: %d rows x %d period cols | %d codes | indicators matched: %d",
  nrow(fin), length(per_cols), uniqueN(fin$code), length(IND))

fl <- melt(fin[指标 %in% unname(IND)], id.vars = c("code", "指标"),
           measure.vars = per_cols, variable.name = "period", value.name = "value",
           variable.factor = FALSE)
setnames(fl, "指标", "ind")
# 指标名 -> 英文 key
fl[, key := names(IND)[match(ind, IND)]]
fl <- fl[!is.na(key)]
fl[, value := suppressWarnings(as.numeric(value))]
fl <- fl[is.finite(value)]
fl[, report_date := as.Date(period, format = "%Y%m%d")]
fl[, mm := as.integer(substr(period, 5L, 6L))]
# ★ PIT 修正：按报告期类型用法定披露截止日
fl[, avail_date := report_date + fcase(
  mm == 12L, 120L,      # 年报  12-31 -> 次年 4-30（原 +90 天会前视 31 天）
  mm ==  6L,  62L,      # 半年报 6-30 -> 8-31
  mm ==  3L,  30L,      # 一季报 3-31 -> 4-30
  default  =  31L)]     # 三季报 9-30 -> 10-31
P("[1] long table: %d rows | avail_date lag: annual +120d, H1 +62d, Q1 +30d, Q3 +31d",
  nrow(fl))

fl <- unique(fl, by = c("code", "report_date", "key"))   # 防 dcast 遇重复键静默聚合
fw <- dcast(fl, code + report_date + avail_date ~ key, value.var = "value")
P("[1] wide table: %d rows x %d cols (%.0f s)", nrow(fw), ncol(fw), el())

## ================= 2. 滚动匹配：每个交易日取"最新已披露"的财务 =================
setkey(fw, code, avail_date)
grid <- CJ(code = unique(fw$code), date = panel_days)
grid[, avail_date := date]
setkey(grid, code, avail_date)
matched <- fw[grid, roll = TRUE]                      # roll=TRUE: 向后取最近一条
matched <- matched[!is.na(report_date)]
P("[2] matched: %d (date,code) rows | coverage vs panel = %.1f%%",
  nrow(matched), 100 * nrow(matched) / nrow(grid))

## ================= 3. 构建因子 =================
## 不复权价（amount/turnover）—— close 是前复权价，不能配 BPS/EPS
pp <- daily[turnover > 0 & date %in% panel_days,
            .(date, code, px = amount / turnover)]
pp <- pp[is.finite(px) & px > 0]
m <- merge(matched, pp, by = c("date", "code"))
P("[3] merged with prices: %d rows", nrow(m))

m[, `:=`(
  value_bp    = bps / px,                             # 账面市值比
  value_ep    = eps / px,                             # 盈利收益率
  quality_roe = roe,
  quality_gpm = gpm,
  quality_npr = npr,
  growth_np   = npg,
  growth_rev  = rg,
  leverage    = -lev                                  # 取负：低杠杆得分高
)]
FACS_EXT <- c("value_bp", "value_ep", "quality_roe", "quality_gpm",
              "quality_npr", "growth_np", "growth_rev", "leverage")

## ---- 价格类因子（复用 build_style_exposures 的口径，但重算以便统一处理）----
dt <- copy(daily)[, .(code, date, close, amount, turnover)]
setorder(dt, code, date)
dt[, ret := close / shift(close) - 1, by = code]
dt[, momentum := shift(close, 21L) / shift(close, 271L) - 1, by = code]
dt[, m1 := frollmean(ret, 60L, align = "right", na.rm = TRUE), by = code]
dt[, m2 := frollmean(ret^2, 60L, align = "right", na.rm = TRUE), by = code]
dt[, nok := frollsum(as.integer(!is.na(ret)), 60L, align = "right"), by = code]
dt[, volatility := fifelse(nok == 60L, sqrt(pmax(m2 - m1^2, 0) * nok / (nok - 1L) * 252), NA_real_),
   by = code]
dt[, liquidity := frollmean(turnover, 60L, align = "right", na.rm = TRUE), by = code]
dt[, size := fifelse(turnover > 0, log(pmax(amount / turnover, 1e-6)), NA_real_)]
pf <- dt[date %in% panel_days, .(date, code, momentum, volatility, size, liquidity)]
FACS_PRICE <- c("momentum", "volatility", "size", "liquidity")

## ---- 合并（保持【宽表】，不 melt —— melt 会变成 1.1 亿行 / 5 GB）----
allf <- merge(m[, c("date", "code", FACS_EXT), with = FALSE], pf,
              by = c("date", "code"), all = FALSE)
dup <- nrow(allf) - uniqueN(allf, by = c("date", "code"))
if (dup > 0L) P("[3] WARN: %d duplicate (date,code) rows removed", dup)
allf <- unique(allf, by = c("date", "code"))
FACS <- c(FACS_EXT, FACS_PRICE)
P("[3] factor table: %d rows x %d cols | %d factors", nrow(allf), ncol(allf), length(FACS))

# 逐列 winsorize(1%/99%) + 按日横截面 z-score（宽表上直接做，无需 melt）
for (f in FACS) {
  allf[, (f) := {
    v <- get(f)
    q <- quantile(v, c(0.01, 0.99), na.rm = TRUE)
    v <- pmin(pmax(v, q[1]), q[2])
    s <- sd(v, na.rm = TRUE)
    if (is.finite(s) && s > 0) (v - mean(v, na.rm = TRUE)) / s else NA_real_
  }, by = date]
}
setorder(allf, date, code)
arrow::write_parquet(allf, "_factor_panel.parquet")
P("[3] -> _factor_panel.parquet (%.1f MB)", file.size("_factor_panel.parquet") / 1e6)

## ================= 4. 因子画像 =================
hr("FACTOR COVERAGE BY YEAR (%)  -- share of stocks with a finite value")
allf[, yr := year(date)]
cov <- rbindlist(lapply(FACS, function(f)
  allf[, .(factor = f, cov = mean(is.finite(get(f)))), by = yr]))
print(dcast(cov, yr ~ factor, value.var = "cov")[
  , lapply(.SD, function(x) round(100 * x, 1)), by = yr])

hr("PAIRWISE CORRELATION OF FACTORS (pooled, complete rows only)")
M <- as.matrix(allf[, ..FACS])
keep <- rowSums(!is.finite(M)) == 0L
P("  rows used = %d / %d", sum(keep), nrow(M))
cm <- cor(M[keep, , drop = FALSE])
print(round(cm, 2))
P("\n  READ: |corr| near 1 between two factors => they carry the same information.")

hr("DONE")
P("Total elapsed %.0f s (%.1f min)", el(), el() / 60)
P("NEXT: use _factor_panel.csv to build a multi-factor score, then re-run the")
P("      active-weight L3 evaluation to see whether conclusions change.")
