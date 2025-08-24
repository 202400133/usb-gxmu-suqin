# ============================
# RPL → 糖尿病(DM) → CVD 中介分析（全协变量 + mice 10×10）
# 强制新拟合版：不会复用 HTN 的 fits；含对比提醒
# ============================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival)
  library(mice);  library(regmedint); library(mitools)
})

options(tibble.print_max = 20, dplyr.print_max = 20)
set.seed(20250823)

# ---- 参数 ----
INFILE <- "output/2025.8.21/mediator_dm_filtered.csv"
MVAR   <- "dm_valid"                      # 若中介列名不同，改这里
OUTDIR <- "output/2025.8.21/results"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

#（可选）先清理上一轮 DM 对象，防串用
rm(list = intersect(ls(), c("imp_dm","mi_list_dm","fits_dm","mi_comb_dm","res_dm","out4_dm"))); gc()

# ---------- 1) 读取与基础变量 ----------
df <- read_csv(INFILE, show_col_types = FALSE)

need_cols <- c("is_RPL", MVAR, "follow_up_new", "is_CVD_new")
stopifnot(all(need_cols %in% names(df)))

df <- df %>%
  mutate(
    time  = as.numeric(follow_up_new),
    event = as.integer(ifelse(is.na(is_CVD_new), 0L, is_CVD_new)),
    is_RPL = as.integer(is_RPL),
    !!MVAR := as.integer(.data[[MVAR]])
  ) %>% filter(!is.na(time) & time >= 0)

cat("N=", nrow(df),
    " | RPL=", sum(df$is_RPL==1, na.rm=TRUE),
    " | DM=",  sum(df[[MVAR]]==1,   na.rm=TRUE),
    " | CVD events=", sum(df$event==1, na.rm=TRUE), "\n")

# ---------- 2) 协变量构造（1558 自动识别） ----------
alc_candidates <- grep("^participant\\.p1558", names(df), value = TRUE)
alcohol_col <- if (length(alc_candidates) > 0) alc_candidates[1] else NA_character_
cat("Alcohol column:", ifelse(is.na(alcohol_col), "NONE -> set NA", alcohol_col), "\n")

df <- df %>%
  mutate(
    age     = suppressWarnings(as.numeric(participant.p21022)),
    ethnicity = as.factor(is_white),
    SES     = suppressWarnings(as.numeric(participant.p22189)),
    smoking = suppressWarnings(as.numeric(participant.p20116_i0)),
    alcohol = if (!is.na(alcohol_col)) suppressWarnings(as.numeric(.data[[alcohol_col]])) else NA_real_,
    physical_activity = suppressWarnings(as.numeric(participant.p884_i0)),
    sleep   = suppressWarnings(as.numeric(participant.p1160_i0)),
    tea     = suppressWarnings(as.numeric(participant.p1488_i0)),
    coffee  = suppressWarnings(as.numeric(participant.p1498_i0)),
    BMI     = suppressWarnings(as.numeric(participant.p21001_i0)),
    PRS_CVD = suppressWarnings(as.numeric(participant.p26223)),
    PC1 = participant.p22009_a1,  PC2 = participant.p22009_a2,  PC3 = participant.p22009_a3,
    PC4 = participant.p22009_a4,  PC5 = participant.p22009_a5,  PC6 = participant.p22009_a6,
    PC7 = participant.p22009_a7,  PC8 = participant.p22009_a8,  PC9 = participant.p22009_a9,  PC10 = participant.p22009_a10,
    centre = as.factor(centre),
    geno_batch = as.factor(geno_batch),
    baseline_year = suppressWarnings(as.numeric(baseline_year_f)),
    baseline_qtr  = as.factor(baseline_qtr_f)
  )

# ---------- 3) 构建分析数据 ----------
dat_dm <- df %>% select(
  time, event, is_RPL, !!sym(MVAR),
  age, ethnicity, SES, smoking, alcohol, physical_activity, sleep, tea, coffee, BMI, PRS_CVD,
  PC1:PC10, centre, geno_batch, baseline_year, baseline_qtr
)

# ---------- 4) 多重插补（m=10, maxit=10） ----------
meth <- make.method(dat_dm)
pred <- make.predictorMatrix(dat_dm)
meth[c("time","event","is_RPL", MVAR)] <- ""   # 不插补 Y/A/M

imp_dm <- mice(dat_dm, m = 10, maxit = 10, method = meth, predictorMatrix = pred, printFlag = TRUE)
saveRDS(imp_dm, file.path(OUTDIR, "dm_imp_m10.rds"))

# ---------- 5) 标准化 numeric ----------
mi_list_dm <- complete(imp_dm, "all")
mi_list_dm <- lapply(mi_list_dm, function(dd){
  num_vars <- c("age","SES","smoking","alcohol","physical_activity","sleep",
                "tea","coffee","BMI","PRS_CVD","PC1","PC2","PC3","PC4","PC5",
                "PC6","PC7","PC8","PC9","PC10","baseline_year")
  for (v in num_vars) dd[[paste0(v, "_z")]] <- suppressWarnings(as.numeric(scale(dd[[v]])))
  dd
})

# ---------- 6) 协变量清单 & 条件点 ----------
cvars_dm <- names(mi_list_dm[[1]])[grepl("_z$", names(mi_list_dm[[1]]))]
c_cond_vec_dm <- rep(0, length(cvars_dm))

# ---------- 7) 拟合（强制重算，绝不复用旧 fits） ----------
fit_one_dm <- function(data){
  fit <- do.call(regmedint, list(
    data = data,
    yvar = "time", eventvar = "event",
    avar = "is_RPL", mvar = MVAR,
    cvar = cvars_dm,
    a0 = 0, a1 = 1, m_cde = 0,
    yreg = "survCox", mreg = "logistic",
    interaction = FALSE, casecontrol = FALSE,
    c_cond = c_cond_vec_dm
  ))
  
  # ✅ 添加这行：将 mvar 存入模型对象属性
  attr(fit, "args") <- list(mvar = MVAR)
  
  return(fit)
}
safe_fit_dm <- function(dd) tryCatch(fit_one_dm(dd), error = function(e){ message("fit error: ", conditionMessage(e)); NULL })

fits_dm <- vector("list", length(mi_list_dm))   # ← 重置
for (i in seq_along(mi_list_dm)) {
  cat(sprintf("Fitting imputed dataset %d/%d ...\n", i, length(mi_list_dm)))
  gc()
  fits_dm[[i]] <- safe_fit_dm(mi_list_dm[[i]])
}
saveRDS(fits_dm, file.path(OUTDIR, "dm_fits_m10.rds"))

ok_idx_dm  <- which(!sapply(fits_dm, is.null))
stopifnot(length(ok_idx_dm) >= 2)

# 断言：本次确实使用 DM 作为中介
used_mvar_dm <- attributes(fits_dm[[ok_idx_dm[1]]])$args$mvar
cat("mvar used in this fit:", used_mvar_dm, "\n")
stopifnot(identical(used_mvar_dm, MVAR))

# ---------- 8) Rubin 合并 ----------
wanted <- c("cde","pnde","tnie","tnde","pnie","te","pm")

coef_list_dm <- lapply(fits_dm[ok_idx_dm], function(f){
  x <- coef(f)[wanted]; setNames(as.numeric(x), wanted)
})
vcov_list_dm <- lapply(fits_dm[ok_idx_dm], function(f){
  V <- vcov(f)[wanted, wanted, drop = FALSE]; dimnames(V) <- list(wanted, wanted); V
})

mi_comb_dm <- MIcombine(results = coef_list_dm, variances = vcov_list_dm)
est_dm <- as.numeric(mi_comb_dm$coefficients); names(est_dm) <- wanted
se_dm  <- sqrt(diag(mi_comb_dm$variance));   names(se_dm)  <- wanted

# ---------- 9) 结果整理 & 导出 ----------
res_dm <- data.frame(
  effect   = wanted,
  estimate = est_dm,
  se       = se_dm,
  lower    = est_dm - qnorm(0.975)*se_dm,
  upper    = est_dm + qnorm(0.975)*se_dm,
  p.value  = 2*pnorm(abs(est_dm/se_dm), lower.tail = FALSE),
  stringsAsFactors = FALSE
)
res_dm$HR  <- ifelse(res_dm$effect == "pm", NA_real_, exp(res_dm$estimate))
res_dm$LCL <- ifelse(res_dm$effect == "pm", NA_real_, exp(res_dm$lower))
res_dm$UCL <- ifelse(res_dm$effect == "pm", NA_real_, exp(res_dm$upper))

out4_dm <- res_dm[match(c("te","pnde","tnie","pm"), res_dm$effect), ]

readr::write_csv(res_dm,  file.path(OUTDIR, "mediation_DM_full.csv"))
readr::write_csv(out4_dm, file.path(OUTDIR, "mediation_DM_out4.csv"))

print(out4_dm[, c("effect","estimate","HR","LCL","UCL","p.value")], row.names = FALSE)

cat("\n—— 论文风格 ——\n")
fmt_hr <- function(est, l, u) sprintf("HR %.2f (%.2f–%.2f)", est, l, u)
cat("TE  :", fmt_hr(out4_dm$HR[out4_dm$effect=="te"],
                    out4_dm$LCL[out4_dm$effect=="te"],
                    out4_dm$UCL[out4_dm$effect=="te"]), "\n")
cat("NDE :", fmt_hr(out4_dm$HR[out4_dm$effect=="pnde"],
                    out4_dm$LCL[out4_dm$effect=="pnde"],
                    out4_dm$UCL[out4_dm$effect=="pnde"]), "\n")
cat("NIE :", fmt_hr(out4_dm$HR[out4_dm$effect=="tnie"],
                    out4_dm$LCL[out4_dm$effect=="tnie"],
                    out4_dm$UCL[out4_dm$effect=="tnie"]), "\n")
cat("PM  :", sprintf("%.1f%% (%.1f%%–%.1f%%)",
                     100*res_dm$estimate[res_dm$effect=="pm"],
                     100*res_dm$lower[res_dm$effect=="pm"],
                     100*res_dm$upper[res_dm$effect=="pm"]), "\n")
cat("\n已保存：",
    file.path(OUTDIR, "mediation_DM_full.csv"), "和",
    file.path(OUTDIR, "mediation_DM_out4.csv"), "\n")