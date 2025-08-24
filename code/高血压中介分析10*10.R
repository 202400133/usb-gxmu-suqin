#结果放一起
OUTDIR <- "output/2025.8.21/results"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================
# RPL → 高血压 → CVD 中介分析（全协变量 + mice 插补 + 10x10）
# ============================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival)
  library(mice);  library(regmedint); library(mitools)
})

set.seed(20250823)

# ---------- 1) 数据读取与处理 ----------
infile <- "output/2025.8.21/mediator_htn_filtered.csv"
df <- read_csv(infile, show_col_types = FALSE)

df <- df %>%
  mutate(
    time  = as.numeric(follow_up_new),
    event = as.integer(ifelse(is.na(is_CVD_new), 0L, is_CVD_new)),
    is_RPL    = as.integer(is_RPL),
    htn_valid = as.integer(htn_valid)
  ) %>%
  filter(!is.na(time) & time >= 0)

# ---------- 2) 提取并重命名协变量 ----------
df <- df %>%
  mutate(
    age     = as.numeric(participant.p21022),
    ethnicity = as.factor(is_white),
    SES     = as.numeric(participant.p22189),
    smoking = as.numeric(participant.p20116_i0),
    alcohol = as.numeric(participant.p1558_i0),
    physical_activity = as.numeric(participant.p884_i0),
    sleep   = as.numeric(participant.p1160_i0),
    tea     = as.numeric(participant.p1488_i0),
    coffee  = as.numeric(participant.p1498_i0),
    BMI     = as.numeric(participant.p21001_i0),
    PRS_CVD = as.numeric(participant.p26223),
    PC1 = participant.p22009_a1,
    PC2 = participant.p22009_a2,
    PC3 = participant.p22009_a3,
    PC4 = participant.p22009_a4,
    PC5 = participant.p22009_a5,
    PC6 = participant.p22009_a6,
    PC7 = participant.p22009_a7,
    PC8 = participant.p22009_a8,
    PC9 = participant.p22009_a9,
    PC10 = participant.p22009_a10,
    centre = as.factor(centre),
    geno_batch = as.factor(geno_batch),
    baseline_year = as.numeric(baseline_year_f),
    baseline_qtr  = as.factor(baseline_qtr_f)
  )

# ---------- 3) 构建分析用数据集 ----------
dat <- df %>% select(
  time, event, is_RPL, htn_valid,
  age, ethnicity, SES, smoking, alcohol, physical_activity, sleep, tea, coffee, BMI, PRS_CVD,
  PC1:PC10, centre, geno_batch, baseline_year, baseline_qtr
)

# ---------- 4) mice 多重插补（m=10, maxit=10） ----------
meth <- make.method(dat)
pred <- make.predictorMatrix(dat)
meth[c("time","event","is_RPL","htn_valid")] <- ""   # 不插补 Y/A/M

imp <- mice(dat, m = 10, maxit = 10, method = meth, predictorMatrix = pred, printFlag = TRUE)

# ---------- 5) 每个插补集内：协变量标准化 ----------
mi_list <- complete(imp, "all")
mi_list <- lapply(mi_list, function(dd){
  # 标准化 numeric 类型变量
  num_vars <- c("age","SES","smoking","alcohol","physical_activity","sleep",
                "tea","coffee","BMI","PRS_CVD","PC1","PC2","PC3","PC4","PC5",
                "PC6","PC7","PC8","PC9","PC10","baseline_year")
  for (v in num_vars) dd[[paste0(v, "_z")]] <- as.numeric(scale(dd[[v]]))
  dd
})

# ---------- 6) 设置协变量及条件 ----------
cvars <- unlist(lapply(names(mi_list[[1]]), function(x) if (grepl("_z$", x)) x else NULL))
c_cond_vec <- rep(0, length(cvars))

# ---------- 7) 中介分析函数 ----------
fit_one <- function(data){
  data <- droplevels(data)
  do.call(regmedint, list(
    data = data,
    yvar = "time", eventvar = "event",
    avar = "is_RPL", mvar = "htn_valid",
    cvar = cvars,
    a0 = 0, a1 = 1, m_cde = 0,
    yreg = "survCox", mreg = "logistic",
    interaction = FALSE, casecontrol = FALSE,
    c_cond = c_cond_vec
  ))
}

safe_fit <- function(dd) tryCatch(fit_one(dd), error = function(e){ attr(e,"msg") <- conditionMessage(e); NULL })
fits <- lapply(mi_list, safe_fit)

# ---------- 8) 合并结果（不依赖 summary 的行名） ----------
wanted <- c("cde","pnde","tnie","tnde","pnie","te","pm")

coef_list <- lapply(fits[ok_idx], function(f){
  x <- coef(f)[wanted]
  setNames(as.numeric(x), wanted)
})

vcov_list <- lapply(fits[ok_idx], function(f){
  V <- vcov(f)[wanted, wanted, drop = FALSE]
  dimnames(V) <- list(wanted, wanted)
  V
})

mi_comb <- MIcombine(results = coef_list, variances = vcov_list)

# 直接从 MIcombine 抽取（顺序就是我们传入的 wanted）
est <- as.numeric(mi_comb$coefficients)
se  <- sqrt(diag(mi_comb$variance))

# 贴上名字，避免行名不匹配的问题
names(est) <- wanted
names(se)  <- wanted

# ---------- 9) 结果整理 ----------
res <- data.frame(
  effect   = wanted,
  estimate = est,
  se       = se,
  lower    = est - qnorm(0.975)*se,
  upper    = est + qnorm(0.975)*se,
  p.value  = 2*pnorm(abs(est/se), lower.tail = FALSE),
  stringsAsFactors = FALSE
)

# 除 pm 外转 HR
res$HR  <- ifelse(res$effect == "pm", NA_real_, exp(res$estimate))
res$LCL <- ifelse(res$effect == "pm", NA_real_, exp(res$lower))
res$UCL <- ifelse(res$effect == "pm", NA_real_, exp(res$upper))

# 仅输出 TE / NDE(pnde) / NIE(tnie) / PM
out4 <- res[match(c("te","pnde","tnie","pm"), res$effect), ]
print(out4[, c("effect","estimate","HR","LCL","UCL","p.value")], row.names = FALSE)

cat("\n—— 论文风格 ——\n")
fmt_hr <- function(est, l, u) sprintf("HR %.2f (%.2f–%.2f)", est, l, u)
cat("TE  :", fmt_hr(out4$HR[out4$effect=="te"],   out4$LCL[out4$effect=="te"],   out4$UCL[out4$effect=="te"]),   "\n")
cat("NDE :", fmt_hr(out4$HR[out4$effect=="pnde"], out4$LCL[out4$effect=="pnde"], out4$UCL[out4$effect=="pnde"]), "\n")
cat("NIE :", fmt_hr(out4$HR[out4$effect=="tnie"], out4$LCL[out4$effect=="tnie"], out4$UCL[out4$effect=="tnie"]), "\n")
cat("PM  :", sprintf("%.1f%% (%.1f%%–%.1f%%)",
                     100*res$estimate[res$effect=="pm"],
                     100*res$lower[res$effect=="pm"],
                     100*res$upper[res$effect=="pm"]), "\n")
# 全量指标（cde/pnde/tnie/tnde/pnie/te/pm，含 estimate、se、HR、95%CI、P值）
readr::write_csv(res,  file.path(OUTDIR, "mediation_HTN_full.csv"))

# 精简四项（TE / NDE / NIE / PM），便于论文三线表
readr::write_csv(out4, file.path(OUTDIR, "mediation_HTN_out4.csv"))

cat("\n已保存：",
    file.path(OUTDIR, "mediation_HTN_full.csv"), "和",
    file.path(OUTDIR, "mediation_HTN_out4.csv"), "\n")
