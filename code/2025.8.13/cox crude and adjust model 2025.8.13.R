# ==== Cox 回归：Crude + Model1/2/3/4（协变量沿用 logistic） ====

# Packages
library(survival)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tidyr)
library(broom)

# ---------- 1) 读取数据 ----------
file_path <- "output/2025.8.11/participant_with_htn_med_flags.csv"
df <- read_csv(file_path)

# 基础预处理（与你此前一致）
df <- df %>%
  mutate(
    is_RPL = as.factor(is_RPL),
    participant.p20116_i0 = factor(participant.p20116_i0,
                                   levels = c(0, 1, 2),
                                   labels = c("Never", "Previous", "Current")),
    participant.p2814_i0 = factor(participant.p2814_i0,
                                  levels = c(0, 1, -1, -3),
                                  labels = c("No", "Yes", "Don't know", "Prefer not"))
  )

# ---------- 2) 结局列表 ----------
cvd_list <- c(
  "is_participant.p131382","is_participant.p131348","is_participant.p131350","is_participant.p131352",
  "is_participant.p131342","is_participant.p131344","is_participant.p131354","is_participant.p131296",
  "is_participant.p131298","is_participant.p131304","is_participant.p131306","is_participant.p131314",
  "is_participant.p131316","is_participant.p131386","is_participant.p131366","is_participant.p131362",
  "is_participant.p131056","is_participant.p131322","is_participant.p131324","is_participant.p131328",
  "is_participant.p131400","is_participant.p131308"
)

# ---------- 3) 协变量集合（沿用你logistic的Model1~4） ----------
model1_cov <- c(
  "participant.p21022",
  "participant.p21001_i0",
  "participant.p20116_i0",
  "participant.p1558_i0...3",
  "participant.p1488_i0",
  "participant.p1160_i0",
  "participant.p1498_i0",
  "participant.p884_i0"
)
model2_cov <- c(model1_cov,
                "participant.p4080_i0_a1","participant.p4079_i0_a1",
                "participant.p30690_i0","participant.p30760_i0","participant.p30710_i0"
)
model3_cov <- c(model2_cov, "hypertension_status")
model4_cov <- c(model3_cov,
                "participant.p3581_i0","participant.p2814_i0","participant.p26223",
                "participant.p22009_a1","participant.p22009_a2","participant.p22009_a3","participant.p22009_a4",
                "participant.p22009_a5","participant.p22009_a6","participant.p22009_a7","participant.p22009_a8",
                "participant.p22009_a9","participant.p22009_a10"
)

# ---------- 4) 构造公式 ----------
get_formula <- function(outcome, covars) {
  base <- paste0("Surv(follow_up_new, ", outcome, ") ~ is_RPL")
  if (length(covars) == 0) {
    as.formula(base)
  } else {
    as.formula(paste0(base, " + ", paste(covars, collapse = " + ")))
  }
}

# ---------- 5) 提取结果（稳健：兼容 is_RPL / is_RPL1；避免 confint 降维） ----------
extract_result <- function(fit, outcome, model_label) {
  if (is.null(fit)) {
    return(data.frame(CVD = outcome, Model = model_label,
                      HR_CI = NA_character_, P_value = NA_character_,
                      stringsAsFactors = FALSE))
  }
  s  <- summary(fit)
  rn <- rownames(s$coefficients)
  idx <- grep("^is_RPL", rn)
  if (length(idx) == 0) {
    return(data.frame(CVD = outcome, Model = model_label,
                      HR_CI = NA_character_, P_value = NA_character_,
                      stringsAsFactors = FALSE))
  }
  
  beta <- coef(fit)[idx][1]
  hr   <- round(exp(beta), 2)
  
  ci_mat <- suppressMessages(confint(fit))
  ci_sel <- ci_mat[idx, , drop = FALSE][1, , drop = FALSE]
  ci_low  <- round(exp(ci_sel[1, 1, drop = TRUE]), 2)
  ci_high <- round(exp(ci_sel[1, 2, drop = TRUE]), 2)
  
  p <- signif(s$coefficients[idx, "Pr(>|z|)"][1], 3)
  
  data.frame(
    CVD = outcome,
    Model = model_label,
    HR_CI = sprintf("%s (%.2f–%.2f)", hr, ci_low, ci_high),
    P_value = format.pval(p, digits = 3, eps = 0.001),
    stringsAsFactors = FALSE
  )
}

# ---------- 6) 单个结局的拟合流程 ----------
na_rows_for_all_models <- function(outcome) {
  dplyr::bind_rows(
    data.frame(CVD = outcome, Model = "Crude",   HR_CI = NA, P_value = NA),
    data.frame(CVD = outcome, Model = "Model 1", HR_CI = NA, P_value = NA),
    data.frame(CVD = outcome, Model = "Model 2", HR_CI = NA, P_value = NA),
    data.frame(CVD = outcome, Model = "Model 3", HR_CI = NA, P_value = NA),
    data.frame(CVD = outcome, Model = "Model 4", HR_CI = NA, P_value = NA)
  )
}

run_one_outcome <- function(outcome) {
  message("Fitting: ", outcome)
  
  needed_cols <- unique(c("follow_up_new", outcome, "is_RPL",
                          model1_cov, model2_cov, model3_cov, model4_cov))
  
  # 显式使用 dplyr:: 前缀，防止与 MASS::select 冲突
  df_use <- df %>%
    dplyr::select(dplyr::any_of(needed_cols)) %>%
    dplyr::filter(
      !is.na(.data[["follow_up_new"]]),
      !is.na(.data[[outcome]]),
      !is.na(.data[["is_RPL"]])
    )
  
  # 若结局没有变异（全0/全1）或样本极少，返回 NA 行
  if (nrow(df_use) <= 50 || length(unique(df_use[[outcome]])) < 2) {
    return(na_rows_for_all_models(outcome))
  }
  
  # 各模型可用协变量（只取当前数据存在的列）
  m1_cov <- intersect(model1_cov, names(df_use))
  m2_cov <- intersect(model2_cov, names(df_use))
  m3_cov <- intersect(model3_cov, names(df_use))
  m4_cov <- intersect(model4_cov, names(df_use))
  
  # 拟合各模型
  fit_crude <- tryCatch(coxph(get_formula(outcome, character(0)), data = df_use), error = function(e) NULL)
  fit_m1    <- tryCatch(coxph(get_formula(outcome, m1_cov),        data = df_use), error = function(e) NULL)
  fit_m2    <- tryCatch(coxph(get_formula(outcome, m2_cov),        data = df_use), error = function(e) NULL)
  fit_m3    <- tryCatch(coxph(get_formula(outcome, m3_cov),        data = df_use), error = function(e) NULL)
  fit_m4    <- tryCatch(coxph(get_formula(outcome, m4_cov),        data = df_use), error = function(e) NULL)
  
  # 收集结果
  dplyr::bind_rows(
    extract_result(fit_crude, outcome, "Crude"),
    extract_result(fit_m1,    outcome, "Model 1"),
    extract_result(fit_m2,    outcome, "Model 2"),
    extract_result(fit_m3,    outcome, "Model 3"),
    extract_result(fit_m4,    outcome, "Model 4")
  )
}

# ---------- 7) 批量运行并保存 ----------
results <- purrr::map_dfr(cvd_list, run_one_outcome)

dir.create("output/2025.8.13", showWarnings = FALSE, recursive = TRUE)
out_path <- "output/2025.8.13/CVD_RPL_Cox_Subgroups_results_Model0to4_logiCov.csv"
readr::write_csv(results, out_path)

message("Done. Saved to: ", out_path)
