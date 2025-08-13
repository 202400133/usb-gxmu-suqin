# --- libraries ---
library(dplyr)
library(readr)
library(broom)
library(MASS)     # stepAIC
# 如果需要Firth作为兜底（完全分离），取消下一行注释并确保已安装 logistf 包
# library(logistf)

# 读取数据
data <- read_csv("output/2025.8.11/participant_with_htn_med_flags.csv")

# 结局（22种CVD）
cvd_vars <- c(
  "is_participant.p131382","is_participant.p131348","is_participant.p131350","is_participant.p131352",
  "is_participant.p131342","is_participant.p131344","is_participant.p131354","is_participant.p131296",
  "is_participant.p131298","is_participant.p131304","is_participant.p131306","is_participant.p131314",
  "is_participant.p131316","is_participant.p131386","is_participant.p131366","is_participant.p131362",
  "is_participant.p131056","is_participant.p131322","is_participant.p131324","is_participant.p131328",
  "is_participant.p131400","is_participant.p131308"
)

# Model 1 候选协变量
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

# 结果容器
results <- tibble()

# 循环各个CVD
for (cvd in cvd_vars) {
  
  # 仅保留当前需要的列并删缺失
  keep_cols <- c(cvd, "is_RPL", model1_cov)
  df <- data %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    tidyr::drop_na()
  
  # 样本太少跳过
  if (nrow(df) <= 50) {
    results <- dplyr::bind_rows(
      results,
      tibble(
        CVD = cvd, Model = "Model 1 (forward AIC)",
        OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
        OR_95CI = NA_character_, P_value = NA_character_,
        Selected_Covariates = NA_character_, AIC = NA_real_
      )
    )
    next
  }
  
  # 起始模型：仅含 is_RPL
  f_null <- as.formula(paste(cvd, "~ is_RPL"))
  m0 <- try(glm(f_null, data = df, family = binomial), silent = TRUE)
  
  if (inherits(m0, "try-error")) {
    results <- dplyr::bind_rows(
      results,
      tibble(
        CVD = cvd, Model = "Model 1 (forward AIC)",
        OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
        OR_95CI = NA_character_, P_value = NA_character_,
        Selected_Covariates = NA_character_, AIC = NA_real_
      )
    )
    next
  }
  
  # 前向逐步：候选项仅来自 model1_cov（is_RPL 固定保留）
  upper_scope <- as.formula(paste("~", paste(c("is_RPL", intersect(model1_cov, names(df))), collapse = " + ")))
  m_step <- try(suppressWarnings(stepAIC(m0, scope = upper_scope, direction = "forward", trace = FALSE)),
                silent = TRUE)
  
  # 如果 stepAIC 失败，就退回 m0
  if (inherits(m_step, "try-error")) m_step <- m0
  
  # 提取 is_RPL 的OR/CI/P
  tidy_fit <- broom::tidy(m_step, conf.int = TRUE, exponentiate = TRUE)
  rpl_row  <- dplyr::filter(tidy_fit, term == "is_RPL")
  
  if (nrow(rpl_row) == 1) {
    OR    <- round(rpl_row$estimate, 2)
    lower <- round(rpl_row$conf.low, 2)
    upper <- round(rpl_row$conf.high, 2)
    pval  <- formatC(rpl_row$p.value, format = "e", digits = 2)
    or_ci <- paste0(OR, " (", lower, "-", upper, ")")
  } else {
    # 若 is_RPL 被异常删除或完全分离，尝试兜底（可选：Firth）
    OR <- lower <- upper <- NA_real_
    or_ci <- NA_character_
    pval  <- NA_character_
    # 可选兜底：Firth
    # fit_firth <- try(logistf(as.formula(paste(cvd, "~ is_RPL +", paste(setdiff(attr(terms(m_step), "term.labels"), "is_RPL"), collapse = " + "))), data = df), silent = TRUE)
    # ...
  }
  
  # 记录被选入的协变量（不含 is_RPL）
  sel_cov <- setdiff(attr(terms(m_step), "term.labels"), "is_RPL")
  sel_cov_str <- if (length(sel_cov) == 0) "(none)" else paste(sel_cov, collapse = " + ")
  
  results <- dplyr::bind_rows(
    results,
    tibble(
      CVD = cvd,
      Model = "Model 1 (forward AIC)",
      OR = OR, Lower_CI = lower, Upper_CI = upper,
      OR_95CI = or_ci, P_value = pval,
      Selected_Covariates = sel_cov_str,
      AIC = AIC(m_step)
    )
  )
}

# 保存
readr::write_csv(results, "output/2025.8.12/CVD_RPL_Model1_forward_step.csv")
