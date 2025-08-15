# --- libraries ---
library(dplyr)
library(readr)
library(broom)
library(MASS)   # stepAIC
# 可选兜底（完全分离）：library(logistf)

# 读取数据
data <- read_csv("output/2025.8.14/merged_filtered_with_Date5.csv")

# 结局（22种 CVD）
cvd_vars <- c(
  "is_participant.p131382","is_participant.p131348","is_participant.p131350","is_participant.p131352",
  "is_participant.p131342","is_participant.p131344","is_participant.p131354","is_participant.p131296",
  "is_participant.p131298","is_participant.p131304","is_participant.p131306","is_participant.p131314",
  "is_participant.p131316","is_participant.p131386","is_participant.p131366","is_participant.p131362",
  "is_participant.p131056","is_participant.p131322","is_participant.p131324","is_participant.p131328",
  "is_participant.p131400","is_participant.p131308"
)

# Model 1 协变量
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

# Model 2 = Model1 + 4080/4079/30690/30760/30710
model2_cov <- c(
  model1_cov,
  "participant.p4080_i0_a1","participant.p4079_i0_a1",
  "participant.p30690_i0","participant.p30760_i0","participant.p30710_i0"
)

# Model 3 = Model2 + 高血压状态
model3_cov <- c(model2_cov, "hypertension_status")

# Model 4 候选协变量 = Model3 + 其余变量
model4_cov <- c(
  model3_cov,
  "participant.p3581_i0","participant.p2814_i0","participant.p26223",
  "participant.p22009_a1","participant.p22009_a2","participant.p22009_a3","participant.p22009_a4",
  "participant.p22009_a5","participant.p22009_a6","participant.p22009_a7","participant.p22009_a8",
  "participant.p22009_a9","participant.p22009_a10"
)

# 结果容器
results_m4 <- tibble()

for (cvd in cvd_vars) {
  
  # 仅保留需要列并去缺失；兼容缺失列名
  keep_cols <- c(cvd, "is_RPL", model4_cov)
  df <- data %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    tidyr::drop_na()
  
  if (nrow(df) <= 50) {
    results_m4 <- bind_rows(
      results_m4,
      tibble(
        CVD = cvd, Model = "Model 4 (forward AIC)",
        OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
        OR_95CI = NA_character_, P_value = NA_character_,
        P_value_decimal = NA_character_,
        Selected_Covariates = NA_character_, AIC = NA_real_
      )
    )
    next
  }
  
  # 起始模型：仅 is_RPL
  f_null <- as.formula(paste(cvd, "~ is_RPL"))
  m0 <- try(glm(f_null, data = df, family = binomial), silent = TRUE)
  if (inherits(m0, "try-error")) {
    results_m4 <- bind_rows(
      results_m4,
      tibble(
        CVD = cvd, Model = "Model 4 (forward AIC)",
        OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
        OR_95CI = NA_character_, P_value = NA_character_,
        P_value_decimal = NA_character_,
        Selected_Covariates = NA_character_, AIC = NA_real_
      )
    )
    next
  }
  
  # 前向逐步的上界：is_RPL + 所有可用的 Model4 候选变量
  usable_cov <- intersect(model4_cov, names(df))
  upper_scope <- as.formula(paste("~", paste(c("is_RPL", usable_cov), collapse = " + ")))
  
  m_step <- try(suppressWarnings(stepAIC(m0, scope = upper_scope, direction = "forward", trace = FALSE)),
                silent = TRUE)
  if (inherits(m_step, "try-error")) m_step <- m0
  
  # 提取 is_RPL 的 OR/CI/P
  tidy_fit <- broom::tidy(m_step, conf.int = TRUE, exponentiate = TRUE)
  rpl_row  <- dplyr::filter(tidy_fit, term == "is_RPL")
  
  if (nrow(rpl_row) == 1) {
    OR    <- round(rpl_row$estimate, 2)
    lower <- round(rpl_row$conf.low,  2)
    upper <- round(rpl_row$conf.high, 2)
    p_e   <- formatC(rpl_row$p.value, format = "e", digits = 2)   # 科学计数
    p_dec <- formatC(rpl_row$p.value, format = "f", digits = 6)   # 十进制
    or_ci <- paste0(OR, " (", lower, "-", upper, ")")
  } else {
    OR <- lower <- upper <- NA_real_
    p_e <- p_dec <- NA_character_
    or_ci <- NA_character_
    # 如遇完全分离，可考虑 Firth：logistf(...)
  }
  
  # 记录被选入协变量（不含 is_RPL）
  sel_cov <- setdiff(attr(terms(m_step), "term.labels"), "is_RPL")
  sel_cov_str <- if (length(sel_cov) == 0) "(none)" else paste(sel_cov, collapse = " + ")
  
  results_m4 <- bind_rows(
    results_m4,
    tibble(
      CVD = cvd,
      Model = "Model 4 (forward AIC)",
      OR = OR, Lower_CI = lower, Upper_CI = upper,
      OR_95CI = or_ci,
      P_value = p_e,
      P_value_decimal = p_dec,
      Selected_Covariates = sel_cov_str,
      AIC = AIC(m_step)
    )
  )
}

# 保存结果
write_csv(results_m4, "output/2025.8.14/CVD_RPL_Model4_forward_step2025.8.14.csv")