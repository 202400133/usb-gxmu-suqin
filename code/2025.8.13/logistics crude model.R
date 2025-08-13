# --- libraries ---
library(dplyr)
library(readr)
library(broom)

# 读取数据
data <- read_csv("output/2025.8.11/participant_with_htn_med_flags.csv")

# 结局（22种 CVD，二分类0/1）
cvd_vars <- c(
  "is_participant.p131382","is_participant.p131348","is_participant.p131350","is_participant.p131352",
  "is_participant.p131342","is_participant.p131344","is_participant.p131354","is_participant.p131296",
  "is_participant.p131298","is_participant.p131304","is_participant.p131306","is_participant.p131314",
  "is_participant.p131316","is_participant.p131386","is_participant.p131366","is_participant.p131362",
  "is_participant.p131056","is_participant.p131322","is_participant.p131324","is_participant.p131328",
  "is_participant.p131400","is_participant.p131308"
)

# 一个安全的小函数：提取 is_RPL 的行（无论因是否被当作因子产生“is_RPL1”等名）
get_rpl_row <- function(tidy_df) {
  r <- tidy_df %>% dplyr::filter(term == "is_RPL" | grepl("^is_RPL", term))
  if (nrow(r) == 0) return(NULL) else return(r[1,])
}

# 结果容器
res_crude <- tibble()

for (cvd in cvd_vars) {
  # 仅保留所需列，去缺失
  keep_cols <- c(cvd, "is_RPL")
  df <- data %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    tidyr::drop_na()
  
  # 样本量与事件数
  n_total  <- nrow(df)
  if (n_total == 0) {
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = 0, Events = NA_integer_, NonEvents = NA_integer_,
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
    next
  }
  
  # 确保结局是0/1，is_RPL保留原始编码（建议是0/1；若是字符/因子也可）
  # （若你的 is_RPL 是数值 0/1，此处无需改；保持即可）
  # df[[cvd]] <- as.integer(df[[cvd]] > 0)  # 如需强制二值化可解注
  
  # 极端小样本跳过
  if (n_total <= 50) {
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = sum(df[[cvd]] == 1, na.rm = TRUE),
      NonEvents = sum(df[[cvd]] == 0, na.rm = TRUE),
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
    next
  }
  
  # 粗模型：结局 ~ is_RPL
  form <- as.formula(paste(cvd, "~ is_RPL"))
  
  fit <- try(glm(formula = form, data = df, family = binomial), silent = TRUE)
  if (inherits(fit, "try-error")) {
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = sum(df[[cvd]] == 1, na.rm = TRUE),
      NonEvents = sum(df[[cvd]] == 0, na.rm = TRUE),
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
    next
  }
  
  tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  rpl_row  <- get_rpl_row(tidy_fit)
  
  if (is.null(rpl_row)) {
    #（少见）若提取不到 is_RPL（比如完全分离），返回 NA
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = sum(df[[cvd]] == 1, na.rm = TRUE),
      NonEvents = sum(df[[cvd]] == 0, na.rm = TRUE),
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
  } else {
    OR    <- round(rpl_row$estimate, 2)
    lower <- round(rpl_row$conf.low,  2)
    upper <- round(rpl_row$conf.high, 2)
    p_e   <- formatC(rpl_row$p.value, format = "e", digits = 2)  # 科学计数法
    p_dec <- formatC(rpl_row$p.value, format = "f", digits = 6)  # 十进制
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = sum(df[[cvd]] == 1, na.rm = TRUE),
      NonEvents = sum(df[[cvd]] == 0, na.rm = TRUE),
      OR = OR, Lower_CI = lower, Upper_CI = upper,
      OR_95CI = paste0(OR, " (", lower, "-", upper, ")"),
      P_value = p_e, P_value_decimal = p_dec
    ))
  }
}

# 保存结果
readr::write_csv(res_crude, "output/2025.8.12/CVD_RPL_Crude_only.csv")
