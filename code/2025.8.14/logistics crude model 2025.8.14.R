# --- libraries ---
library(dplyr)
library(readr)
library(broom)
library(tidyr)
library(stringr)

# ===== 路径与输入 =====
input_file  <- "output/2025.8.14/merged_filtered_with_Date5.csv"
output_dir  <- "output/2025.8.14"
output_file <- file.path(output_dir, "CVD_RPL_Crude_with_counts.csv")

if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 读取数据
data <- read_csv(input_file, show_col_types = FALSE)

# 结局（22种 CVD，二分类0/1）
cvd_vars <- c(
  "is_participant.p131382","is_participant.p131348","is_participant.p131350","is_participant.p131352",
  "is_participant.p131342","is_participant.p131344","is_participant.p131354","is_participant.p131296",
  "is_participant.p131298","is_participant.p131304","is_participant.p131306","is_participant.p131314",
  "is_participant.p131316","is_participant.p131386","is_participant.p131366","is_participant.p131362",
  "is_participant.p131056","is_participant.p131322","is_participant.p131324","is_participant.p131328",
  "is_participant.p131400","is_participant.p131308"
)

# 提取 is_RPL 的行（无论是否被当作因子产生“is_RPL1”等名）
get_rpl_row <- function(tidy_df) {
  r <- tidy_df %>% dplyr::filter(term == "is_RPL" | grepl("^is_RPL", term))
  if (nrow(r) == 0) return(NULL) else return(r[1,])
}

fmt_prop <- function(x) formatC(x, format = "f", digits = 3)

res_crude <- tibble()

for (cvd in cvd_vars) {
  keep_cols <- c(cvd, "is_RPL")
  df <- data %>%
    dplyr::select(dplyr::any_of(keep_cols)) %>%
    tidyr::drop_na()
  
  n_total <- nrow(df)
  
  # 预计算分组病例数与基数
  # 假设 is_RPL：0=Non-RPL, 1=RPL
  grp <- df %>%
    group_by(is_RPL) %>%
    summarise(
      n_group = n(),
      cases   = sum(.data[[cvd]] == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 提取两组（若你的 is_RPL 编码不同，这两行改一下）
  rpl_n      <- grp %>% filter(is_RPL == 1) %>% pull(n_group)
  rpl_cases  <- grp %>% filter(is_RPL == 1) %>% pull(cases)
  nonrpl_n   <- grp %>% filter(is_RPL == 0) %>% pull(n_group)
  nonrpl_cases <- grp %>% filter(is_RPL == 0) %>% pull(cases)
  
  # 若存在缺组，置0
  if (length(rpl_n) == 0)        rpl_n <- 0
  if (length(rpl_cases) == 0)    rpl_cases <- 0
  if (length(nonrpl_n) == 0)     nonrpl_n <- 0
  if (length(nonrpl_cases) == 0) nonrpl_cases <- 0
  
  # 比例
  rpl_prop    <- ifelse(rpl_n > 0, rpl_cases / rpl_n, NA_real_)
  nonrpl_prop <- ifelse(nonrpl_n > 0, nonrpl_cases / nonrpl_n, NA_real_)
  
  # 极端小样本：仅返回计数信息
  if (n_total == 0 || n_total <= 50) {
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = ifelse(n_total == 0, NA_integer_, rpl_cases + nonrpl_cases),
      NonEvents = ifelse(n_total == 0, NA_integer_, n_total - (rpl_cases + nonrpl_cases)),
      RPL_cases = rpl_cases,
      RPL_n = rpl_n,
      RPL_cases_prop = ifelse(is.na(rpl_prop), NA_character_, paste0(rpl_cases, " (", fmt_prop(rpl_prop), ")")),
      NonRPL_cases = nonrpl_cases,
      NonRPL_n = nonrpl_n,
      NonRPL_cases_prop = ifelse(is.na(nonrpl_prop), NA_character_, paste0(nonrpl_cases, " (", fmt_prop(nonrpl_prop), ")")),
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
      Events = rpl_cases + nonrpl_cases,
      NonEvents = n_total - (rpl_cases + nonrpl_cases),
      RPL_cases = rpl_cases,
      RPL_n = rpl_n,
      RPL_cases_prop = ifelse(is.na(rpl_prop), NA_character_, paste0(rpl_cases, " (", fmt_prop(rpl_prop), ")")),
      NonRPL_cases = nonrpl_cases,
      NonRPL_n = nonrpl_n,
      NonRPL_cases_prop = ifelse(is.na(nonrpl_prop), NA_character_, paste0(nonrpl_cases, " (", fmt_prop(nonrpl_prop), ")")),
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
    next
  }
  
  tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
  rpl_row  <- get_rpl_row(tidy_fit)
  
  if (is.null(rpl_row)) {
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = rpl_cases + nonrpl_cases,
      NonEvents = n_total - (rpl_cases + nonrpl_cases),
      RPL_cases = rpl_cases,
      RPL_n = rpl_n,
      RPL_cases_prop = ifelse(is.na(rpl_prop), NA_character_, paste0(rpl_cases, " (", fmt_prop(rpl_prop), ")")),
      NonRPL_cases = nonrpl_cases,
      NonRPL_n = nonrpl_n,
      NonRPL_cases_prop = ifelse(is.na(nonrpl_prop), NA_character_, paste0(nonrpl_cases, " (", fmt_prop(nonrpl_prop), ")")),
      OR = NA_real_, Lower_CI = NA_real_, Upper_CI = NA_real_,
      OR_95CI = NA_character_,
      P_value = NA_character_, P_value_decimal = NA_character_
    ))
  } else {
    OR    <- round(rpl_row$estimate, 2)
    lower <- round(rpl_row$conf.low,  2)
    upper <- round(rpl_row$conf.high, 2)
    p_e   <- formatC(rpl_row$p.value, format = "e", digits = 2)  # 科学计数
    p_dec <- formatC(rpl_row$p.value, format = "f", digits = 6)  # 十进制
    
    res_crude <- bind_rows(res_crude, tibble(
      CVD = cvd, Model = "Crude",
      N = n_total,
      Events = rpl_cases + nonrpl_cases,
      NonEvents = n_total - (rpl_cases + nonrpl_cases),
      RPL_cases = rpl_cases,
      RPL_n = rpl_n,
      RPL_cases_prop = paste0(rpl_cases, " (", fmt_prop(rpl_prop), ")"),
      NonRPL_cases = nonrpl_cases,
      NonRPL_n = nonrpl_n,
      NonRPL_cases_prop = paste0(nonrpl_cases, " (", fmt_prop(nonrpl_prop), ")"),
      OR = OR, Lower_CI = lower, Upper_CI = upper,
      OR_95CI = paste0(OR, " (", lower, "-", upper, ")"),
      P_value = p_e, P_value_decimal = p_dec
    ))
  }
}
# 保存结果
write_csv(res_crude, output_file)

cat("已完成：", output_file, "\n")
