# ==== 准备包 ====
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)

# ==== 读取 ====
df <- read_csv("output/2025.8.11/participant_with_htn_med_flags.csv",
               show_col_types = FALSE) %>%
  mutate(is_RPL = factor(is_RPL, levels = c(0,1), labels = c("非RPL组","RPL组")))

# ==== 工具函数 ====
fmt_pct <- function(x, denom, digits = 1){
  pct <- ifelse(denom > 0, 100 * x / denom, NA_real_)
  sprintf("%d (%.1f%%)", x, pct)
}
fmt_p <- function(p){
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

# ==== 1) 饮酒频次（χ²） ====
alcohol_var <- "participant.p1558_i0...3"
alc_map <- function(v){
  case_when(
    v %in% c(1,2) ~ "≥3 times/week",
    v %in% c(3,4,5) ~ "≤2 times/week",
    v == 6 ~ "Never",
    v == -3 ~ NA_character_,
    TRUE ~ NA_character_
  )
}

alc_df <- df %>%
  transmute(is_RPL,
            alc_cat = alc_map(.data[[alcohol_var]])) %>%
  filter(!is.na(alc_cat))

# 分组内（各自非缺失为分母）计数
denom_by_grp <- alc_df %>% count(is_RPL, name = "denom")

alc_cnt <- alc_df %>%
  count(is_RPL, alc_cat, name = "n") %>%
  left_join(denom_by_grp, by = "is_RPL") %>%
  mutate(cell = fmt_pct(n, denom)) %>%
  select(is_RPL, alc_cat, cell) %>%
  pivot_wider(names_from = is_RPL, values_from = cell) %>%
  # 行顺序与显示中文变量名
  mutate(alc_cat = factor(alc_cat,
                          levels = c("Never","≤2 times/week","≥3 times/week"))) %>%
  arrange(alc_cat)

# χ² 检验
alc_tab <- xtabs(~ is_RPL + alc_cat, data = alc_df)
alc_p <- suppressWarnings(chisq.test(alc_tab)$p.value)

alcohol_table <- alc_cnt %>%
  mutate(变量 = c("饮酒频率，n (%)","", "")) %>%
  relocate(变量, alc_cat, .before = 1) %>%
  rename(`分类` = alc_cat,
         `非RPL组` = `非RPL组`, `RPL组` = `RPL组`) %>%
  mutate(`p值` = c(fmt_p(alc_p), "", ""))

# ==== 2) 连续变量（Wilcoxon） ====
cont_vars <- c(
  "participant.p1488_i0",  # 茶（杯/天）
  "participant.p1498_i0",  # 咖啡（杯/天）
  "participant.p1160_i0",  # 睡眠（小时/天）
  "participant.p884_i0"    # 每周≥10分钟中等强度活动天数
)

clean_cont <- function(x, var){
  x[x %in% c(-1,-3)] <- NA_real_
  if (var %in% c("participant.p1488_i0","participant.p1498_i0")){
    x[x == -10] <- 0.5
  }
  as.numeric(x)
}

df_cont <- df %>%
  mutate(across(all_of(cont_vars), ~clean_cont(., cur_column())))

summarise_one <- function(v, label_cn){
  tmp <- df_cont %>% select(is_RPL, !!sym(v)) %>% drop_na()
  # 组别汇总
  summ <- tmp %>%
    group_by(is_RPL) %>%
    summarise(
      n = n(),
      med = median(.data[[v]], na.rm = TRUE),
      q1  = quantile(.data[[v]], 0.25, na.rm = TRUE),
      q3  = quantile(.data[[v]], 0.75, na.rm = TRUE),
      mn  = mean(.data[[v]], na.rm = TRUE),
      sdv = sd(.data[[v]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(`中位数[IQR]` = sprintf("%.2f [%.2f, %.2f]", med, q1, q3),
           `均值±SD`    = sprintf("%.2f ± %.2f", mn, sdv)) %>%
    select(is_RPL, n, `中位数[IQR]`, `均值±SD`)
  
  # Wilcoxon
  p <- NA_real_
  if (nrow(tmp %>% count(is_RPL)) == 2) {
    p <- suppressWarnings(wilcox.test(tmp[[v]] ~ tmp$is_RPL, exact = FALSE)$p.value)
  }
  
  tibble(
    变量 = label_cn,
    `非RPL组` = paste0(
      summ$n[summ$is_RPL=="非RPL组"], "；",
      summ$`中位数[IQR]`[summ$is_RPL=="非RPL组"],
      "；", summ$`均值±SD`[summ$is_RPL=="非RPL组"]
    ),
    `RPL组` = paste0(
      summ$n[summ$is_RPL=="RPL组"], "；",
      summ$`中位数[IQR]`[summ$is_RPL=="RPL组"],
      "；", summ$`均值±SD`[summ$is_RPL=="RPL组"]
    ),
    `p值` = fmt_p(p)
  )
}

cont_block <- bind_rows(
  summarise_one("participant.p1488_i0", "茶摄入量（杯/天）"),
  summarise_one("participant.p1498_i0", "咖啡摄入量（杯/天）"),
  summarise_one("participant.p1160_i0", "睡眠时间（小时/天）"),
  summarise_one("participant.p884_i0",  "每周≥10分钟中等强度活动天数")
)

# 为了更像论文三线表：把“非RPL组 (n=…) / RPL组 (n=…)”加在列名里
# 这里用各自总体样本量（不剔除缺失）；你若希望“变量特异的非缺失 n”，可改为：对每列单独统计。
n_nonrpl <- sum(df$is_RPL=="非RPL组", na.rm=TRUE)
n_rpl    <- sum(df$is_RPL=="RPL组", na.rm=TRUE)

names(alcohol_table)[names(alcohol_table)=="非RPL组"] <- sprintf("非RPL组 (n=%d)", n_nonrpl)
names(alcohol_table)[names(alcohol_table)=="RPL组"]   <- sprintf("RPL组 (n=%d)", n_rpl)
names(cont_block)[names(cont_block)=="非RPL组"]       <- sprintf("非RPL组 (n=%d)", n_nonrpl)
names(cont_block)[names(cont_block)=="RPL组"]         <- sprintf("RPL组 (n=%d)", n_rpl)

# ==== 3) 合并成一张“论文三线表” ====
three_line_table <- bind_rows(
  alcohol_table,
  tibble(变量=""),   # 空行分隔
  cont_block
)

# 查看 & 导出
print(three_line_table, n = Inf)
write.csv(three_line_table,
          "output/2025.8.11/table_RPL_vs_nonRPL_alcohol_lifestyle.csv",
          row.names = FALSE)

# 如需单独导出饮酒分层计数表，取消下面注释：
write.csv(alcohol_table, "output/2025.8.11/table_alcohol_by_RPL.csv", row.names = FALSE)

