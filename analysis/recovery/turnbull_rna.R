library(dplyr)

#### interval censored non-parametric ------------------------------------------

rna_ch <- readRDS("data_clean/ch/rna_ch.rds") %>% select(id, date_rna, rna, time_diff)
rna_rsa <- readRDS("data_clean/rsa/rna_rsa.rds") %>% select(id, date_rna, rna, time_diff)
rna <- rbind(rna_ch, rna_rsa) %>% mutate(time_diff = ifelse(time_diff < 0, 0, time_diff))
art_both <- readRDS("data_clean/art.rds") %>% select(id, art_start_date, fup_time, cohort)
rna <- left_join(rna, art_both, by = "id") %>% mutate(recovered = ifelse(rna <400,1,0))

rna_interval <- rna %>%
  group_by(id) %>%
  arrange(time_diff, .by_group = TRUE) %>%
  mutate(left = lag(time_diff, default = 0),
         right = if_else(row_number() == n(), NA, time_diff),
         censoring = if_else(row_number() == n(), 1, 0)) %>%
  ungroup() %>% 
  select(id, cohort, left, right, censoring, recovered) %>% 
  group_by(id) %>%
  # Summarize to keep only one row per id
  summarize(across(c(cohort, left, right, censoring, recovered), 
                   ~ if (any(recovered == 1)) {
                     .[which.max(recovered == 1)]
                   } else {
                     last(.)
                   },
                   .names = "{.col}")) %>%
  ungroup()

library(icenReg)
fit <- ic_np(cbind(left, right) ~ cohort, data = rna_interval %>% sample_frac(.1))
summary(fit)
png("results/plot_turnbull_rna.png", width = 16, height = 12, units = "cm", res = 300)
plot(fit)  # Recreate the plot
dev.off()

#### non-interval censored non-parametric ------------------------------------------

rna_right <- rna %>%
  group_by(id) %>%
  arrange(time_diff, .by_group = TRUE) %>%
  filter(if (any(recovered == 1)) row_number() == which.max(recovered == 1) else row_number() == 1) %>%
  ungroup() %>% 
  mutate(time = ifelse(recovered == 1, time_diff, fup_time)) %>%
  select(id, cohort, time, recovered)

fit_km <- survival::survfit(Surv(time, recovered) ~ cohort, data = rna_right %>% sample_frac(.1))
plot(fit_km, ylim = c(0,1), )