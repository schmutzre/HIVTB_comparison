library(dplyr)

#### interval censored non-parametric ------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds") %>% select(id, date_cd4, cd4, time_diff)
cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% select(id, date_cd4, cd4, time_diff)
cd4 <- rbind(cd4_ch, cd4_rsa) %>% mutate(time_diff = ifelse(time_diff < 0, 0, time_diff))
art_both <- readRDS("data_clean/art.rds") %>% select(id, art_start_date, fup_time, cohort)
cd4 <- left_join(cd4, art_both, by = "id") %>% mutate(recovered = ifelse(cd4 >350,1,0))

cd4_interval <- cd4 %>%
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
fit <- ic_np(cbind(left, right) ~ cohort, data = cd4_interval %>% sample_frac(.1))
summary(fit)
png("results/plot_turnbull.png", width = 16, height = 12, units = "cm", res = 300)
plot(fit)  # Recreate the plot
dev.off()

#### non-interval censored non-parametric ------------------------------------------

cd4_right <- cd4 %>%
  group_by(id) %>%
  arrange(time_diff, .by_group = TRUE) %>%
  filter(if (any(recovered == 1)) row_number() == which.max(recovered == 1) else row_number() == 1) %>%
  ungroup() %>% 
  mutate(time = ifelse(recovered == 1, time_diff, fup_time)) %>%
  select(id, cohort, time, recovered)

fit_km <- survival::survfit(Surv(time, recovered) ~ cohort, data = cd4_right %>% sample_frac(.1))
plot(fit_km, ylim = c(0,1))
