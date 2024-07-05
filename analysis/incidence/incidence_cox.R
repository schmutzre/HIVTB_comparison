#### libraries -----------------------------------------------------------------

library(tidyverse)
library(cmprsk)
library(survival)
library(finalfit)
library(brms)

##### data preparation ---------------------------------------------------------

custom_breaks <- c(16, 34, 44, 100)

df <- readRDS("data_clean/art_noTB.rds")   %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         fup_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           TRUE ~ fup_time),
         fup_years = fup_days / 360,
         cohort = fct_relevel(cohort, "RSA"),
         agegroup = cut(age_at_art_start, breaks = custom_breaks, labels = c("16-34","35-44", "45+"), include.lowest = TRUE),
         bcd4_tr = sqrt(cd4_baseline)) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    death == 1 ~ 2, # only for cases that are incident = 0 (thats how case_when works)
    TRUE ~ 0))) %>%  # Loss to fup is not considered a competing risk
  dplyr::select(id, cohort, art_start_date, incident_tb, fup_years, cd4_group, age_at_art_start, gender, agegroup, cd4_baseline, event_type, bcd4_tr) %>% 
  mutate(sex = fct_drop(gender),
         cd4_group = fct_drop(cd4_group),
         cd4_group = fct_relevel(cd4_group, "350+")) %>% 
  rename(age = age_at_art_start)

df_cox <- df %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  select(-id) 

df_cox_ch <- df %>% 
  filter(cohort == "CH")
  
df_cox_rsa <- df %>% 
  filter(cohort == "RSA")

#### fine-gray model -----------------------------------------------------------

## complete records analysis ##

cov1 <- model.matrix(~factor(cohort) + age +
                       factor(sex) + bcd4_tr,
                     data = df_cox)[, -1]
# Fit the competing risks model
fit_complete <- crr(ftime = df_cox$fup_years, 
                    fstatus = df_cox$event_type, 
                    cov1 = cov1, 
                    failcode = 1)

fit_complete <- fit_complete %>% 
  fit2df(estimate_name = "HR (complete records)", exp = TRUE) 

## Imputed analysis ##

both.imputed <- readRDS("data_clean/imputed/both_imp.rds")

fit_imputed_raw <- with(both.imputed, crr(ftime=fup_years, fstatus=event_type, cov1=cbind(cohort,age, sex, bcd4_tr), failcode=1)) #cov1=cbind()
pooled_fit_imputed  <- pool(fit_imputed_raw)
fit_imputed <- pooled_fit_imputed  %>% 
  fit2df(estimate_name = "HR (multiple imputation)", exp = TRUE)

combined_fit <- cbind(fit_complete, fit_imputed, last_merge = T) %>% select(-explanatory, -last_merge)
combined_fit %>% gt::gt() %>% gt::gtsave("results/incidenceTB/HR.docx")
