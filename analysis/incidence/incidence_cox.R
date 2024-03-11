#### libraries -----------------------------------------------------------------

install.packages("tidyverse")
remove.packages("tidyverse")
remove.packages("ggplot2")

library(tidyverse)
library(cmprsk)
library(survival)
library(finalfit)

##### data preparation ---------------------------------------------------------

custom_breaks <- c(16, 34, 44, 100)

df <- readRDS("data_clean/art_noTB.rds")  %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           incident_tb == 0 ~ last_persontime),
         persontime_years = persontime_days / 360,
         persontime_death_days = case_when(
           !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
           is.na(exitdate) ~ last_persontime),
         cohort = fct_relevel(cohort, "RSA"),
         agegroup = cut(age_at_art_start, breaks = custom_breaks, labels = c("16-34","35-44", "45+"), include.lowest = TRUE)) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~ 0 # Loss to follow-up isnt considered a competing risk
  ))) %>% 
  dplyr::select(cohort, art_start_date, incident_tb, persontime_years, cd4_group, age_at_art_start, gender, agegroup, cd4_baseline, event_type) %>% 
  filter(persontime_years > 0) %>% 
  filter_all(all_vars(!is.na(.))) %>% 
  mutate(gender = fct_drop(gender),
         cd4_group = fct_drop(cd4_group),
         cd4_group = fct_relevel(cd4_group, "350+"))

#### fine-gray model -----------------------------------------------------------

explanatory   <- c("cohort","agegroup", "gender", "cd4_group")
dependent_crr <- "Surv(persontime_years, event_type)"

df %>%
  # Summary table
  summary_factorlist(dependent_crr, explanatory, 
                     column = TRUE, fit_id = TRUE) %>% 
  # Fine and Gray competing risks regression
  ff_merge(
    df %>%
      crrmulti(dependent_crr, explanatory) %>%
      fit2df(estimate_suffix = " (competing risks multivariable)")
  ) %>% 
  select(-fit_id, -index) %>%
  dependent_label(df, "Incident TB") 