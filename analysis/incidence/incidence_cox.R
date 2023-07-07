##### Libraries #####

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  survival,
  gridExtra
)

##### data import #####

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### prepare dataframe for analysis #####

cox_ch <- readRDS("data_clean/df_inc_ch.rds")

cox_test <- cox_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m)) %>% 
  ungroup()

cox_sa <- #...
  
cox <- cox_test #... #join them together 

# Create a survival object with your time and event variables
cox.surv_obj <- Surv(cox$persontime_years, cox$case_incident_2m)

#Model
cox_model <- coxph(cox.surv_obj ~ cohort, data = cox)

#Results
summary(cox_model)

#Check model assumptions; A significant p-value (less than 0.05) indicates that the proportional hazards assumption has been violated
cox.zph(cox_model)


