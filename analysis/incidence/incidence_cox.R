##### Intro ----
#The Cox model is a multiple regression model for survival data that assumes that the effects of the predictor variables upon survival are multiplicative and do not change over time. It estimates how different factors influence the risk (hazard) of experiencing the event at any given time. The output of the Cox model is a set of hazard ratios, one for each level of each factor included in the model, showing how each factor affects survival risk.
#To summarize, the Kaplan-Meier method describes how survival changes over time, while the Cox model describes how different factors affect survival. While a Kaplan-Meier plot shows raw, unadjusted survival curves for different groups, a Cox model's survival curves show survival probabilities adjusted for other covariates in the model.
#In this specific context, if the cohorts have other differences besides the location (for example, age distribution, sex ratio, etc.), the survival curves from the Cox model would provide a more accurate comparison between the cohorts, because it adjusts for these potential confounders. The Kaplan-Meier curves might show different survival probabilities, but it's hard to tell how much of that is due to location, and how much is due to other differences between the cohorts.
##### Libraries ----

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

##### data import ----

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### preprocessing ----

cox_ch <- readRDS("data_clean/df_inc_ch.rds")

cox_test <- cox_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m)) %>% 
  ungroup()

cox_sa <- #...
  
cox <- cox_test #... #join them together 

#### model ----

# Create a survival object with your time and event variables
cox.surv_obj <- Surv(cox$persontime_years, cox$case_incident_2m)

#Model
cox_model <- coxph(cox.surv_obj ~ cohort, data = cox)

#### results/plotting ----
summary(cox_model)

# Compute predicted survival probabilities for the Cox model
cox_surv <- survfit(cox_model, newdata = cohort_df)

 # Create the new data  
cohort_df <- data.frame(cohort = as.factor(c("CH", "SA")))

surv_plot.cox <- ggsurvplot(
  cox_surv,
  data = cohort_df,
  conf.int = TRUE,
  legend.labs = c("Switzerland", "South Africa"),
  palette = c("#FFA07A", "#B0C4DE"),
  ggtheme = theme_bw(),
  xlab = "Time",
  ylab = "Probability of No Incident"
)

print(surv_plot.cox)

zoom_plot.cox <- ggsurvplot(
  cox_surv,
  data = cohort_df,
  conf.int = TRUE,
  legend.labs = c("Switzerland", "South Africa"),
  palette = c("#FFA07A", "#B0C4DE"),
  ggtheme = theme_bw(),
  xlab = "Time",
  ylab = "Probability of No Incident",
  ylim = c(0.95,1)
)

print(zoom_plot.cox)

ggsave(filename = "results/cox_full.png", plot = surv_plot.cox$plot, width = 10, height = 8, dpi = 300)
ggsave(filename = "results/cox_zoom.png", plot = zoom_plot.cox$plot, width = 10, height = 8, dpi = 300)

##### assumptions ----

#A significant p-value (less than 0.05) indicates that the proportional hazards assumption has been violated
#the smooth line in the plot should be more or less a horizontal line y = 0. 
cox.fit <- cox.zph(cox_model)
print(cox.fit)
plot(cox.fit)



