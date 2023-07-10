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
cox_model <- coxph(cox.surv_obj ~ cohort + age_art_start + sex, data = cox)

#### results/plotting ----
#The exponentiated coefficients in the second column of the first panel (and in the first column of the second panel) of the output are interpretable as multiplicative effects on the hazard. 
#Thus, for example, holding the other covariates constant, an additional year of age reduces the yearly hazard of incident TB by a factor of eb2 = 0.9956 on average.
sum_cox <- summary(cox_model)

# Find the row corresponding to the predictor of interest
row_index <- which(rownames(sum_cox$coefficients) == "cohortSA")

# Extract the exponentiated coefficient (hazard ratio)
exp_coef <- sum_cox$conf.int[row_index, "exp(coef)"]
lower_ci <- sum_cox$conf.int[row_index, "lower .95"]
upper_ci <- sum_cox$conf.int[row_index, "upper .95"]

# Create a data frame for plotting
plot_data <- data.frame(
  Predictor = "cohort",
  Hazard_Ratio = exp_coef,
  Lower_CI = lower_ci,
  Upper_CI = upper_ci
)

# Create the plot
hazard_ratio <- ggplot(plot_data, aes(x = Hazard_Ratio, y = "", xmin = Lower_CI, xmax = Upper_CI)) +
  geom_point(size = 2, color = "black") +
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#FFA07A", alpha = 0.3) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#B0C4DE", alpha = 0.3) +
  geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.1, linewidth = 1,  color = "black") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5,  3), limits = c(0, 3), expand = c(0,0)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkblue", linewidth = 1) +  # Add a dashed line at HR = 1 (no effect)
  geom_text(aes(x = 0.5, y = 1.5, label = "Switzerland"), hjust = 0.5, vjust = 1, size = 5, fontface = "bold") +  
  geom_text(aes(x = 2, y = 1.5, label = "South Africa"), hjust = 0.5, vjust = 1, size = 5, fontface = "bold") +
  theme_bw() +
  labs(x = "Hazard ratio", y = "")

# Print the plot
print(hazard_ratio)

ggsave(filename = "results/hazard_ratio.png", plot = hazard_ratio, width = 10, height = 8, dpi = 300)

##### assumptions ----

#A significant p-value (less than 0.05) indicates that the proportional hazards assumption has been violated
#the smooth line in the plot should be more or less a horizontal line y = 0. 
cox.fit <- cox.zph(cox_model)
print(cox.fit)
plot(cox.fit)



