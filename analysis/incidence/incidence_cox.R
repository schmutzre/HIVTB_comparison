##### Intro ----
#The Cox model is a multiple regression model for survival data that assumes that the effects of the predictor variables upon survival are multiplicative and do not change over time. It estimates how different factors influence the risk (hazard) of experiencing the event at any given time. The output of the Cox model is a set of hazard ratios, one for each level of each factor included in the model, showing how each factor affects survival risk.
#To summarize, the Kaplan-Meier method describes how survival changes over time, while the Cox model describes how different factors affect survival. While a Kaplan-Meier plot shows raw, unadjusted survival curves for different groups, a Cox model's survival curves show survival probabilities adjusted for other covariates in the model.
#In this specific context, if the cohorts have other differences besides the location (for example, age distribution, sex ratio, etc.), the survival curves from the Cox model would provide a more accurate comparison between the cohorts, because it adjusts for these potential confounders. The Kaplan-Meier curves might show different survival probabilities, but it's hard to tell how much of that is due to location, and how much is due to other differences between the cohorts.

#### libraries -----------------------------------------------------------------

library(tidyverse)
library(cmprsk)
library(survival)

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  haven,
  survival,
  gridExtra,
  cmprsk,
  sjPlot
)
library(sjPlot)
##### data preparation ---------------------------------------------------------

custom_breaks <- c(16, 24, 34, 44, 100)

df <- readRDS("data_clean/art_noTB.rds") %>% 
  mutate(incident_tb = as.numeric(as.character(incident_tb)),
         persontime_days = case_when(
           incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
           incident_tb == 0 ~ last_persontime),
         persontime_years = persontime_days / 360,
         persontime_death_days = case_when(
           !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
           is.na(exitdate) ~ last_persontime),
         cohort = fct_relevel(cohort, "RSA")) %>% 
  dplyr::select(cohort, art_start_date, incident_tb, persontime_years, cd4_group, age_at_art_start, gender) %>% 
  filter(persontime_years > 0) %>% 
  mutate(event_type = as.factor(case_when(
    incident_tb == 1 ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~ 0 # Loss to follow-up isnt considered a competing risk
  )))

#### model ---------------------------------------------------------------------

surv_obj <- Surv(df$persontime_years, 
                     df$incident_tb)

model <- coxph(cox.surv_obj ~ cohort + age_at_art_start + gender + cd4_group, 
               data = df)

#### results -------------------------------------------------------------------

summary(model)

plot_model(model)
hr <- exp(coef(model.cox))
ci <- exp(confint(model.cox))

#plot hazard ratio
plot <- plot_model(model, 
                   vline.color = "red",
                   show.values = TRUE, 
                   value.offset = .3,
                   group.terms = c(1, 2, 3, 4, 4, 4)) +
  theme_bw()+
  labs(title = "Hazard Ratios for Factors Associated with TB Incidence",
       y = "Hazard ratios") + 
  theme(plot.title = element_text(hjust = 0.5))

plot

ggsave(plot = plot, filename = "results/incidence/hazard.png")

##### assumptions ----

#A significant p-value (less than 0.05) indicates that the proportional hazards assumption has been violated
#the smooth line in the plot should be more or less a horizontal line y = 0. 
cox.fit <- cox.zph(model.cox)
print(cox.fit)
plot(cox.fit)

#### Competing risk model (Fine-Gray) ----
# Create factor with three levels

gray <- cox %>% 
  filter(persontime_years > 0)

# Define factors for model spec
factors <- c("cohort", "sex", "agegroup", "baselineRNA", "baselineCD4")
gray$sex <- relevel(gray$sex, ref = "Male")

model_spec <- reformulate(
  paste0("cohort+", paste(factors, collapse = "+"))
) # ~age + sex + income

covariates_matrix <- model.matrix(
  model_spec,
  data = gray,
  contrasts.arg = lapply(gray[factors], contrasts)
)[, -1] # first column is constant intercept (1)

head(covariates_matrix, 3)

model.finegray <- crr(ftime = gray$persontime_years, 
                      fstatus = gray$event_type, 
                      cov1 = covariates_matrix,
                      failcode = 1,
                      cencode = 0)

summary(model.finegray)

# Extract coefficients and confidence intervals for subdistribution hazards
coef <- model.finegray$coef
se <- sqrt(diag(model.finegray$var))
HR <- exp(coef)
lower_CI <- exp(coef - 1.96 * se)
upper_CI <- exp(coef + 1.96 * se)

results <- data.frame(HazardRatio = HR, LowerCI = lower_CI, UpperCI = upper_CI)
results$cov <- rownames(results)
rownames(results) <- NULL

results <- results %>%
  mutate(cov = c("South Africa (cohort)", "Female (sex)", "24-34 (age)", "34-44 (age)", "44 + (age)", "1,000 - 9,999 (Baseline VL)", "10,000+ (Baseline VL)", "NA (Baseline VL)", "100 - 349 (Baseline CD4)", "350+ (Baseline CD4)", "NA (Baseline CD4)"))

results$cov <- factor(results$cov, levels = rev(c("South Africa (cohort)", "Female (sex)", "24-34 (age)", "34-44 (age)", "44 + (age)", "1,000 - 9,999 (Baseline VL)", "10,000+ (Baseline VL)", "NA (Baseline VL)", "100 - 349 (Baseline CD4)", "350+ (Baseline CD4)", "NA (Baseline CD4)")))

plot.gray <- results %>% 
  ggplot(aes(x = HazardRatio, y = cov, color = HazardRatio > 1)) +
  geom_vline(xintercept = 1, color = "gray75") +
  geom_linerange(aes(xmin = lower_CI, xmax = upper_CI), size = 1, alpha = 0.5) +
  geom_point(size = 2) +
  geom_text(aes(label = sprintf("%.2f", HazardRatio)), vjust = -0.9, size = 3) +
  theme_bw() +
  scale_color_manual(values = c("green4", "red3"), guide = "none") +
  labs(title = "Estimated HR using a fine-gray model", y = NULL,
       x = "Hazard ratio estimate (95% CI)") +
  theme(axis.text.y = element_text(hjust = 0, size = 10),
        plot.title = element_text(hjust = 0.5))  # This line centers the title

plot.gray

ggsave(plot = plot.gray, filename = "results/incidence/hazard-gray.png")

