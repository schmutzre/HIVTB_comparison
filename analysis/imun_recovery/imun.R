##### set-up -------------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  lmerTest,
  mgcv,
  gridExtra, 
  wesanderson,
  ggpubr,
  survival,
  survminer,
  scam, 
  gamm4
)

source("utils/plot.R")

# Set seed for reproducibility

set.seed(123)

#### data ---------------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH")) %>% 
  select(-region)

cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4 <- rbind(cd4_ch %>% dplyr::select(-timepoint),
             cd4_rsa %>% dplyr::select(-timepoint)) %>% 
  filter(time_diff >= -60 & time_diff < 360) %>% 
  mutate(
    cd4_trans = sqrt(cd4))

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

#### model cd4 (interaction) ---------------------------------------------------

cd4$interaction_term <- interaction(cd4$cohort, cd4$presenting_tb)

model_cd4_count <- gamm4(cd4_trans ~ cohort + presenting_tb  +
                           s(time_diff, by = interaction_term, bs = "cr", k = 4),
                                            random = ~ (1 | id),
                                            data = cd4)
summary(model_cd4_count$gam)

pred <- pred_trend(model_cd4_count, data = cd4)

trend_cd4 <- ggplot(pred, aes(x = time_diff)) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 
  
trend_cd4

ggsave(trend_cd4, filename = "results/cd4_predEff.png",
       width = 20, height = 12.7, unit = "cm")

#### model cd4 (separate for presenting TB) -------------------------------------

#non-presenting
model_cd4_count_npres <- gamm4(cd4_trans ~ cohort +
                                s(time_diff, by = cohort, bs = "cr", k = 4),
                              random = ~ (1 | id),
                              data = cd4_npres)

pred_npres <- pred_trend_uni(model_cd4_count_npres, 
                   data = cd4_npres)

ggplot(pred_npres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 600)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

#presenting
model_cd4_count_pres <- gamm4(cd4_trans ~ cohort +
                           s(time_diff, by = cohort, bs = "cr", k = 4),
                         random = ~ (1 | id),
                         data = cd4_pres)

pred_pres <- pred_trend_uni(model_cd4_count_pres, 
                         data = cd4_pres)

ggplot(pred_pres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

#### simpler model

# Load necessary library
library(lme4)

# Fit a linear mixed model with the polynomial terms
model_cd4_count_simple <- lmer(cd4_trans ~ cohort + presenting_tb + 
                                 time_diff + I(time_diff^2) +  I(time_diff^3)+
                                 (1 | id),
                               data = cd4)

# Check the summary of the model
summary(model_cd4_count_simple)


# Load necessary libraries
library(lme4)
library(merTools)
library(ggplot2)

# Assuming you have already fitted the model: model_cd4_count_simple

# Prepare a new data frame for predictions
# This frame should cover the range of time_diff you are interested in
new_data <- expand.grid(
  time_diff = seq(min(cd4$time_diff), max(cd4$time_diff), length.out = 20),  # Increase length.out
  cohort = levels(cd4$cohort),
  presenting_tb = levels(cd4$presenting_tb),
  id = factor(1)  # Assuming id is a factor in your original model
)



# Add polynomial terms
new_data$time_diff_squared <- new_data$time_diff^2
new_data$time_diff_cubed <- new_data$time_diff^3

# Get predictions and confidence intervals
library(merTools)
predictions <- predictInterval(model_cd4_count_simple, newdata = new_data, level = 0.95, type = "linear.prediction")
# Combine predictions with new_data
new_data$fit <- predictions
new_data$lower <- predictions$lwr
new_data$upper <- predictions$upr

# Plot the predictions with confidence intervals
ggplot(new_data, aes(x = time_diff, y = fit$fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) + # Confidence interval
  geom_line() +  # Fitted line
  labs(x = "Time Difference", y = "CD4 Trans", title = "Trend in CD4 Trans with Confidence Intervals") +
  theme_minimal()