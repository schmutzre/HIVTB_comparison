##### Libraries #####

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  ggplot2, # for plotting
  haven,
  lme4,
  sjPlot,
  pROC,
  caret
)

##### data import #####

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### prepare dataframe for analysis #####

logistic_ch <- readRDS("data_clean/df_inc_ch.rds")

logistic_test <- logistic_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m)) %>% 
  ungroup()

logistic_sa <- #...
  
logistic <- logistic_test #... #join them together 

#### analysis

logistic_main <- glmer(case_incident_2m ~ cohort + sex + age_art_start +(1|cd4_group) + (1|rna_group), data = logistic, family = "binomial")

# Print model summary
summary(logistic_main)

# Odds Ratios
exp(fixef(logistic_main))

# Random effects
ranef(logistic_main)

#plot odds ratio
plot_model(logistic_main, type = 'est', title = 'Fixed Effects')

## model vs observed
# Get the predicted probabilities
predicted_probabilities <- predict(logistic_main, type = "response")

# Create the ROC curve
roc_obj <- roc(response = logistic$case_incident_2m, predictor = predicted_probabilities)

# Print the ROC curve
print(roc_obj)

# Plot the ROC curve
plot(roc_obj, print.auc = TRUE)
#It is equivalent to the probability that a randomly chosen positive instance is ranked higher than a randomly chosen negative instance, i.e. it is equivalent to the two sample Wilcoxon rank-sum statistic.

# Calculate the optimal threshold
roc_obj_optimal <- coords(roc_obj, "best")

# The 'best' method uses the Youden's J statistic (sensitivity + specificity - 1)
optimal_threshold <- roc_obj_optimal["threshold"]

# Use this threshold to classify probabilities into binary outcomes
predicted_class <- logistic %>% 
  mutate(prediction = predicted_probabilities,
         classification = as.factor(ifelse(prediction > optimal_threshold$threshold, 1, 0)),
         case_incident_2m = as.factor(case_incident_2m))

# Generate the confusion matrix with the new threshold
confusionMatrix(predicted_class$classification, predicted_class$case_incident_2m)
