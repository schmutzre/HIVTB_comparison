##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  sjPlot,
  pROC,
  caret,
  car
)


##### data import ----
table(logistic$cohort)
table(logistic$sex)
table(logistic$region_born)
table(logistic$rna_group)
table(logistic$cd4_group)




logistic_ch <- readRDS("data_clean/art_ch.rds")

logistic_test <- logistic_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m == 1)) %>% 
  ungroup()

logistic_sa <- #...
  
logistic <- logistic_test #... #join them together 

#### model ----

logistic_main <- glmer(case_incident_2m ~ cohort + sex + region_born + age_at_ART_start + rna_group + cd4_group + (1| rna_group), data = logistic, family = "binomial")

#### results ----

# Print model summary
summary(logistic_main)

# Odds Ratios
# Interpretation:
# OR = 1: No association between exposure and outcome.
# OR > 1: The exposure is associated with higher odds of outcome.
# OR < 1: The exposure is associated with lower odds of outcome.
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
#An AUC of 1 means the model has perfect discrimination, while an AUC of 0.5 means the model's discrimination is no better than random chance.

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
