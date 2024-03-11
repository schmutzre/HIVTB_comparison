##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load( # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  sjPlot,
  pROC, #ROC curve
  caret, #confusionmatrix
  car,
  dplyr,
  reshape2
)


##### data import / preprocessing ----

custom_breaks <- c(16, 24, 34, 44, 100)

ch <- readRDS("data_clean/supress_df.rds") %>% 
  mutate(agegroup = cut(age_at_ART_start, breaks = custom_breaks, include.lowest = TRUE),
         agegroup = as.factor(agegroup),
         incidence = as.factor(case_incident_2m),
         baselineCD4 = as.factor(cd4_group),
         baselineRNA = as.factor(rna_group)) %>% 
  filter(baselineRNA != "NA",
         baselineCD4 != "NA") %>% 
  dplyr::select(id, incidence.sup , sex, cohort, agegroup, baselineCD4, baselineRNA)


sa <- #...
  
logistic <- ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) #... #join them together 

#### model ----

#logistic_main <- glmer(incidence ~ cohort + sex + born + agegroup + baselineRNA + baselineCD4 + (1| baselineRNA), data = logistic, family = "binomial")

model.log <- glm(incidence.sup ~ cohort + sex + agegroup + baselineRNA + baselineCD4, 
                 data = logistic, 
                 family = "binomial")

#### Results ----

# Print model summary
summary(model.log)

# Odds Ratios
# Interpretation:
# OR = 1: No association between exposure and outcome.
# OR > 1: The exposure is associated with higher odds of outcome.
# OR < 1: The exposure is associated with lower odds of outcome.
# Odds Ratios
or <- exp(coef(model.log))

# Confidence Intervals
ci <- exp(confint(model.log))

ref_categories <- sapply(logistic[, c('cohort', 'sex', 'agegroup', 'baselineRNA', 'baselineCD4')], function(x) levels(x)[1])
print(ref_categories)

n <- length(ref_categories)

reference_data <- data.frame(
  OR = rep(1, n),
  LowerCI = rep(1, n),
  UpperCI = rep(1, n),
  row.names = ref_categories
)

# Combining results
results <- data.frame(OR = round(or,3), 
                      LowerCI = round(ci[,1],3), 
                      UpperCI = round(ci[,2],3))

final_results <- rbind(reference_data, results)

#plot odds ratio
plot <- plot_model(model.log, 
                   vline.color = "red",
                   show.values = TRUE, 
                   value.offset = .3,
                   group.terms = c(1, 2, 3, 3, 3, 4, 4, 5, 5)) +
  theme_bw()+
  labs(title = "Risk factors for viral suppression") + 
  theme(plot.title = element_text(hjust = 0.5))
plot

ggsave(plot = plot, filename = "results/sup/sup-odds.png")

#### Diagnostics ----

#reference categories:
levels(logistic$agegroup)
levels(logistic$baselineCD4)
levels(logistic$baselineRNA)

# Get the predicted probabilities from the glm model
predicted_probabilities <- predict(model.log, type = "response", newdata = logistic)

# Create the ROC curve
roc_obj <- roc(response = logistic$incidence.sup, predictor = predicted_probabilities)  # changed `case_incident_2m` to `incidence` 

# Print the ROC curve
print(roc_obj)

# Plot the ROC curve
plotROC <- plot(roc_obj, print.auc = TRUE)

# Calculate the optimal threshold
roc_obj_optimal <- coords(roc_obj, "best")

# Use this threshold to classify probabilities into binary outcomes
predicted_class <- logistic %>% 
  mutate(prediction = predicted_probabilities,
         classification = as.factor(ifelse(prediction > as.numeric(roc_obj_optimal["threshold"]), 1, 0)),
         incidence = as.factor(incidence.sup))

# Generate the confusion matrix with the new threshold
confusionMatrix(predicted_class$classification, predicted_class$incidence, positive = "1")

#### Plots (diagnostic) ----
## Confusion matrix

# Generate the confusion matrix
cm <- confusionMatrix(predicted_class$classification, predicted_class$incidence, positive = "1")

# Convert the confusion matrix to a table format
cm_table <- as.table(cm$table)

# Melt the table for ggplot
cm_melted <- melt(cm_table) %>% 
  mutate(Reference = as.factor(Reference),
         Prediction = as.factor(Prediction))

# Plot using ggplot2
plot_confusion <- ggplot(data=cm_melted, aes(x=Prediction, y=Reference)) +
  geom_tile(aes(fill=value), color='white') +
  geom_text(aes(label=sprintf("%d", value)), vjust=1) +
  scale_fill_gradient(low="white", high="blue") +
  scale_x_discrete(limits=c("0", "1")) +   # Set discrete x-axis scale
  scale_y_discrete(limits=c("0", "1")) +   # Set discrete y-axis scale
  theme_bw() +
  labs(title="Confusion Matrix", x="Prediction", y="Observed", fill="Frequency")+
  theme(plot.title = element_text(hjust = 0.5))


## Key metrics
cm2 <- confusionMatrix(predicted_class$classification, predicted_class$incidence, positive = "1")
key_metrics <- c(Accuracy = cm2$overall['Accuracy'],
                 Sensitivity = cm2$byClass['Sensitivity'],
                 Specificity = cm2$byClass['Specificity'],
                 PPV = cm2$byClass['Pos Pred Value'],
                 NPV = cm2$byClass['Neg Pred Value'])

# Convert to data frame for ggplot
key_metrics_df <- data.frame(Metric = names(key_metrics), Value = as.numeric(key_metrics))

key_metrics_df$Metric <- sub("\\..*", "", key_metrics_df$Metric)

# Plot
plot_metrics <- ggplot(key_metrics_df, aes(x = Metric, y = Value)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  geom_text(aes(label=sprintf("%.2f", Value)), vjust=2, color = "white") +
  ylim(0, 1) +
  labs(title="Key Metrics from Confusion Matrix", y="", x ="") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))
