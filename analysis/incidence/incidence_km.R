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
  survminer,
  gridExtra
)

##### data import ----

kaplan_ch <- readRDS("data_clean/df_inc_ch.rds")

kaplan_test <- kaplan_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m)) %>% 
  ungroup()

kaplan_sa <- #...
  
kaplan <- kaplan_test #... #join them together 

#### model (inc. results/plotting) ----

# Create a survival object with your time and event variables
km.surv_obj <- Surv(kaplan$persontime_years, kaplan$case_incident_2m)

# Use the survfit() function to calculate survival probabilities
km.fit <- survfit(km.surv_obj ~ kaplan$cohort)  

# Generate Kaplan-Meier plot
kaplan_plot <- ggsurvplot(
  km.fit, 
  data = kaplan, 
  pval = TRUE, 
  xlab = "Time after starting ART in years",
  ylab = "Probability of No Incident",
  risk.table = TRUE,
  conf.int = TRUE,
  palette = c("#FFA07A", "#B0C4DE"),
  legend.labs = 
    c("Switzerland", "South Africa")) 

print(kaplan_plot)

# Create a zoomed in plot
zoom_plot <- kaplan_plot$plot + 
  coord_cartesian(ylim = c(0.95, 1)) +
  labs(y = "Probability of No Incident") +
  theme_bw()

print(zoom_plot)

km_full <- grid.arrange(kaplan_plot$plot, kaplan_plot$table, ncol = 1)

ggsave(filename = "results/km_full.png", plot = kaplan_plot$plot, width = 10, height = 8, dpi = 300)
ggsave(filename = "results/kmrisk_full.png", plot = km_full, width = 10, height = 8, dpi = 300)
ggsave(filename = "results/km_zoom.png", plot = zoom_plot, width = 10, height = 8, dpi = 300)

#### log-rank test ----

# Run the log-rank test
log_rank_test <- survdiff(km.surv_obj ~ kaplan$cohort)  

# Print the results
print(log_rank_test)

#p<.05 --> significant differents between the two cohorts. 

