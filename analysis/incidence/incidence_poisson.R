##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  epitools
  )

##### data import ----

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### preprocessing ----

poisson_ch <- readRDS("data_clean/df_inc_ch.rds")

poisson_test <- poisson_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA"))) %>% 
  group_by(cohort) %>% 
  mutate(cases_p_cohort = sum(case_incident_2m)) %>% 
  ungroup()

poisson_sa <- #...

poisson <- poisson_test #... #join them together 

##### model ----

### first model

poisson_df <- poisson %>% 
  group_by(cohort) %>% 
  summarise(tb_incidence = sum(case_incident_2m == 1), 
            person_years = sum(persontime_years)/100000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95)) 

irr <- poisson_df$pois$rate[poisson_df$cohort == "SA"] / poisson_df$pois$rate[poisson_df$cohort == "CH"]

print(irr)

### regression model

model_main <- glm(case_incident_2m ~ cohort, offset=log(persontime_years), family="poisson", data=poisson)
# Run a Poisson regression model with 'cohort' as the predictor and 'case_incident_2m' as the outcome. The log of 'person_time_years' is used as an offset.

# Calculate exponentiated coefficients (IRR) along with CI
exp(coef(summary(model_main)))

# CI
exp(confint(model_main))

## modeled outcome vs observed outcome
# add predictions to data

poisson_compare <- poisson %>%
  mutate(predicted_counts_rate = predict(model_main, type = "response"),
         predicted_counts_abs = format((predicted_counts_rate * persontime_years), scientific = FALSE))

ggplot(poisson_compare, aes(x = case_incident_2m, y = predicted_counts_abs)) +
  geom_jitter() +
  labs(x = "Observed", y = "Predicted", title = "Observed vs Predicted")

#add 2x2 matrix

##### plots ----

### For main analysis
results_main <- broom::tidy(model_main, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(model = "main",
         term = factor(term, levels = c("cohort"), labels = c("Cohort"))) %>% 
  mutate(term = c("Intercept", "Effect"))

# Create the plot
incidence_rate_ratio <- ggplot(results_main %>% filter(term == "Effect"), aes(x = estimate, y = "")) +
  annotate("rect", xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#FFA07A", alpha = 0.3) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#B0C4DE", alpha = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkblue", linewidth = 1) +
  geom_point(size = 2, color = "black") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.1, linewidth = 1,  color = "black") +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5,  3), limits = c(0, 3), expand = c(0,0)) +
  geom_text(aes(x = 0.5, y = 1.5, label = "Switzerland"), hjust = 0.5, vjust = 1, size = 5, fontface = "bold") +  
  geom_text(aes(x = 2, y = 1.5, label = "South Africa"), hjust = 0.5, vjust = 1, size = 5, fontface = "bold") +  
  labs(x = "Incidence Rate Ratio (IRR)", y = "") +
  theme_bw()

 incidence_rate_ratio

# Save plot as a jpeg (check dimensions again after plotting)
ggsave(filename = "results/incidence_primary.png", plot = incidence_rate_ratio, width = 10, height = 8, dpi = 300)

#additional plots concerning incidence
inc_ch <- ch %>% 
  filter(case_incident_2m == 1) %>% 
  group_by(rna_group, cd4_group) %>% 
  summarise(n = n())

inc_n_ch <- flextable(inc_ch)

# Save your table as an image
save_as_image(inc_n_ch, "results/inc_n_ch.png", expand = 10)

#### assumption ----

## checking for overdispersion - the Poisson model assumes that the mean and variance of the outcome are equal. If the variance is greater than the mean, the data are overdispersed.
disp_test <- dispersiontest(model_main, trafo=1)  # 'trafo = 1' for log-transformed data
print(disp_test)


