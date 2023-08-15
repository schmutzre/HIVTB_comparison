##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  epitools,
  sjPlot
)

##### data import/ prep ----
ch <- readRDS("data_clean/NOsupress_df.rds") %>% 
  filter(persontime_years.NOsuppression > 0)

#sa <- readRDS("data_clean/art_sa")

test <- ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

combine <- test

##### model ----

### regression model

model_pois <- glm(incidence.nosup ~ cohort, offset=log(persontime_years.NOsuppression), family="poisson", data=combine)
# Run a Poisson regression model with 'cohort' as the predictor and 'incidence.nosup' as the outcome. The log of 'person_time_years' is used as an offset.

# Calculate exponentiated coefficients (IRR) along with CI
exp(coef(summary(model_pois)))

# CI
exp(confint(model_pois))

#double checking manually calculating it

poisson_df <- combine %>% 
  group_by(cohort) %>% 
  summarise(incidence = sum(incidence.nosup == 1), 
            person_years = sum(persontime_years.NOsuppression)/1000) %>% 
  mutate(pois = pois.exact(x = incidence, pt = person_years, conf.level = 0.95)) 

irr <- poisson_df$pois$rate[poisson_df$cohort == "SA"] / poisson_df$pois$rate[poisson_df$cohort == "CH"]

print(irr)

##### main plot ----

plot_IRR <- plot_model(model_pois,
                       colors = "blue",
                       vline.color = "red",
                       show.values = TRUE, 
                       value.offset = .3) +
  theme_bw() +
  theme(plot.title = element_blank())

plot_IRR

ggsave(filename = "results/non_sup/non_sup.IRR.png", plot = plot_IRR, width = 10, height = 8, dpi = 300)

#### assumption ----

## checking for overdispersion - the Poisson model assumes that the mean and variance of the outcome are equal. If the variance is greater than the mean, the data are overdispersed.
disp_test <- dispersiontest(model_pois, trafo=1)  # 'trafo = 1' for log-transformed data
print(disp_test)

#### secondary plots ----

#incidence of viral non-suppression
rna_supression_treshold <- 400
cd4_suppression_treshold <- 350

incidence <- ch %>% 
  summarise(sup_incidence = sum(incidence.nosup), 
          person_years = sum(persontime_years.NOsuppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))


