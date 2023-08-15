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

##### data import/ prep ----

custom_breaks <- c(16, 24, 34, 44, 100)

ch <- readRDS("data_clean/art_ch.rds") %>% 
  mutate(agegroup = cut(age_at_ART_start, breaks = custom_breaks, include.lowest = TRUE),
         agegroup = as.factor(agegroup),
         incidence = case_incident_2m,
         baselineCD4 = as.factor(cd4_group),
         baselineRNA = as.factor(rna_group),
         born = as.factor(region_born),
         persontime_years = case_when(
           case_incident_2m == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")/360),
           case_incident_2m == 0 ~ last_persontime/360
         )
         ) %>% 
  filter(persontime_years > 0) %>% 
  dplyr::select(id, incidence, sex, cohort, born, agegroup, baselineCD4, baselineRNA, persontime_years)


#sa <- readRDS("data_clean/art_sa")

##### preprocessing ----

poisson_test <- ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

poisson_sa <- #...

poisson <- poisson_test #... #join them together 

##### model ----

### regression model

model_pois <- glm(incidence ~ cohort, offset=log(persontime_years), family="poisson", data=poisson)
# Run a Poisson regression model with 'cohort' as the predictor and 'case_incident_2m' as the outcome. The log of 'person_time_years' is used as an offset.

# Calculate exponentiated coefficients (IRR) along with CI
exp(coef(summary(model_pois)))

# CI
exp(confint(model_pois))

#double checking manually calculating it

poisson_df <- poisson %>% 
  group_by(cohort) %>% 
  summarise(tb_incidence = sum(incidence == 1), 
            person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95)) 

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


ggsave(filename = "results/incidence/incidence_IRR.png", plot = plot_IRR, width = 10, height = 8, dpi = 300)

#### assumption ----

## checking for overdispersion - the Poisson model assumes that the mean and variance of the outcome are equal. If the variance is greater than the mean, the data are overdispersed.
disp_test <- dispersiontest(model_main, trafo=1)  # 'trafo = 1' for log-transformed data
print(disp_test)

#### secondary plots ----



