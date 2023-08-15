##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  lubridate, # for date handling
  rlang, # for '!!' inside the functions
  AER, # for dispersiontest
  ggplot2, # for plotting
  haven,
  epitools,
  sjPlot
  )

##### data import/ prep ----
#SHCS
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
  dplyr::select(id, incidence, sex, cohort, born, agegroup, baselineCD4, baselineRNA, persontime_years, date_tb, art_start_date)

#rate 
pois_ch <- ch %>% 
  summarise(tb_incidence = sum(incidence == 1), 
            person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95))

#Western Cape

#sa <- readRDS("data_clean/art_sa")

##### preprocessing ----

poisson_test <- ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

poisson_sa <- #...

  
## Complete dataset
  
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

# Overall TB incidence rate
indidence.total <- ch %>% 
  summarise(tb_incidence = sum(incidence == 1), 
            person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95)) 

# Incidence calculations by each rna_group
incidence_by_rna <- ch %>% 
  group_by(baselineRNA) %>%
  summarise(tb_incidence = sum(incidence == 1), 
            person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95))

# Incidence calculations by each cd4_group
incidence_by_cd4 <- ch %>% 
  group_by(baselineCD4) %>%
  summarise(tb_incidence = sum(incidence == 1), 
            person_years = sum(persontime_years)/1000) %>% 
  mutate(pois = pois.exact(x = tb_incidence, pt = person_years, conf.level = 0.95))

#cd4plot

incidence_cd4 <- incidence_by_cd4 %>% 
  ggplot() +
  geom_point(aes(x = baselineCD4, y = pois$rate)) +
  geom_errorbar(aes(x= baselineCD4, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=overall_tb_incidence$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = overall_tb_incidence$pois$lower, 
                    ymax = overall_tb_incidence$pois$upper), width = 0.2) +
  labs(x= "Baseline CD4 cell count", y = "TB incidence per 1,000 person-years") +
  theme_bw()

incidence_cd4

ggsave(plot = incidence_cd4, filename = "results/incidence/incidence_cd4.png", width = width_descr, height = height_descr)

#rnaplot

incidence_rna <- incidence_by_rna %>% 
  ggplot() +
  geom_point(aes(x = baselineRNA, y = pois$rate)) +
  geom_errorbar(aes(x= baselineRNA, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=overall_tb_incidence$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = overall_tb_incidence$pois$lower, 
                    ymax = overall_tb_incidence$pois$upper), width = 0.2) +
  labs(x= "Baseline HIV RNA viral load", y = "TB incidence per 1,000 person-years") +
  theme_bw()

incidence_rna

ggsave(plot = incidence_rna, filename = "results/incidence/incidence_rna.png", width = width_descr, height = height_descr)

#time to tb
time_to_tb <- ch %>% 
  filter(incidence ==1) %>% 
  mutate(time_to_tb = as.numeric(date_tb - (art_start_date))) %>% 
  summarise(median = median(time_to_tb, na.rm = TRUE),
            "25%" = quantile(time_to_tb, 0.25, na.rm = TRUE),
            "75%" = quantile(time_to_tb, 0.75, na.rm = TRUE))
time_to_tb

