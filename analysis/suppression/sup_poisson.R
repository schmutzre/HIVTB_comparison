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
ch <- readRDS("data_clean/supress_df.rds") %>% 
  filter(persontime_years.suppression > 0)

#sa <- readRDS("data_clean/art_sa")

test <- ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

combine <- test

##### model ----

### regression model

model_pois <- glm(incidence.sup ~ cohort, offset=log(persontime_years.suppression), family="poisson", data=combine)
# Run a Poisson regression model with 'cohort' as the predictor and 'incidence.nosup' as the outcome. The log of 'person_time_years' is used as an offset.

# Calculate exponentiated coefficients (IRR) along with CI
exp(coef(summary(model_pois)))

# CI
exp(confint(model_pois))

#double checking manually calculating it

poisson_df <- combine %>% 
  group_by(cohort) %>% 
  summarise(incidence = sum(incidence.sup == 1), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
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

ggsave(filename = "results/sup/sup.IRR.png", plot = plot_IRR, width = 10, height = 8, dpi = 300)

#### assumption ----

## checking for overdispersion - the Poisson model assumes that the mean and variance of the outcome are equal. If the variance is greater than the mean, the data are overdispersed.
disp_test <- dispersiontest(model_pois, trafo=1)  # 'trafo = 1' for log-transformed data
print(disp_test)

#### secondary plots ----

#incidence of viral non-suppression
rna_supression_treshold <- 50

suppression <- ch %>% 
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

# Incidence calculations by each rna_group
suppression_by_rna <- ch %>% 
  filter(rna_group != "NA") %>% 
  group_by(rna_group) %>%
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

# Incidence calculations by each cd4_group
suppression_by_cd4 <- ch %>% 
  filter(cd4_group != "NA") %>% 
  group_by(cd4_group) %>%
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

suppression_cd4 <- suppression_by_cd4 %>% 
  ggplot() +
  geom_point(aes(x = cd4_group, y = pois$rate)) +
  geom_errorbar(aes(x= cd4_group, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=suppression$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = suppression$pois$lower, 
                    ymax = suppression$pois$upper), width = 0.2) +
  labs(x= "Baseline CD4 cell count", y = "Viral suppression per 1,000 person-years") +
  theme_bw()

suppression_cd4

ggsave(plot = incidence_cd4, filename = "results/sup/sup_cd4.png", width = width_descr, height = height_descr)

#rnaplot

suppression_rna <- incidence_by_rna %>% 
  ggplot() +
  geom_point(aes(x = rna_group, y = pois$rate)) +
  geom_errorbar(aes(x= rna_group, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=suppression$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = suppression$pois$lower, 
                    ymax = suppression$pois$upper), width = 0.2) +
  labs(x= "Baseline HIV RNA viral load", y = "Viral suppression per 1,000 person-years") +
  theme_bw()

suppression_rna

ggsave(plot = incidence_rna, filename = "results/sup/sup_rna.png", width = width_descr, height = height_descr)

suppression.both <- grid.arrange(suppression_cd4, suppression_rna, ncol = 2)

ggsave(plot = suppression.both, file = "results/sup/sup_both.png", width = width_descr*1.5, height = height_descr)
