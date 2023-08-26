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
  filter(persontime_years.suppression >= 0)

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
rna_supression_treshold <- 400
cd4_suppression_treshold <- 350

incidence <- ch %>% 
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

# Incidence calculations by each rna_group
incidence_by_rna <- ch %>% 
  group_by(rna_group) %>%
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

# Incidence calculations by each cd4_group
incidence_by_cd4 <- ch %>% 
  group_by(cd4_group) %>%
  summarise(sup_incidence = sum(incidence.sup), 
            person_years = sum(persontime_years.suppression)/1000) %>% 
  mutate(pois = pois.exact(x = sup_incidence, pt = person_years, conf.level = 0.95))

incidence_cd4 <- incidence_by_cd4 %>% 
  ggplot() +
  geom_point(aes(x = cd4_group, y = pois$rate)) +
  geom_errorbar(aes(x= cd4_group, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=incidence$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = incidence$pois$lower, 
                    ymax = incidence$pois$upper), width = 0.2) +
  labs(x= "Baseline CD4 cell count", y = "TB incidence per 1,000 person-years") +
  theme_bw()

incidence_cd4

ggsave(plot = incidence_cd4, filename = "results/incidence/incidence_cd4.png", width = width_descr, height = height_descr)

#rnaplot

incidence_rna <- incidence_by_rna %>% 
  ggplot() +
  geom_point(aes(x = baselineRNA, y = pois$rate)) +
  geom_errorbar(aes(x= baselineRNA, ymin = pois$lower, ymax = pois$upper), width = 0.2) +
  geom_point(aes(x = "overall", y=incidence.total$pois$rate)) +
  geom_errorbar(aes(x= "overall", ymin = incidence.total$pois$lower, 
                    ymax = incidence.total$pois$upper), width = 0.2) +
  labs(x= "Baseline HIV RNA viral load", y = "") +
  scale_y_continuous(limits = c(0,3), expand = c(0,0))+
  theme_bw()

incidence_rna

ggsave(plot = incidence_rna, filename = "results/incidence/incidence_rna.png", width = width_descr, height = height_descr)

incidence.both <- grid.arrange(incidence_cd4, incidence_rna, ncol = 2)

ggsave(plot = incidence.both, file = "results/incidence/incidence_both.png", width = width_descr*1.5, height = height_descr)