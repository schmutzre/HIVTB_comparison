##### Libraries #####

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  brms,
  haven,
  tidyverse,
  bayesplot
)

##### data import #####

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### prepare dataframe for analysis #####

baypoisson_ch <- ch %>% 
  mutate(art_start_date_2m = art_start_date + months(2), # Adding 2 months to the start date of antiretroviral therapy (ART) as incidence is defined as 2 months after art start, so the time in between shoouldn't be counted.
         persontime_years = ifelse(case_incident_2m == 1, # Computing the person-time in years. If the event of interest occurred, it's the time from ART start to the event. 
                                   as.numeric(difftime(date_tb, art_start_date_2m, units = "days")) / 360, #convert to years
                                   ifelse(!is.na(exitdate), # If the event did not occur, it's either the time from ART start to the exit date or to the end of the study.
                                          as.numeric(difftime(exitdate, art_start_date_2m, units = "days")) / 360,  
                                          as.numeric(difftime(as.Date("2022-12-31"), art_start_date_2m, units = "days")) /360
                                   ))) %>% 
  select(id, sex, date_tb, case_incident_2m, cd4_baseline, rna_baseline, cd4_group, rna_group, persontime_years, exitdate, art_start_date, art_start_date_2m) %>% 
  filter(persontime_years > 0) %>% # Excluding patients with less than 2 months follow-up (now only 5005 patients)
  mutate(cohort = as.factor("CH")) # Adding a cohort indicator

baypoisson_test <- poisson_ch %>% 
  mutate(cohort = factor(ifelse(row_number() <= 2000, "CH", "SA")))

baypoisson_sa <- #...
  
baypoisson <- poisson_test #... #join them together 

baypoisson[] <- lapply(poisson, zap_labels)

##### analysis #####

##main model
model_bayes_main <- brm(case_incident_2m ~ cohort + offset(log(persontime_years)), 
                        family = poisson(), data = baypoisson,
                        prior = NULL, 
                        control = list(adapt_delta = 0.99))

# get the posterior_summary
post_summ <- posterior_summary(model_bayes_main) %>% 
  as.data.frame() %>% 
  rownames_to_column("Effect") %>% 
  mutate(Estimate = round(exp(Estimate), 6),
         'Est.Error' = round(exp(`Est.Error`),6),
         'Lower_CI' = round(exp(`Q2.5`),6),
         'Upper_CI' = round(exp(`Q97.5`),6)) %>% 
  select(-c(Q2.5,Q97.5))

# print the updated summary
print(post_summ)

# Calculate the probability that the IRR of `cohortSA` is greater than 1
posterior_samples(model_bayes_main, pars = "cohortSA") %>%
  mutate(cohortSA = exp(cohortSA)) %>%
  summarise(prob = mean(cohortSA > 1))

##hierarchical model

model_bayes_hierarchical <- brm(case_incident_2m ~ cohort + cd4_group + rna_group + 
                                  (1 | cd4_group) + (1 | rna_group) + 
                                  offset(log(persontime_years)), 
                                family = poisson(), data = baypoisson,
                                control = list(adapt_delta = 0.99),
                                prior = NULL)

# get the posterior_summary
post_summ_hie <- posterior_summary(model_bayes_hierarchical) %>% 
  as.data.frame() %>% 
  rownames_to_column("Effect") %>% 
  mutate(Estimate = case_when(
    str_detect(Effect, "^b") ~ format(round(exp(Estimate), 6), scientific = FALSE),
    TRUE ~ format(round(Estimate,6) , scientific = FALSE)
  ),
  Est.Error = case_when(
    str_detect(Effect, "^b") ~ format(round(exp(`Est.Error`), 6), scientific = FALSE),
    TRUE ~ format(round(`Est.Error`,6) , scientific = FALSE)
  ),
  Q2.5 = case_when(
    str_detect(Effect, "^b") ~ format(round(exp(`Q2.5`), 6), scientific = FALSE),
    TRUE ~ format(round(`Q2.5`,6) , scientific = FALSE)
  ),
  Q97.5 = case_when(
    str_detect(Effect, "^b") ~ format(round(exp(`Q97.5`), 6), scientific = FALSE),
    TRUE ~ format(round(`Q97.5`,6) , scientific = FALSE)
  )
  )

# print the summary
print(post_summ_hie)

# posterior distribution
summary(model_bayes_hierarchical)

conditional_effects(model_bayes_hierarchical)

mcmc_plot(model_bayes_hierarchical, type = "trace")
stanplot(model_bayes_hierarchical, type = "dens_overlay")

pp_check(model_bayes_hierarchical)


