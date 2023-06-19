##### Libraries #####

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr,
  lubridate,
  rlang
)

##### data import #####

ch <- readRDS("data_clean/art_ch.rds")
#sa <- readRDS("data_clean/art_sa")

##### prepare dataframe for analysis #####

poisson_ch <- ch %>% 
  mutate(art_start_date_2m = art_start_date + months(2), 
         persontime_years = ifelse(case_incident_2m == 1,
                             as.numeric(difftime(date_tb, art_start_date_2m, units = "days")) / 360, #convert to years
                             ifelse(!is.na(exitdate),
                                    as.numeric(difftime(exitdate, art_start_date_2m, units = "days")) / 360,  
                                    as.numeric(difftime(as.Date("2022-12-31"), art_start_date_2m, units = "days")) /360  
                             ))) %>% 
  select(id, sex, date_tb, case_incident_2m, cd4_baseline, rna_baseline, cd4_group, rna_group, persontime_years, exitdate, art_start_date, art_start_date_2m) %>% 
  filter(persontime_years > 0) %>% #patients which exited or the study ended inside two months after ART.
  mutate(cohort = "CH")

poisson_sa <- #...

poisson <- #... #join them together 
  
##### analysis #####

### main analysis ----

model_main <- glm(case_incident_2m ~ cohort, offset=log(person_time_years), family="poisson", data=poisson)

# Print the results
print("Main Analysis")
print(summary(model_main))
print(exp(coef(model_main)))  # Incidence rate ratios

### subgroup analysis ----

# Run a separate model for each CD4 group
for(cd4_group in unique(poisson$cd4_group)) {
  # Subset the data
  df_subgroup <- poisson %>% filter(cd4_group == !!cd4_group)
  
  # Run the model
  model <- glm(case_incident_2m ~ cohort, offset=log(person_time_years), family="poisson", data=df_subgroup)
  
  # Print the results
  print(paste("CD4 Group:", cd4_group))
  print(summary(model))
  print(exp(coef(model)))  # Incidence rate ratios
}

# Run a separate model for each RNA group 
for(rna_group in unique(poisson$rna_group)) {
  # Subset the data
  df_subgroup <- poisson %>% filter(rna_group == !!rna_group) # include only the rows where the rna_group value is the same as the current value of rna_group in the loop, The !! before rna_group inside the filter function is a way to force R to use the current value of the rna_group variable in the loop. 
  
  # Run the model
  model <- glm(case_incident_2m ~ cohort, offset=log(person_time_years), family="poisson", data=df_subgroup)
  
  # Print the results
  print(paste("RNA Group:", rna_group))
  print(summary(model))
  print(exp(coef(model)))  # Incidence rate ratios
}

# Run a separate model for each combination of CD4 group and RNA group
for(cd4_group in unique(poisson$cd4_group)) {
  for(rna_group in unique(poisson$rna_group)) {
    # Subset the data
    df_subgroup <- poisson %>% filter(cd4_group == !!cd4_group, rna_group == !!rna_group)
    
    # Run the model
    model <- glm(case_incident_2m ~ cohort, offset=log(person_time_years), family="poisson", data=df_subgroup)
    
    # Print the results
    print(paste("CD4 Group:", cd4_group, "; RNA Group:", rna_group))
    print(summary(model))
    print(exp(coef(model)))  # Incidence rate ratios
  }
}
