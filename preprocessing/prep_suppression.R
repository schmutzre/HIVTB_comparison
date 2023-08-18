library(tidyverse)

#### SHCS ####
ch <- readRDS("data_clean/art_ch.rds")
ch.long <- readRDS("data_clean/art_ch.long.rds")
#sa.long <- 
long <- ch.long

#incidence of viral non-suppression
rna_supression_treshold <- 400
cd4_suppression_treshold <- 350

# Sorting the data by id and date
lab_long <- long %>% 
  arrange(id, labdate) %>% 
  filter(labdate > art_start_date)

# Creating the next_rna_value, next_cd4_value, and next_labdate columns
lab_long_next <- long %>% 
  group_by(id) %>% 
  filter(labdate > art_start_date) %>% 
  mutate(next_rna = lead(rna, order_by = labdate),
         next_cd4 = lead(cd4, order_by = labdate),
         next_labdate = lead(labdate, order_by = labdate)) # add the next_labdate column

# Filter the rows where both the current and next rna_value are >400 and cd4_value are <350
non_suppression_dates <- lab_long_next %>%
  group_by(id) %>%
  filter(rna > rna_supression_treshold & next_rna > rna_supression_treshold & cd4 < cd4_suppression_treshold & next_cd4 < cd4_suppression_treshold) %>%
  slice_min(next_labdate, n = 1) %>% # Use slice_min on next_labdate
  dplyr::select(id, next_labdate)  # Select next_labdate instead of labdate

# Print non_suppression_dates to see the results
suppression_dates <- lab_long_next %>%
  group_by(id) %>%
  filter(rna < rna_supression_treshold & next_rna < rna_supression_treshold & cd4 > cd4_suppression_treshold & next_cd4 > cd4_suppression_treshold) %>%
  slice_min(next_labdate, n = 1) %>% 
  dplyr::select(id, next_labdate)  # Select next_labdate instead of labdate

# Joining the suppression_dates with hiv
hiv_nonsup <- left_join(ch, non_suppression_dates, by = "id") %>% 
  rename(viral_non_suppression = next_labdate) %>% 
  dplyr::select(id, case_incident_2m, sex, cohort, born, art_start_date, cd4_group, rna_group, viral_non_suppression, last_persontime, exitdate, age_at_ART_start)  # Select next_labdate instead of labdate

hiv_sup <- left_join(ch, suppression_dates, by = "id") %>% 
  rename(viral_suppression = next_labdate) %>% 
  dplyr::select(id, case_incident_2m, sex, cohort, born, art_start_date, cd4_group, rna_group, viral_suppression, last_persontime, exitdate, age_at_ART_start)  # Select next_labdate instead of labdate

viral_non_suppression_df <- hiv_nonsup %>%
  mutate(persontime_years.NOsuppression = case_when(
    !is.na(viral_non_suppression)~ as.numeric(difftime(viral_non_suppression, art_start_date, units = "days")/360),
    is.na(viral_non_suppression) ~ last_persontime/360)) %>%
  mutate(event_type = case_when(
  !is.na(viral_non_suppression) ~ 1,
  !is.na(exitdate) ~ 2,# Loss to follow-up isnt considered a competing risk
  TRUE ~0),
  incidence.nosup = case_when(!is.na(viral_non_suppression) ~ 1,
                            TRUE ~ 0))

saveRDS(object = viral_non_suppression_df, file = "data_clean/NOsupress_df.rds")

viral_suppression_df <- hiv_sup %>%
  mutate(persontime_years.suppression = case_when(
    !is.na(viral_suppression)~ as.numeric(difftime(viral_suppression, art_start_date, units = "days")/360),
    is.na(viral_suppression) ~ last_persontime/360)) %>% 
  mutate(event_type = case_when(
    !is.na(viral_suppression) ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~0),
    incidence.sup = case_when(!is.na(viral_suppression) ~ 1,
                              TRUE ~ 0)) 
  
saveRDS(object = viral_suppression_df, file = "data_clean/supress_df.rds")
