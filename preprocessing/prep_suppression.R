library(tidyverse)

#### SHCS ####

## Suppression ##

ch <- readRDS("data_clean/art_ch_lab.rds")
ch.long <- readRDS("data_clean/art_ch_lablong.rds")
#sa.long <- 
long <- ch.long

#incidence of viral non-suppression
rna_suppression_treshold <- 50

# Sorting the data by id and date
lab_long <- long %>% 
  arrange(id, labdate) %>% 
  filter(labdate > art_start_date,
         !is.na(rna))

# Creating the next_rna_value and next_labdate columns
lab_long_next <- long %>% 
  group_by(id) %>% 
  filter(labdate > art_start_date) %>% 
  mutate(next_rna = lead(rna, order_by = labdate),
         next_labdate = lead(labdate, order_by = labdate)) # add the next_labdate column

# Print suppression_dates to see the results
suppression_dates <- lab_long_next %>%
  group_by(id) %>%
  filter(rna < rna_suppression_treshold & next_rna < rna_suppression_treshold) %>%
  slice_min(next_labdate, n = 1) %>% 
  dplyr::select(id, next_labdate)  # Select next_labdate instead of labdate

hiv_sup <- left_join(ch, suppression_dates, by = "id") %>% 
  rename(viral_suppression = next_labdate) %>% 
  dplyr::select(id, case_incident_2m, sex, cohort, born, art_start_date, cd4_group, rna_group, viral_suppression, last_persontime, exitdate, age_at_ART_start)  # Select next_labdate instead of labdate

viral_suppression_df <- hiv_sup %>%
  mutate(persontime_years.suppression = case_when(
    !is.na(viral_suppression)~ as.numeric(difftime(viral_suppression, art_start_date, units = "days")/360),
    is.na(viral_suppression) ~ last_persontime/360),
    last_date_supp = case_when( 
      !is.na(viral_suppression)~ viral_suppression,
      is.na(viral_suppression) ~ art_start_date + last_persontime)) %>% 
  mutate(event_type = case_when(
    !is.na(viral_suppression) ~ 1,
    !is.na(exitdate) ~ 2,
    TRUE ~0),
    incidence.sup = case_when(!is.na(viral_suppression) ~ 1,
                              TRUE ~ 0)) 
  
saveRDS(object = viral_suppression_df, file = "data_clean/supress_df.rds")

## Rebound ##

lab_long_rebound <- lab_long %>% 
  inner_join(suppression_dates, by = "id") %>%   # Join by id
  filter(labdate > next_labdate,
         !is.na(rna)) %>% 
  rename(suppression = next_labdate) %>% 
  group_by(id) %>% 
  mutate(rebound_rna = lead(rna, order_by = labdate),
         rebound = lead(labdate, order_by = labdate)) %>% 
  filter(rebound_rna > rna_suppression_treshold) %>% 
  arrange(id, labdate) %>% 
  slice_head(n=1)

new_rows <- viral_suppression_df %>% 
  filter(incidence.sup == 1) %>% 
  anti_join(lab_long_rebound, by = "id")

# Combining the dataframes
lab_long_rebound <- bind_rows(lab_long_rebound, new_rows) %>%
  mutate(supp_until_rebound_or_lasfup = case_when(!is.na(rebound) ~ as.numeric(difftime(rebound, suppression, units = "days")/360),
                                        TRUE ~ last_persontime /360),
         indicator_rebound = case_when(!is.na(rebound) ~ 1,
                   TRUE ~ 0))

saveRDS(object = lab_long_rebound, file = "data_clean/rebound.rds")
