#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven,
  here
)

#### Data ####

shcs_509_hivall <- read_dta(here("data_raw/shcs_509_hivall.dta"))
lab <- read_dta(here("data_raw/lab.dta"))
dis <- read_dta(here("data_raw/dis.dta"))
var_disease <- read_dta(here("data_raw/var_disease.dta"))


#### Preprocessing ####

# shcs_509_hivall is already a good basis for the study population, I'll remove some columns first. 

shcs_509_hivall <- shcs_509_hivall %>% 
  select(id:regdate, caseyear_tb:newdate_tbc, disease_tex, newdate_tex, 
         tbd_year_diag:tbd_outcome, case_incident_1m:date_com, exitdate, 
         exit_why, riskgroup:current_art, cd4_dx_date:rna_dx, eligibility_art)

## HIV patients

studypopulation_hiv <- shcs_509_hivall %>% 
  filter(regdate <=as.Date("2022-12-31"),
         year(art_start_date) - born >= 16,
         art_start_date <=as.Date("2022-12-31"),
         art_start_date >=as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         is.na(exitdate) | exitdate >= art_start_date) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born)

saveRDS(studypopulation_hiv, here("data_clean/studypopulation_hiv.rds"))

# the number of patients which can be included (5044) varies from the number of patients which are eligibility_ART = 1 (5058) according to the dataset I received. I will now try to identify them and find out what the reasons are.

eligibility_art <- shcs_509_hivall %>%
  filter(eligibility_art ==1) #those are all the patients which are art eligible according to the dataset.

ids_only_in_df1 <- setdiff(eligibility_art$id, studypopulation_hiv$id) 
ids_only_in_df1 <- unlist(ids_only_in_df1) #those are all the patiens who are art eligible but not according to the procedure above.

check_ids <- shcs_509_hivall %>% # I will now apply one filter after another to check where they fall out of the filtering system.
  filter(id %in% ids_only_in_df1,
         regdate <=as.Date("2022-12-31"),#Two were registred in 2023
         year(art_start_date) - born >= 16) #The rest was removed here, this means I will keep 5,044 as the number of patients included in the study.

#The study includes a total of 5,044 patients. This number is different from the previous count of ART eligible patients (5,058) because it excludes 
#patients from the year 2023 and those who are under the age of 16 years.

## HIV/TB patients 

studypopulation_tb <-  shcs_509_hivall %>% 
  filter(disease_tb == 1,
         regdate <= as.Date("2022-12-31"),
         year(art_start_date) - born >= 16, #two patients excluded
         art_start_date <=as.Date("2022-12-31"),
         art_start_date >=as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         is.na(exitdate) | exitdate >= art_start_date,
         date_tb >= art_start_date) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born)

saveRDS(studypopulation_tb, here("data_clean/studypopulation_tb.rds"))

# The study includes a total of 75 patients which were HIV/TB co-infected.

## HIV without TB patients

studypopulation_without_tb <- studypopulation_hiv %>% 
  filter(disease_tb == 0)

saveRDS(studypopulation_without_tb, here("data_clean/studypopulation_without_tb.rds"))

## Laberatory data

laboratory <- lab %>% 
  select(id:cd4date, cd4, rna, weight, rna_limit, rna_method) %>% 
  filter(id %in% studypopulation_hiv$id)

saveRDS(laboratory, "data_clean/laboratory.rds")

## disease data

disease <- dis %>% 
  filter(id %in% studypopulation_hiv$id) %>% 
  left_join(var_disease, by = "disease") %>% 
  select(id:disease, disease_id, cdc_group) %>% 
  mutate(cdc_group = ifelse(cdc_group =="D", "C", cdc_group))  
  
saveRDS(disease, "data_clean/disease.rds")

## cd4 and rna mesaurements for each id

# Arrange the lab data by id and labdate
lab <- lab %>% arrange(id, labdate)

# Function to reshape lab data for dates and values
reshape_lab_data <- function(lab_data, column){
  lab_data %>%
    filter(!is.na(!!sym(column))) %>% #!! is used to unquote the column variable which is a input in the function
    arrange(id, labdate) %>%
    group_by(id, labdate) %>%
    summarise(!!column := first(!!sym(column)), .groups='drop') %>%
    group_by(id) %>%
    mutate(row = row_number()) %>%
    pivot_wider(names_from = row, 
                values_from = c(labdate, !!sym(column)), 
                names_prefix = paste0(column, "_"),
                names_sep = "")
}

# Apply the function to rna and cd4
lab_rna_pivot <- reshape_lab_data(lab, "rna") %>%
  select(id, sort(names(.)[-1]))

lab_cd4_pivot <- reshape_lab_data(lab, "cd4") %>%
  select(id, sort(names(.)[-1]))

# Create two separate dataframes for rna and cd4
hiv_rna <- studypopulation_hiv %>%
  distinct(id) %>%
  left_join(lab_rna_pivot, by = "id") 

hiv_cd4 <- studypopulation_hiv %>%
  distinct(id) %>%
  left_join(lab_cd4_pivot, by = "id") 

# Create vectors for labdate and rnarna columns
labdaterna_cols <- paste0("labdaterna_", 1:93)
rna_cols <- paste0("rnarna_", 1:93)

# Create vectors for labdate and rnarna columns
labdatecd4_cols <- paste0("labdatecd4_", 1:95)
cd4_cols <- paste0("cd4cd4_", 1:95)
  
# Interleave labdate and rnarna columns
column_order_rna <- c("id", unlist(purrr::transpose(list(labdaterna_cols, rna_cols))))

# Interleave labdate and cd4cd4 columns
column_order_cd4 <- c("id", unlist(purrr::transpose(list(labdatecd4_cols, cd4_cols))))

# Reorder the columns in dataframe
hiv_rna <- hiv_rna %>%
  select(all_of(column_order_rna)) %>% 
  left_join(studypopulation_hiv, by ="id") %>% 
  select(id:rnarna_93, date_tb, art_start_date) %>% 
  rename_with(~str_replace(., "rnarna", "rna"), starts_with("rnarna"))


# Reorder the columns in dataframe
hiv_cd4 <- hiv_cd4 %>%
  select(all_of(column_order_cd4)) %>% 
  left_join(studypopulation_hiv, by ="id") %>% 
  select(id:cd4cd4_95, date_tb, art_start_date) %>% 
  rename_with(~str_replace(., "cd4cd4", "cd4"), starts_with("cd4cd4"))

saveRDS(hiv_rna, "data_clean/lab_rna.rds") #Attention: labdates before art start aren't filterd out yet here
saveRDS(hiv_cd4, "data_clean/lab_cd4.rds")

# create long format (easier for plotting)
# Reshape the labdate columns
rna_labdate <- hiv_rna %>%
  pivot_longer(cols = starts_with("labdate"),
               names_to = "labdate",
               values_to = "date") %>%
  mutate(labdate_index = parse_number(labdate)) %>% 
  select(id, art_start_date, labdate, date, labdate_index)

cd4_labdate <- hiv_cd4 %>%
  pivot_longer(cols = starts_with("labdate"),
               names_to = "labdate",
               values_to = "date") %>%
  mutate(labdate_index = as.numeric(stringr::str_extract(labdate, "(?<=_)[0-9]+"))) %>% 
  select(id, art_start_date, labdate, date, labdate_index)

# Reshape the rnarna columns
rna_rna <- hiv_rna %>%
  pivot_longer(cols = starts_with("rna"),
               names_to = "rna",
               values_to = "rna_value") %>%
  mutate(rna_index = parse_number(rna)) %>% 
  select(id, art_start_date, rna, rna_value, rna_index)

cd4_cd4 <- hiv_cd4 %>%
  pivot_longer(cols = starts_with("cd4"),
               names_to = "cd4",
               values_to = "cd4_value") %>%
  mutate(cd4_index = as.numeric(stringr::str_extract(cd4, "(?<=_)[0-9]+"))) %>% 
  select(id, art_start_date, cd4, cd4_value, cd4_index)

# Join the two reshaped dataframes together
rna_long <- inner_join(rna_labdate, rna_rna, by = c("id", "labdate_index" = "rna_index")) %>% 
  arrange(id, labdate_index) %>% 
  group_by(id) %>% 
  mutate(time_diff = date - art_start_date.y) %>%
  filter(time_diff > 0) %>% 
  select(id, labdate, date, art_start_date.y, rna_value, time_diff) %>% 
  mutate(labdate = as.numeric(str_extract(labdate, "\\d+$")),
         time_diff = as.numeric(time_diff))

saveRDS(rna_long, "data_clean/lab_rna_long.rds")

cd4_long <- inner_join(cd4_labdate, cd4_cd4, by = c("id", "labdate_index" = "cd4_index")) %>% 
  arrange(id, labdate_index) %>% 
  group_by(id) %>% 
  mutate(time_diff = date - art_start_date.y) %>%
  filter(time_diff > 0) %>% 
  select(id, labdate, date, art_start_date.y, cd4_value, time_diff) %>% 
  mutate(labdate = as.numeric(str_extract(labdate, "\\d+$")),
         time_diff = as.numeric(time_diff))

lab_both_long <- merge(rna_long, cd4_long, by = c("id", "time_diff")) %>% 
  dplyr::select(id, date.x, art_start_date.y.x, rna_value, cd4_value, time_diff) %>% 
  arrange(id, time_diff)

saveRDS(lab_both_long, "data_clean/lab_both_long.rds") #labdates before art start are filtered out

lab_both_long_tb <- lab_both_long %>% 
  filter(id %in% studypopulation_tb$id) %>% 
  left_join(studypopulation_tb, by = "id") %>% 
  select(id, art_start_date.y.x, date.x, rna_value, cd4_value, time_diff, date_tb) %>% 
  rename(art_start_date = art_start_date.y.x, date = date.x) %>% 
  mutate(time_diff_tb = as.numeric(date_tb- art_start_date))

saveRDS(lab_both_long_tb, "data_clean/lab_both_long_tb.rds")
