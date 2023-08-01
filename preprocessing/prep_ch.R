#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven
)

#### Data


lab2 <- read_dta("data_raw/lab.dta")
complete_ch <- read_dta("data_raw/shcs_509_hivall.dta")
var_region <- read_dta("data_raw/var_region.dta") %>% 
  mutate(region = as.numeric(region))

studypopulation_ch <- complete_ch %>% 
  filter(art_start_date <= as.Date("2022-12-31"),
         art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born)

#filter studypopulation (function)

process_data <- function(file_name, output_name) {
  # Paths to lab, dis and var_disease data are constant
  lab <- read_dta("data_raw/lab.dta")
  dis <- read_dta("data_raw/dis.dta")
  var_disease <- read_dta("data_raw/var_disease.dta")
  
  #add infection source
  
  numbers_risk <- c(0:7, 9)
  mapping_risk <- c("other sources", "homosexual contacts", "heterosexual contacts", "i.v.drug use", 
                    "i.v. drugs/sexual contacts", "clotting factors against hemophilia", "other blood products",
                    "perinatal transmission", "unknown/inconclusive")
  
  studypopulation_ch <- file_name %>% 
    select(id:regdate, caseyear_tb:newdate_tbc, disease_tex, newdate_tex, 
           tbd_year_diag:tbd_outcome, case_incident_1m:date_com, exitdate, 
           exit_why, riskgroup:current_art, cd4_dx_date:rna_dx, eligibility_art, risk, region) %>% 
    mutate(risk2 = factor(risk, levels = numbers_risk, labels = mapping_risk)) %>% 
    mutate(risk2 = ifelse(is.na(as.character(risk2)), "unknown/inconclusive", as.character(risk2))) %>% 
    select(- c(height, diagnosis_tbc: newdate_tbc, disease_pcp:date_com, rna_last_val:rna_last_date, haart_start_date))
  
  #add region and gender

  studypopulation_ch <- studypopulation_ch %>% 
    mutate(sex = factor(sex, levels = c(1, 2), labels = c("male", "female")),
           region = factor(region, levels = var_region$region, labels = var_region$var_desc))
  
  #add cd4 and viral load baseline
  
  lab <- lab %>% 
    select(id:cd4date, cd4, rna) %>% 
    filter(id %in% studypopulation_ch$id)
  
  baseline_cd4_ch <- studypopulation_ch %>% 
    left_join(lab, by = "id") %>% 
    mutate(cd4_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), cd4, NA)) %>% 
    filter(!is.na(cd4_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup() 
  
  baseline_rna_ch <- studypopulation_ch %>% 
    left_join(lab, by = "id") %>% 
    mutate(rna_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), rna, NA)) %>% 
    filter(!is.na(rna_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup()
  
  colnames_x <- paste0(setdiff(colnames(studypopulation_ch), "id"), ".x")
  
  studypopulation_ch_joined <- studypopulation_ch %>% 
    left_join(baseline_cd4_ch %>% select(id,cd4_baseline, labdate), by = "id") %>% 
    left_join(baseline_rna_ch %>% select(id,rna_baseline, labdate), by = "id") %>% 
    select(c(id, all_of(colnames(studypopulation_ch)), cd4_baseline, labdate.x, rna_baseline, labdate.y)) %>% 
    rename(labdate_cd4 = labdate.x, labdate_rna = labdate.y) %>% 
    rename_with(~ str_replace_all(., "\\..$", ""))
  
  # add groups of baseline values
  studypopulation_ch_joined <- studypopulation_ch_joined %>%
    mutate(rna_group = case_when(
      rna_baseline >= 0 & rna_baseline <= 999 ~ "0-999",
      rna_baseline >= 1000 & rna_baseline <= 9999 ~ "1000-9999",
      rna_baseline >= 10000 ~ "10000+",
      TRUE ~ "NA"
    ),
    cd4_group = case_when(
      cd4_baseline >= 0 & cd4_baseline <= 99 ~ "0-99",
      cd4_baseline >= 100 & cd4_baseline <= 349 ~ "100-349",
      cd4_baseline >= 350 ~ "350+",
      TRUE ~ "NA"
    )) 
  
  #add baseline WHO stage
  
  dis <- dis %>% 
    filter(id %in% studypopulation_ch$id) %>% 
    left_join(var_disease, by = "disease") %>% 
    select(id:disease, disease_id, cdc_group) %>% 
    mutate(cdc_group = ifelse(cdc_group =="D", "C", cdc_group))  
  
  who_baseline_ch <- studypopulation_ch_joined %>% 
    left_join(dis, by = "id") %>% 
    mutate(cdc_group = ifelse(cdc_group =="", NA, cdc_group)) %>%
    group_by(id) %>%
    mutate(cdc_group =  
             ifelse(newdate <= (art_start_date + 30), cdc_group, NA)) %>%
    arrange(abs(art_start_date - newdate)) %>%
    filter(!is.na(cdc_group)) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    mutate(who_stage = case_when(
      cdc_group == "A" ~ "1",
      cdc_group == "B" ~ "2/3",
      cdc_group == "C" ~ "4"
    )) %>% 
    ungroup() 
  
  studypopulation_ch_joined2 <- studypopulation_ch_joined %>% 
    left_join(who_baseline_ch %>% select(id, who_stage), by = "id") 
  
  #add cd4 and viral load at tb diagnosis
  
  # Create tb_diag_cd4 dataset
  tb_cd4_ch <- studypopulation_ch_joined2 %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_cd4 = ifelse(labdate >= (date_tb - 180) & labdate <= (date_tb + 15), cd4, NA)
    ) %>%
    filter(!is.na(tb_diag_cd4)) %>%
    arrange(id, abs(date_tb - labdate)) %>%
    group_by(id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup()
  
  # Create tb_diag_rna dataset
  tb_rna_ch <- studypopulation_ch_joined2 %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_rna = ifelse(labdate >= (date_tb - 180) & labdate <= (date_tb + 15), rna, NA)
    ) %>%
    filter(!is.na(tb_diag_rna)) %>%
    arrange(id, abs(date_tb - labdate)) %>%
    group_by(id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup()
  
  # Join the two datasets together
  tb_cd4_rna_ch <- tb_cd4_ch %>%
    select(id, tb_diag_cd4) %>%
    full_join(tb_rna_ch %>% select(id, tb_diag_rna), by = "id")
  
  studypopulation_ch_joined3 <- studypopulation_ch_joined2 %>%
    left_join(tb_cd4_rna_ch, by = "id")
    
# At the end, save the data
saveRDS(studypopulation_ch_joined3, file = paste0("data_clean/", output_name, ".rds"))
}

process_data(studypopulation_ch, "art_ch")

##### adding lab data ####

art_ch <- readRDS("data_clean/art_ch.rds")

lab2 <- lab2 %>% 
  select(id:cd4date, cd4, rna) %>% 
  filter(id %in% art_ch$id)

# First, sort the lab dataframe by id and labdate
lab2 <- lab2 %>% arrange(id, labdate)

# Add a column to count the occurrence of each id
lab2 <- lab2 %>% 
  group_by(id) %>% 
  mutate(row = row_number())

# Now reshape lab dataframe using pivot_wider
lab_wide <- lab2 %>%
  pivot_wider(names_from = row, 
              values_from = c(labdate, cd4, rna), 
              names_sep = "_") 

reshape_lab_data_ordered <- function(lab_data){
  lab_data <- lab_data %>% 
    arrange(id, labdate) %>%
    group_by(id) %>%
    mutate(row = row_number())
  
  # Split the data into separate dataframes for each measurement
  labdate_wide <- lab_data %>%
    select(id, row, labdate) %>%
    pivot_wider(names_from = row, 
                values_from = labdate, 
                names_prefix = "labdate_")
  
  cd4_wide <- lab_data %>%
    select(id, row, cd4) %>%
    pivot_wider(names_from = row, 
                values_from = cd4, 
                names_prefix = "cd4_")
  
  rna_wide <- lab_data %>%
    select(id, row, rna) %>%
    pivot_wider(names_from = row, 
                values_from = rna, 
                names_prefix = "rna_")
  
  # Combine the separate measurements back into one dataframe
  combined_wide <- reduce(list(labdate_wide, cd4_wide, rna_wide), left_join, by = "id")
  
  return(combined_wide)
}

# Use the function to reshape lab data
lab_wide_ordered <- reshape_lab_data_ordered(lab2)

# Now join with art_ch
art_ch_lab_ordered <- left_join(art_ch, lab_wide_ordered, by = "id") 

saveRDS(art_ch_lab_ordered, "data_clean/art_lab_ch.rds")
# Long format for plotting, adding time from art_start for each labdate

# Separate labdate, cd4, and rna columns
labdate_df <- lab_wide_ordered %>% 
  select(id, starts_with("labdate")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "labdate")

cd4_df <- lab_wide_ordered %>% 
  select(id, starts_with("cd4")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "cd4")

rna_df <- lab_wide_ordered %>% 
  select(id, starts_with("rna")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "rna")

cd4_df$timepoint <- gsub("cd4_", "", cd4_df$timepoint)
labdate_df$timepoint <- gsub("labdate_", "", labdate_df$timepoint)
rna_df$timepoint <- gsub("rna_", "", rna_df$timepoint)

# Join the labdate, cd4, and rna columns into one dataframe
lab_long <- labdate_df %>%
  full_join(cd4_df, by = c("id", "timepoint")) %>%
  full_join(rna_df, by = c("id", "timepoint"))

# Add a time_diff column

lab_long <- lab_long %>% 
  left_join(art_ch_lab_ordered, by = "id") %>% 
  select(id:sex, date_tb, art_start_date, cd4_group, rna_group, who_stage) %>% 
  mutate(time_diff = as.numeric(labdate - art_start_date, units = "days"))

saveRDS(lab_long, "data_clean/art_lab_ch_long.rds")

#### preparing incidence data ####

ch_inc <- readRDS("data_clean/art_ch.rds")  %>%
  mutate(art_start_date_2m = art_start_date + months(2), # Adding 2 months to the start date of antiretroviral therapy (ART) as incidence is defined as 2 months after art start, so the time in between shoouldn't be counted.
         persontime_years = ifelse(case_incident_2m == 1, # Computing the person-time in years. If the event of interest occurred, it's the time from ART start to the event. 
                                   as.numeric(difftime(date_tb, art_start_date_2m, units = "days")) / 360, #convert to years
                                   ifelse(!is.na(exitdate), # If the event did not occur, it's either the time from ART start to the exit date or to the end of the study.
                                          as.numeric(difftime(exitdate, art_start_date_2m, units = "days")) / 360,  
                                          as.numeric(difftime(as.Date("2022-12-31"), art_start_date_2m, units = "days")) /360
                                   )),
         age_art_start = year(art_start_date) - born) %>% 
  select(id, sex, date_tb, case_incident_2m, cd4_baseline, rna_baseline, cd4_group, rna_group, persontime_years, exitdate, art_start_date, art_start_date_2m, age_art_start) %>% 
  filter(persontime_years > 0) %>% # Excluding patients with less than 2 months follow-up (now only 5005 patients)
  mutate(cohort = as.factor("CH"), cases = sum(ch$case_incident_2m)) # Adding a cohort indicator

saveRDS(ch_inc, "data_clean/df_inc_ch.rds")  

