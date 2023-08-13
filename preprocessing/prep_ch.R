#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven
)

#### Data ----
treatment3 <- read_dta("data_raw/med_drug_code.dta") %>% 
  mutate(drug = drug_code)
lab <- read_dta("data_raw/lab.dta")
complete_ch <- read_dta("data_raw/shcs_509_hivall.dta")
var_region <- read_dta("data_raw/var_region.dta") %>% 
  mutate(region = as.numeric(region))
dis <- read_dta("data_raw/dis.dta")
var_disease <- read_dta("data_raw/var_disease.dta")
admin <- read_dta("data_raw/admin.dta")
modif <- read_dta("data_raw/modif_art.dta")

# The followng df are the study populations. Once patients starting ART between 2010-2022 and once patients having a TB diagnosis 
# between 2010-2022. These dataframes are mostly overlapping. 

filteredBOTH <- complete_ch %>%
  filter(
    (art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01")) |
      (date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"))
  ) %>% 
  select(id, born, sex, regdate, risk, art_start_date, art_start_cd4, virus_type, pretreat, disease_tb, type_tb_shcs, date_tb, disease_tbc, tbd_pat_birth, region, case_incident_1m, case_incident_2m, exitdate, exit_why, current_art, eligibility_art, tbd_drug_resist___1:tbd_drug_resist_others) %>% 
  mutate(cohort = as.factor("CH")) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born,
         prevalent_TB = case_when(
           (date_tb < art_start_date+ 60 & date_tb > art_start_date - 60) ~ 1,
           TRUE ~ 0),
         recent_TB = case_when(
           (date_tb < art_start_date - 60 & date_tb > art_start_date - 360) ~ 1,
           TRUE ~ 0),
         presenting_TB = case_when(
           prevalent_TB == 1 | recent_TB ==1 ~ 1,
           TRUE ~ 0) 
         )
         
      
#### Prep ----

## Type of infection
numbers_risk <- c(0:7, 9)
mapping_risk <- c("other sources", "homosexual contacts", "heterosexual contacts", "i.v.drug use", 
                    "i.v. drugs/sexual contacts", "clotting factors against hemophilia", "other blood products",
                    "perinatal transmission", "unknown/inconclusive")
  
filteredBOTH <- filteredBOTH %>% 
    mutate(risk2 = factor(risk, levels = numbers_risk, labels = mapping_risk)) %>% 
    mutate(risk2 = ifelse(is.na(as.character(risk2)), "unknown/inconclusive", as.character(risk2)))
  
## add birth_country
  
filteredBOTH <- as_factor(filteredBOTH, levels="labels")
  
# Define a function to map countries to continents
get_continent <- function(country) {
    africa_countries <- c("Algeria", "Angola", "Cameroon", "Côte dIvoire", "Ethiopia", "Eritrea", "Gambia", "Ghana",
                          "Guinea", "Kenya", "Kongo (Brazzaville)", "Kongo (Kinshasa)", "Mauritania", "Nigeria",
                          "Somalia", "South Africa", "Togo", "Zimbabwe", "Uganda")
    
    if (country %in% africa_countries) {
      return("Sub-Saharan Africa")
    } else if (country %in% c("Egypt", "Morocco")) {
      return("North Africa")
    } else if (country %in% c("Russia", "Austria", "Romania", "Italy", "Spain", "Portugal", "Switzerland", "Russian Federation")) {
      return("Europe")
    } else if (country %in% c("United States", "Canada")) {
      return("North America")
    } else if (country %in% c("Argentina", "Chile", "Brasil", "Peru", "Bolivia")) {
      return("South America")
    } else if (country %in% c("Thailand", "Indonesia", "Cambodia", "Vietnam")) {
      return("Asia")
    }
    else {
      return("Unknown")
    } 
  }
  
# A second function for the imputation if the birth country isnt available. 
map_to_continent <- function(region) {
  case_when(
    region %in% c("Western Africa", "Eastern Africa", "Middle Africa", "Southern Africa") ~ "Sub-Saharan Africa",
    region == "Northern Africa" ~ "North Africa",
    region %in% c("Southern Europe", "Eastern Europe", "Northern Europe", "Western Europe") ~ "Europe",
    region %in% c("Eastern Asia", "Southern Asia", "South-Eastern Asia", "Western Asia", "Central Asia") ~ "Asia",
    region == "Northern America" ~ "North America",
    region %in% c("South America", "Latin America and the Caribbean") ~ "South America",
    region == "Oceania" ~ "Oceania",
    region == "Central America" ~ "Central America",
    TRUE ~ "Unknown"
  )
}

 ## Use mutate to create the 'region_born' column
filteredBOTH <- filteredBOTH %>%
  mutate(region_born = sapply(tbd_pat_birth, get_continent),
         region_born = ifelse(
           region_born != "Unknown",
           region_born,
           map_chr(region, map_to_continent)),
         region_born = ifelse(
           region_born != "Unknown",
           region_born,
           NA)) %>% 
  mutate(case_incident_2m = case_when(case_incident_2m == "Incident TB" ~ 1,
                                      TRUE ~ 0),
         case_incident_1m = case_when(case_incident_2m == "Incident TB" ~ 1,
                                      TRUE ~ 0))

## add cd4 and viral load baseline
  
lab.filtered <- lab %>% 
    select(id:cd4date, cd4, rna) %>% 
    filter(id %in% filteredBOTH$id)
  
baseline_cd4 <- filteredBOTH %>% 
    left_join(lab, by = "id") %>% 
    mutate(cd4_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), cd4, NA)) %>% 
    filter(!is.na(cd4_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup() 
  
baseline_rna <- filteredBOTH %>% 
    left_join(lab, by = "id") %>% 
    mutate(rna_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), rna, NA)) %>% 
    filter(!is.na(rna_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup()
  
colnames_x <- paste0(setdiff(colnames(filteredBOTH), "id"), ".x")
  
filteredBOTH.lab <- filteredBOTH %>% 
    left_join(baseline_cd4 %>% select(id,cd4_baseline, labdate), by = "id") %>% 
    left_join(baseline_rna %>% select(id,rna_baseline, labdate), by = "id") %>% 
    select(c(id, all_of(colnames(filteredBOTH)), cd4_baseline, labdate.x, rna_baseline, labdate.y)) %>% 
    rename(labdate_cd4 = labdate.x, labdate_rna = labdate.y) %>% 
    rename_with(~ str_replace_all(., "\\..$", ""))
  
## add groups of baseline values
filteredBOTH.lab <- filteredBOTH.lab %>%
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
  
## add baseline WHO stage
  
dis2 <- dis %>% 
    filter(id %in% filteredBOTH.lab$id) %>% 
    left_join(var_disease, by = "disease") %>% 
    select(id:disease, disease_id, cdc_group) %>% 
    mutate(cdc_group = ifelse(cdc_group =="D", "C", cdc_group))  
  
who_stages <- filteredBOTH.lab %>% 
    left_join(dis2, by = "id") %>% 
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
  
filteredBOTH.lab <- filteredBOTH.lab %>% 
    left_join(who_stages %>% select(id, who_stage), by = "id") 
  
## add cd4 and viral load at tb diagnosis
  
# Create tb_diag_cd4 dataset
tb_cd4 <- filteredBOTH.lab %>%
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
tb_rna <- filteredBOTH.lab %>%
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
tb_cd4_rna <- tb_cd4 %>%
    select(id, tb_diag_cd4) %>%
    full_join(tb_rna %>% select(id, tb_diag_rna), by = "id")
  
filteredBOTH.lab <- filteredBOTH.lab %>%
    left_join(tb_cd4_rna, by = "id")

filteredBOTH.lab2 <- filteredBOTH.lab %>%
  left_join(admin %>% select(id, last_fup_date), by = "id") %>% 
  mutate(last_persontime = case_when(
    !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
    TRUE ~ as.numeric(difftime(last_fup_date, art_start_date, units = "days"))
  ))
#last persontime is defined as persontime until death or loss to follow up. For any incidencde models (eg. TB or viral suppression)
#the person years have to be re-calculated as I havent incorporated persontime until incidence of any sort.

## add treatment at art start
modif2 <- modif %>% 
left_join(filteredBOTH %>% select(id, art_start_date), by = "id") %>% 
  filter(treatment != "",
         id %in% filteredBOTH$id)%>% 
  mutate(time_diff = moddate - art_start_date) %>% 
  group_by(id) %>% 
  arrange(time_diff) %>% 
  slice_min(time_diff, n = 1)

filteredBOTH.regimen <- filteredBOTH.lab2 %>% 
  left_join(modif2, by = "id") %>% 
  mutate(
    regimen = case_when(
      is.na(num_inti) & is.na(num_nnrti) & is.na(num_pi) ~ "unknown",
      num_pi != 0 ~ "PI-based",
      num_nnrti != 0 & num_inti != 0 ~ "Other",
      num_nnrti != 0 ~ "NNRTI-based",
      num_inti != 0 ~ "INSTI-based",
      TRUE ~ "Other"
    )
  ) %>% 
  select(-(num_art:art_start_date.y))

combinations_count <- filteredBOTH.lab2 %>%
  group_by(num_nrti, num_nnrti, num_pi, num_ntrti, num_fi, num_inti) %>%
  summarise(count = n())

# table of regimens
filteredBOTH.regimen %>%
  group_by(regimen) %>%
  summarise(count = n()) %>%
  mutate(percent = (count / sum(count)) * 100) 

# table of drug combinations
filteredBOTH.regimen %>% 
  group_by(treatment) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  mutate(percent = (count / sum(count)) * 100)

#17 are NNRTI and INSTI based
# Note: there are some id's for which moddate = enddate. Will have to check with Lukas what that means.

# Add TB resistance
filteredBOTH.resist <- filteredBOTH.regimen %>% 
  rowwise() %>%
  mutate(
    tb_drug_resist = case_when(
      # Check if any of the columns have a value of 1
      any(c_across(tbd_drug_resist___1:tbd_drug_resist_others) == "Checked") ~ 1,
      
      # Check if all columns have a value of 0
      all(c_across(tbd_drug_resist___1:tbd_drug_resist___99) == "Unchecked") ~ 0,
      
      # Return NA for all other cases
      TRUE ~ NA_integer_
    )
  ) %>%
  ungroup() %>%   # To remove the rowwise grouping
  select(-c(tbd_drug_resist___1:tbd_drug_resist_others))

# Add TB treatment

tb <- filteredBOTH.resist %>% 
  filter(disease_tb == 1)

drugs_tb <- read_dta("data_raw/med_drug.dta")  %>% 
  left_join(tb, by = "id") %>% 
  filter(disease_tb == 1) %>% 
  filter(startdate > date_tb - 5,
         startdate < date_tb + 60) %>% 
  mutate(tb_regimen = case_when(
    drug %in% c("INH", "RIF", "EMB", "PYZ") ~ "standard",
    TRUE ~ "other")) %>% 
  select(id, tb_regimen, date_tb, startdate) %>% 
  group_by(id) %>%
  slice(if (any(tb_regimen == "standard")) which(tb_regimen == "standard") else sample(1:n(), 1)) %>%
  ungroup() %>%
  # Ensuring no "standard" id appears in "other"
  distinct(id, .keep_all = TRUE)

filteredBOTH.tbtreatment <- filteredBOTH.resist %>% 
  left_join(drugs_tb %>% select(id, tb_regimen), by = "id")

##### adding lab data ----

# First, sort the lab dataframe by id and labdate
lab.filtered <- lab.filtered %>% 
  arrange(id, labdate)

# Add a column to count the occurrence of each id
lab.filtered <- lab.filtered %>% 
  group_by(id) %>% 
  mutate(row = row_number()) %>% 
  ungroup

# Now reshape lab dataframe using pivot_wider
lab_wide <- lab.filtered %>%
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
lab_wide_ordered <- reshape_lab_data_ordered(lab.filtered)

# Now join with filteredBOTH
filteredBOTH.labwide <- left_join(filteredBOTH.tbtreatment, lab_wide_ordered, by = "id") 

## Now i will also retransform in long format, as this is more conveniant for some tasks.

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

filteredBOTH.lablong <- lab_long %>% 
  left_join(filteredBOTH.labwide, by = "id") %>% 
  mutate(time_diff = as.numeric(labdate - art_start_date.x, units = "days")) %>% 
  select(id:tb_regimen)

#### Selecting study population and saving the clean files ----
#selecting the relevant columns for the last time
filteredBOTH.labwidefinal <- filteredBOTH.labwide %>% 
  select(-c(virus_type, pretreat, tbd_pat_birth, region, current_art, eligibility_art)) %>% 
  rename(art_start_date = art_start_date.x)

filteredBOTH.lablongfinal <- filteredBOTH.lablong %>% 
  select(-c(virus_type, pretreat, tbd_pat_birth, region, current_art, eligibility_art)) %>% 
  rename(art_start_date = art_start_date.x)

#wide

art_ch <- filteredBOTH.labwidefinal %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 

saveRDS(art_ch, "data_clean/art_ch.rds")

tb_ch <- filteredBOTH.labwidefinal %>% 
  filter(date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 
  
saveRDS(tb_ch, "data_clean/tb_ch.rds") 

#long
art_ch.long <- filteredBOTH.lablongfinal %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 

saveRDS(art_ch.long, "data_clean/art_ch.long.rds")

tb_ch.long <- filteredBOTH.lablongfinal %>% 
  filter(date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 

saveRDS(tb_ch.long, "data_clean/tb_ch.long.rds") 

  