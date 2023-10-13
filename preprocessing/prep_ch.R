#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven,
  janitor,
  gridExtra
)


checkTB <- read_dta("data_raw/shcs_509_tbcases.dta")

#### Data ----------------------------------------------------------------------

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

width_descr <- 11 / cm(1)
height_descr <- 8 / cm(1)

# The followng df are the study populations. Once patients starting ART between 2010-2022 and once patients having a TB diagnosis 
# between 2010-2022. These dataframes are mostly overlapping. 

#### Filter studypopulation (Union of HIV- and TB-dataset) ---------------------

filteredBOTH <- complete_ch %>%
  # Step 1: Filtering
  filter(
    between(art_start_date, as.Date("2010-01-01"), as.Date("2022-12-31")) |
      between(date_tb, as.Date("2010-01-01"), as.Date("2022-12-31"))
  ) %>%
  # Step 2: Selecting columns
  select(
    id, born, sex, risk, art_start_date, art_start_cd4, virus_type, 
    disease_tb, type_tb_shcs, date_tb, disease_tbc, tbd_pat_birth, region, 
    case_incident_2m, exitdate, exit_why, current_art, eligibility_art, 
    starts_with("tbd_drug_resist_"), starts_with("tbd_antim_resist_")
  ) %>%
  # Step 3: Creating new variables
  mutate(
    cohort = "CH",
    age_at_ART_start = year(art_start_date) - born,
    prevalent_TB = ifelse(date_tb < art_start_date + 60 & date_tb > art_start_date - 60, 1, 0),
    recent_TB = ifelse(date_tb < art_start_date - 60 & date_tb > art_start_date - 360, 1, 0),
    presenting_TB = ifelse(prevalent_TB == 1 | recent_TB == 1, 1, 0)
  ) %>%
  # Step 4: Merging with var_region dataset
  mutate(region = as.numeric(region)) %>%
  left_join(var_region, by = "region") %>%
  select(-region) %>%
  rename(region = var_desc) %>%
  # Step 5: Reclassifying region
  mutate(
    region = trimws(region),
    region = case_when(
      region %in% c("Southern Europe", "Western Europe", "Northern Europe", "Eastern Europe") ~ "Europe",
      region %in% c("Eastern Africa", "Middle Africa", "Southern Africa", "Western Africa") ~ "Sub-Saharan Africa",
      region %in% c("Eastern Asia", "South-Eastern Asia", "Southern Asia", "Western Asia") ~ "Asia",
      TRUE ~ region
    )
  )

#### Type of infection ---------------------------------------------------------

numbers_risk <- c(0:7, 9)
mapping_risk <- c(
  "other sources", "homosexual contacts", "heterosexual contacts", "i.v.drug use",
  "i.v. drugs/sexual contacts", "clotting factors against hemophilia", "other blood products",
  "perinatal transmission", "unknown/inconclusive"
)

filteredBOTH <- filteredBOTH %>% 
  mutate(
    risk2 = factor(risk, levels = numbers_risk, labels = mapping_risk),
    risk2 = replace_na(as.character(risk2), "unknown/inconclusive")
  )

  
#### Birth country / nationality -----------------------------------------------
  
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
      return("South/Latin America")
    } else if (country %in% c("Thailand", "Indonesia", "Cambodia", "Vietnam")) {
      return("Asia")
    }
    else {
      return("Unknown")
    } 
  }
  
# A second function if the birth country isnt available. 
map_to_continent <- function(region) {
  case_when(
    region %in% c("Western Africa", "Eastern Africa", "Middle Africa", "Southern Africa", "Sub-Saharan Africa") ~ "Sub-Saharan Africa",
    region == "Northern Africa" ~ "North Africa",
    region %in% c("Southern Europe", "Eastern Europe", "Northern Europe", "Western Europe", "Europe") ~ "Europe",
    region %in% c("Eastern Asia", "Southern Asia", "South-Eastern Asia", "Western Asia", "Central Asia", "Asia") ~ "Asia",
    region == "Northern America" ~ "North America",
    region %in% c("South America", "Latin America and the Caribbean", "Central America") ~ "South/Latin America",
    region == "Oceania" ~ "Oceania",
    TRUE ~ "Unknown"
  )
}

#' Note, the birth country is mostly available for TB patients, but not for non-TB patients

filteredBOTH_region <- filteredBOTH %>%
  mutate(region_born = case_when(disease_tb == 1 ~ sapply(tbd_pat_birth, get_continent),
                                 TRUE ~ sapply(tbd_pat_birth, get_continent)),
          region_born = case_when(
          region_born != "Unknown" ~ region_born,
          TRUE ~ map_chr(region, map_to_continent)),
         region_born = case_when(
           region_born != "Unknown" ~ region_born,
           TRUE ~ NA),
         case_incident_2m = case_when(case_incident_2m == "Incident TB" ~ 1,
                                      TRUE ~ 0),
         case_incident_1m = case_when(case_incident_2m == "Incident TB" ~ 1,
                                      TRUE ~ 0)) %>% 
  select(-region, -tbd_pat_birth) %>% 
  rename(region = region_born)

#### CD4 and VL baseline -------------------------------------------------------

lab.filtered <- lab %>% 
    select(id:cd4date, cd4, rna) %>% 
    filter(id %in% filteredBOTH$id)
  
baseline_cd4 <- filteredBOTH_region %>% 
    left_join(lab, by = "id") %>% 
    mutate(cd4_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), cd4, NA)) %>% 
    filter(!is.na(cd4_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup() 
  
baseline_rna <- filteredBOTH_region %>% 
    left_join(lab, by = "id") %>% 
    mutate(rna_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30), rna, NA)) %>% 
    filter(!is.na(rna_baseline)) %>% 
    arrange(abs(art_start_date - labdate)) %>% 
    group_by(id) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    ungroup()
  
colnames_x <- paste0(setdiff(colnames(filteredBOTH), "id"), ".x")
  
filteredBOTH.lab <- filteredBOTH_region %>% 
    left_join(baseline_cd4 %>% select(id,cd4_baseline, labdate), by = "id") %>% 
    left_join(baseline_rna %>% select(id,rna_baseline, labdate), by = "id") %>% 
    select(c(id, all_of(colnames(filteredBOTH_region)), cd4_baseline, labdate.x, rna_baseline, labdate.y)) %>% 
    rename(labdate_cd4 = labdate.x, labdate_rna = labdate.y) %>% 
    rename_with(~ str_replace_all(., "\\..$", ""))
  
## add groups of baseline values

filteredBOTH.lab <- filteredBOTH.lab %>%
    mutate(
      rna_group = case_when(
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
  
#### WHO stage -----------------------------------------------------------------
  
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
  
filteredBOTH.who <- filteredBOTH.lab %>% 
    left_join(who_stages %>% select(id, who_stage), by = "id") 
  
#### CD4 and VL @ TB-date ------------------------------------------------------

# Create tb_cd4 dataset
tb_cd4 <- filteredBOTH.who %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_cd4 = ifelse(labdate >= (date_tb - 180) & labdate <= (date_tb + 15), cd4, NA)
    ) %>%
    filter(!is.na(tb_diag_cd4)) %>%
    arrange(id, abs(date_tb - labdate)) %>%
    group_by(id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup()
  
# Create tb_rna dataset
tb_rna <- filteredBOTH.who %>%
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
    full_join(tb_rna %>% 
                select(id, tb_diag_rna), by = "id")
  
filteredBOTH.tb <- filteredBOTH.who %>%
    left_join(tb_cd4_rna, by = "id")

filteredBOTH.person <- filteredBOTH.tb %>%
  left_join(admin %>% select(id, last_fup_date), by = "id") %>% 
  mutate(last_persontime = case_when(
    !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
    TRUE ~ as.numeric(difftime(last_fup_date, art_start_date, units = "days"))
  ))
#' last persontime is defined as persontime until death or loss to follow up. For any incidence models (eg. TB or viral suppression)
#' the person years have to be re-calculated as I haven't incorporated persontime until incidence of any sort.

#### Add ART-treatment ---------------------------------------------------------

modif2 <- modif %>% 
left_join(filteredBOTH %>% select(id, art_start_date), by = "id") %>% 
  filter(treatment != "",
         id %in% filteredBOTH$id)%>% 
  mutate(time_diff_ART = moddate - art_start_date,
         time_diff_STOP = enddate - moddate) %>% 
  group_by(id) %>% 
  arrange(time_diff_ART) %>% 
  slice_min(time_diff_ART, n = 1)

filteredBOTH.regimen <- filteredBOTH.person %>% 
  left_join(modif2, by = "id") %>% 
  mutate(
    regimen = case_when(
      num_pi != 0 ~ "PI-based",
      num_nnrti != 0 & num_inti != 0 ~ "Other",
      num_nnrti != 0 ~ "NNRTI-based",
      num_inti != 0 ~ "INSTI-based",
      TRUE ~ "Other"
    )
  ) %>% 
  select(-(num_art:art_start_date.y))

filteredBOTH.regimen <- filteredBOTH.regimen %>% 
  mutate(current_art = case_when(current_art == "" ~ NA,
                                 TRUE ~ current_art))

# Generate frequency table
freq_table <- tabyl(filteredBOTH.regimen$current_art, show_na = FALSE)
freq_table2 <- tabyl(filteredBOTH.regimen$treatment, show_na = FALSE)
# Sort by proportion and get top 5
top_5 <- freq_table[order(-freq_table$percent), ][1:5, ]
top_52 <- freq_table[order(-freq_table2$percent), ][1:5, ]


#17 are NNRTI and INSTI based
# Note: there are some id's for which moddate = enddate. Will have to check with Lukas what that means.

#### TB resistance -------------------------------------------------------------

library(labelled)
df_tb <- complete_ch %>% 
  filter(disease_tb == 1)
var_labels <- labelled::var_label(df_tb)

cols_to_rename_res <- names(df_tb)[which(names(df_tb) %in% paste0("tbd_drug_resist___", 1:99))]
new_names_res <- unlist(var_labels[cols_to_rename_res])
names_vector_res <- setNames(cols_to_rename_res, new_names_res)

df_tb <- df_tb %>%
  dplyr::rename(!!!names_vector_res)

tb_resi <- df_tb %>%
  select(id, date_tb, Rifampicin_R:Others_R) %>%
  pivot_longer(cols = Rifampicin_R:Others_R, names_to = "Drug", values_to = "value") %>% 
  filter(value == 1) %>%
  group_by(id) %>%
  summarise(
    TB_resistance = paste(gsub("_R$", "", Drug), collapse = ", ")
  ) %>%
  ungroup()

filteredBOTH.resist <- filteredBOTH.regimen %>% 
  left_join(tb_resi , by = "id") %>% 
  select(-c(tbd_drug_resist___1:tbd_drug_resist_others)) %>% 
  mutate(TB_resistance = case_when(!is.na(TB_resistance) ~ 1,
                   TRUE ~ 0))

#### TB treatment --------------------------------------------------------------
test <- complete_ch %>% 
  select(tbd_antim_resist___1:tbd_antim_resist_others)
cols_to_rename <- names(df_tb)[which(names(df_tb) %in% paste0("tbd_antim_resist___", 1:99))]
new_names <- unlist(var_labels[cols_to_rename])
names_vector <- setNames(cols_to_rename, new_names)

df_tb <- df_tb %>%
  dplyr::rename(!!!names_vector)

tb <- df_tb %>%
  select(id, date_tb, RHZE_TX:Others_TX) %>% 
  pivot_longer(cols = RHZE_TX:Others_TX, names_to = "Drug", values_to = "value") %>%
  filter(value == 1) %>%
  group_by(id) %>%
  summarise(
    TB_regimen = ifelse(any(Drug == "RHZE_TX"), "standard", paste(gsub("_TX$", "", Drug), collapse = ", "))
  ) %>%
  ungroup()

filteredBOTH.tbtreatment <- filteredBOTH.resist %>% 
  left_join(tb , by = "id") %>% 
  select(-c(tbd_antim_resist___1:tbd_antim_resist_others))

tableTB <- tabyl(filteredBOTH.tbtreatment$TB_regimen, show_na = FALSE)

filteredBOTH.tbtreatment <- filteredBOTH.tbtreatment %>% 
  mutate(TB_regimen_group = 
           case_when(TB_regimen %in% c(tableTB[20,1], tableTB[13,1],tableTB[12,1],tableTB[10,1],tableTB[9,1],tableTB[7,1]) ~ "standard",
                                        TB_regimen %in% c(tableTB[19,1],tableTB[15,1],tableTB[14,1],tableTB[11,1],tableTB[10,1],tableTB[8,1],tableTB[4,1]) ~ "HRZE (Rif, pyraz, ison, ethamb) plus at least one quinolone",
                                        TB_regimen %in% c(tableTB[18,1],tableTB[17,1],tableTB[16,1],tableTB[6,1],tableTB[5,1],tableTB[4,1],tableTB[3,1],tableTB[1,1]) ~ "Rifabutin-based regimen, plus at least HZE +/- another drug",
                                        TRUE ~ NA))

##### adding lab data ----------------------------------------------------------

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
  mutate(time_diff = as.numeric(labdate - art_start_date.x, units = "days"))

#### Selecting study population and saving the clean files ---------------------

filteredBOTH.labwidefinal <- filteredBOTH.labwide %>% 
  select(-labdate_cd4, -labdate_rna, -current_art) %>% 
  rename(art_start_date = art_start_date.x)

filteredBOTH.lablongfinal <- filteredBOTH.lablong %>% 
  select(-labdate_cd4, -labdate_rna, -current_art) %>% 
  rename(art_start_date = art_start_date.x)

#wide

art_ch <- filteredBOTH.labwidefinal %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date,
         eligibility_art ==1) 

saveRDS(art_ch, "data_clean/art_ch_lab.rds")

art_ch_noLab <- art_ch %>% 
  select(id:TB_regimen)

saveRDS(art_ch_noLab, "data_clean/art_ch.rds")

tb_ch <- filteredBOTH.labwidefinal %>% 
  filter(date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 
  
saveRDS(tb_ch, "data_clean/tb_ch_lab.rds") 

tb_ch_noLab <- tb_ch %>% 
  select(id:TB_regimen)

saveRDS(tb_ch_noLab, "data_clean/tb_ch.rds") 

#long
art_ch.long <- filteredBOTH.lablongfinal %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 

saveRDS(art_ch.long, "data_clean/art_ch_lablong.rds")

tb_ch.long <- filteredBOTH.lablongfinal %>% 
  filter(date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) 

saveRDS(tb_ch.long, "data_clean/tb_ch_lablong.rds") 

### overview TB

tb_overview <- tb_ch %>%
  group_by(TB_resistance, TB_regimen_group) %>%
  summarise(count = n()) %>%
  mutate(proportion = round(count / sum(count),2))
