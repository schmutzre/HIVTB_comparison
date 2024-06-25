#### Library -------------------------------------------------------------------
library(dplyr)
library(lubridate)
library(janitor)
library(stringr)
library(flextable)
library(data.table)
library(haven)
library(labelled)

source("utils/functions.R")

#### Cleaning the separate data files and preparing to then join them ----------

#### Data ----------------------------------------------------------------------

treat <- read_dta("data_raw/CH/med_product.dta")
treatment3 <- read_dta("data_raw/CH/med_drug_code.dta") %>% 
  mutate(drug = drug_code)
lab <- read_dta("data_raw/CH/lab.dta")
complete_ch <- read_dta("data_raw/CH/shcs_509_hivall.dta")
var_region <- read_dta("data_raw/CH/var_region.dta") %>% 
  mutate(region = as.numeric(region))
dis <- read_dta("data_raw/CH/dis.dta")
var_disease <- read_dta("data_raw/CH/var_disease.dta")
admin <- read_dta("data_raw/CH/admin.dta")
modif <- read_dta("data_raw/CH/modif_art.dta")

#' The followng df are the study populations. 
#' Once patients starting ART between 2010-2022 and once patients having a TB diagnosis 
#' between 2010-2022. These dataframes are mostly overlapping. 

#### Filter studypopulation (Union of HIV- and TB-dataset) ---------------------

#' Vorgehen bei region: für TB Population habe ich Geburtsland genommen, bei HIV Citizenship

filteredBOTH <- complete_ch %>%
  # Step 1: Filtering
  filter(
    between(haart_start_date, as.Date("2010-01-01"), as.Date("2022-12-31")) |
      between(date_tb, as.Date("2010-01-01"), as.Date("2022-12-31"))
  ) %>%
  # Step 2: Selecting columns
  dplyr::select(
    id, born, sex, risk, haart_start_date, art_start_date, art_start_cd4, virus_type, 
    disease_tb, tbd_outcome, type_tb_shcs, date_tb, disease_tbc, tbd_pat_birth, region, 
    case_incident_2m, exitdate, current_art, eligibility_art, 
    starts_with("tbd_drug_resist_"), starts_with("tbd_antim_resist_"), regdate) %>%
  mutate(art_start_date = case_when(!is.na(haart_start_date) ~ haart_start_date,
                                    TRUE ~ art_start_date)) %>% 
  dplyr::select(-haart_start_date) %>% 
  # Step 3: Creating new variables
  mutate(
    cohort = as.factor("CH"),
    age_at_art_start = year(art_start_date) - born,
    presenting_tb = as_factor(case_when(
      date_tb < art_start_date + 60 & date_tb > art_start_date - 60 ~ 1,
      TRUE ~ 0)),
    recent_tb = as_factor(case_when(
      date_tb < art_start_date - 60 & date_tb > art_start_date - 360 ~ 1,
      TRUE ~ 0))) %>% 
  rename(outcome_tb = tbd_outcome) %>% 
  # Step 4: Merging with var_region dataset
  mutate(region = as.numeric(region)) %>%
  left_join(var_region, by = "region") %>%
  dplyr::select(-region) %>%
  rename(region = var_desc) %>%
  # Step 5: Reclassifying region
  mutate(
    region = trimws(region),
    region = as.factor(case_when(
      region %in% c("Southern Europe", "Western Europe", "Northern Europe", "Eastern Europe", "Northern America") ~ "Europe/Northern America",
      region %in% c("Eastern Africa", "Middle Africa", "Southern Africa", "Western Africa", "Northern Africa") ~ "Africa",
      region %in% c("Eastern Asia", "South-Eastern Asia", "Southern Asia", "Western Asia", "Oceania", "Central America", "Latin America and the Caribbean", "South America", "South/Latin America") ~ "South/Latin America/Asia/Oceania",
      region %in% "Unknown/World" ~ NA,
      TRUE ~ region
    )))

#### Type of infection ---------------------------------------------------------

filteredBOTH <- as_factor(filteredBOTH, levels="labels") %>% 
  mutate(risk = as.factor(case_when(risk == "Unknown" ~ NA,
                          TRUE ~ risk)))

#### Birth country / nationality -----------------------------------------------

get_continent <- function(country) {
    
  africa_countries <- c("Algeria", "Angola", "Cameroon", "Côte dIvoire", "Ethiopia", "Eritrea", "Gambia", "Ghana",
                          "Guinea", "Kenya", "Kongo (Brazzaville)", "Kongo (Kinshasa)", "Mauritania", "Nigeria",
                          "Somalia", "South Africa", "Togo", "Zimbabwe", "Uganda", "Egypt", "Morocco")
   
  if (country %in% africa_countries) {
      return("Africa")
    } else if (country %in% c("Russia", "Austria", "Romania", "Italy", "Spain", "Portugal", "Switzerland", "Russian Federation", "United States", "Canada")) {
      return("Europe/Northern America")
    } else if (country %in% c("Argentina", "Chile", "Brasil", "Peru", "Bolivia", "South America", "Latin America and the Caribbean", "Thailand", "Indonesia", "Cambodia", "Vietnam")) {
      return("South/Latin America/Asia/Oceania")
    }
    else {
      return(NA)
    } 
  }
  
#' Note, the birth country is mostly available for TB patients, but not for non-TB patients

filteredBOTH_region <- filteredBOTH %>%
  mutate(region_born = sapply(tbd_pat_birth, get_continent),
          region_born = case_when(
          !is.na(region_born) ~ region_born,
          TRUE ~ region),
         case_incident_2m = case_when(case_incident_2m == "Incident TB" ~ 1,
                                      TRUE ~ 0)) %>% 
  dplyr::select(-region, -tbd_pat_birth) %>% 
  mutate(region = as.factor(region_born)) %>% 
  dplyr::select(-region_born)

tabyl(filteredBOTH_region$region)

#### CD4 and VL baseline -------------------------------------------------------

lab.filtered <- lab %>% 
    dplyr::select(id:cd4date, cd4, rna) %>% 
    filter(id %in% filteredBOTH$id)
  
baseline_cd4 <- filteredBOTH_region %>% 
    left_join(lab, by = "id") %>% 
    mutate(cd4_baseline = ifelse(labdate >= (art_start_date - 180) & labdate <= (art_start_date + 30) & cd4 <= 2000, cd4, NA)) %>% 
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
    left_join(baseline_cd4 %>% 
                dplyr::select(id, cd4_baseline, labdate), by = "id") %>% 
    left_join(baseline_rna %>% 
                dplyr::select(id, rna_baseline, labdate), by = "id") %>% 
    dplyr::select(c(id, all_of(colnames(filteredBOTH_region)), cd4_baseline, labdate.x, rna_baseline, labdate.y)) %>% 
    rename(labdate_cd4 = labdate.x, labdate_rna = labdate.y) %>% 
    rename_with(~ str_replace_all(., "\\..$", ""))
  
## add groups of baseline values

filteredBOTH.lab <- filteredBOTH.lab %>%
    mutate(
      rna_group = as.factor(case_when(
      rna_baseline >= 0 & rna_baseline <= 999 ~ "0-999",
      rna_baseline >= 1000 & rna_baseline <= 9999 ~ "1000-9999",
      rna_baseline >= 10000 ~ "10000+",
      TRUE ~ "NA"
    )),
    cd4_group = as.factor(case_when(
      cd4_baseline >= 0 & cd4_baseline <= 99 ~ "0-99",
      cd4_baseline >= 100 & cd4_baseline <= 349 ~ "100-349",
      cd4_baseline >= 350 ~ "350+",
      TRUE ~ "NA"
    )))
  
#### WHO stage -----------------------------------------------------------------
  
dis2 <- dis %>% 
    filter(id %in% filteredBOTH.lab$id) %>% 
    left_join(var_disease, by = "disease") %>% 
    dplyr::select(id:disease, disease_id, cdc_group) %>% 
    mutate(cdc_group = ifelse(cdc_group =="D", "C", cdc_group))  
  
who_stages <- filteredBOTH.lab %>% 
  left_join(dis2, by = "id") %>% 
    mutate(cdc_group = ifelse(cdc_group =="", NA, cdc_group)) %>%
    group_by(id) %>%
    mutate(cdc_group =  
             ifelse(newdate <= (art_start_date + 180) & newdate >= (art_start_date - 180), cdc_group, NA)) %>%
    arrange(abs(art_start_date - newdate)) %>%
    filter(!is.na(cdc_group)) %>% 
    distinct(id, .keep_all = TRUE) %>% 
    mutate(who_stage = as.factor(case_when(
      presenting_tb == 1 ~ 3,
      cdc_group == "A" ~ 1,
      (cdc_group == "B" & cd4_baseline >= 350) ~ 2,
      (cdc_group == "B" & cd4_baseline < 350) ~ 3,
      cdc_group == "C" ~ 4,
      TRUE ~ NA)),
    who_stage = as.factor(case_when(who_stage %in% c(1,2) ~ "1/2", 
                                    who_stage %in% c(3,4) ~ "3/4",
                                    TRUE ~ NA))) %>% 
    ungroup() 

filteredBOTH.who <- filteredBOTH.lab %>% 
    left_join(who_stages %>% dplyr::select(id, who_stage, cdc_group), by = "id") 
  
#### CD4 and VL @ TB-date ------------------------------------------------------

# Create tb_cd4 dataset
tb_cd4 <- filteredBOTH.who %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_cd4 = ifelse(labdate >= (date_tb - 60) & labdate <= (date_tb + 60), cd4, NA)) %>%
    filter(!is.na(tb_diag_cd4)) %>%
    arrange(id, abs(date_tb - labdate)) %>%
    group_by(id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup()
  
# Create tb_rna dataset
tb_rna <- filteredBOTH.who %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_rna = ifelse(labdate >= (date_tb - 60) & labdate <= (date_tb + 60), rna, NA)
    ) %>%
    filter(!is.na(tb_diag_rna)) %>%
    arrange(id, abs(date_tb - labdate)) %>%
    group_by(id) %>%
    distinct(id, .keep_all = TRUE) %>%
    ungroup()
  
# Join the two datasets together
tb_cd4_rna <- tb_cd4 %>%
    dplyr::select(id, tb_diag_cd4) %>%
    full_join(tb_rna %>% 
                dplyr::select(id, tb_diag_rna), by = "id")
  
filteredBOTH.tb <- filteredBOTH.who %>%
    left_join(tb_cd4_rna, by = "id")

filteredBOTH.person <- filteredBOTH.tb %>%
  left_join(admin %>% dplyr::select(id, last_fup_date, stop), by = "id") %>% 
  mutate(fup_time = case_when(
    !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
    !is.na(stop) ~ as.numeric(min(difftime(last_fup_date, art_start_date, units = "days"),
                                  difftime(as.Date("2022-12-31"), art_start_date, units = "days"))),
    TRUE ~ as.numeric(difftime(as.Date("2022-12-31"), art_start_date, units = "days"))),
    ltfu = ifelse(is.na(exitdate) & !is.na(stop), 1, 0),
    death = ifelse(!is.na(exitdate),1,0))

#### Add ART-treatment ---------------------------------------------------------

modif2 <- modif %>% 
left_join(filteredBOTH %>% dplyr::select(id, art_start_date), by = "id") %>% 
  filter(treatment != "",
         id %in% filteredBOTH$id)%>% 
  mutate(time_diff_ART = moddate - art_start_date,
         time_diff_STOP = enddate - moddate) %>% 
  group_by(id) %>% 
  arrange(time_diff_ART) %>% 
  slice_min(time_diff_ART, n = 1)

filteredBOTH.regimen <- filteredBOTH.person %>% 
  left_join(modif2, by = "id") %>% 
  mutate(treatment = str_replace_all(treatment, "ETC", "FTC")) %>%
  mutate(
    treatment = case_when(
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "EFV") | str_detect(treatment, "NVP")) ~ "TDF + 3TC/FTC + EFV/NVP",
      str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "AZT + 3TC/FTC + DTG",
      str_detect(treatment, "TDF") & str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "TDF + AZT + 3TC/FTC + DTG",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "TDF + 3TC/FTC + DTG",
      str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "LPV") ~ "AZT + 3TC/FTC + LPV",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "LPV") ~ "TDF + 3TC/FTC + LPV",
      (str_detect(treatment, "AZT") | str_detect(treatment, "TDF")) & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "LPV") | str_detect(treatment, "ATV") | str_detect(treatment, "DTG")) ~ "AZT/TDF + 3TC/FTC + LPV/r or ATV/r or DTG",
      TRUE ~ str_replace_all(treatment, " ", " + ")  # Replace spaces with plus signs
),
    regimen = case_when(
      num_inti != 0 ~ "INSTI-based",
      num_pi != 0 ~ "PI-based",
      num_nnrti != 0 ~ "NNRTI-based",
      TRUE ~ "Other")) %>% 
  dplyr::select(-(num_art:art_start_date.y)) %>% 
  rename(art_start_date = art_start_date.x)

#17 are NNRTI and INSTI based

#### TB resistance -------------------------------------------------------------

df_tb <- complete_ch %>% 
  filter(disease_tb == 1)
var_labels <- labelled::var_label(df_tb)

cols_to_rename_res <- names(df_tb)[which(names(df_tb) %in% paste0("tbd_drug_resist___", 1:99))]
new_names_res <- unlist(var_labels[cols_to_rename_res])
names_vector_res <- setNames(cols_to_rename_res, new_names_res)

df_tb <- df_tb %>%
  dplyr::rename(!!!names_vector_res)

tb_resi <- df_tb %>%
  dplyr::select(id, date_tb, Rifampicin_R:Others_R) %>%
  pivot_longer(cols = Rifampicin_R:Others_R, names_to = "Drug", values_to = "value") %>% 
  filter(value == 1) %>%
  group_by(id) %>%
  summarise(
    resistance_tb = paste(gsub("_R$", "", Drug), collapse = ", ")) %>%
  ungroup()

filteredBOTH.resist <- filteredBOTH.regimen %>% 
  left_join(tb_resi , by = "id") %>% 
  dplyr::select(-c(tbd_drug_resist___1:tbd_drug_resist_others)) %>% 
  mutate(resistance_tb_any = case_when(!is.na(resistance_tb) ~ 1,
                   TRUE ~ 0),
         resistance_tb_mdr = as.factor(case_when(
           (str_detect(resistance_tb, "Rifampicin")| str_detect(resistance_tb, "Rifabutin")) & str_detect(resistance_tb, "Isoniazid") ~ 1,
           TRUE ~ 0
         )))

#### TB treatment --------------------------------------------------------------

cols_to_rename <- names(df_tb)[which(names(df_tb) %in% paste0("tbd_antim_resist___", 1:99))]
new_names <- unlist(var_labels[cols_to_rename])
names_vector <- setNames(cols_to_rename, new_names)

df_tb <- df_tb %>%
  dplyr::rename(!!!names_vector)

tb <- df_tb %>%
  dplyr::select(id, date_tb, RHZE_TX:Others_TX) %>% 
  pivot_longer(cols = RHZE_TX:Others_TX, names_to = "Drug", values_to = "value") %>%
  filter(value == 1) %>%
  group_by(id) %>%
  summarise(
    regimen_tb = as.factor(ifelse(any(Drug == "RHZE_TX"), "Standard", paste(gsub("_TX$", "", Drug), collapse = ", "))
  )) %>%
  ungroup()

filteredBOTH.tbtreatment <- filteredBOTH.resist %>% 
  left_join(tb , by = "id") %>% 
  dplyr::select(-c(tbd_antim_resist___1:tbd_antim_resist_others))

tableTB <- tabyl(filteredBOTH.tbtreatment$regimen_tb, show_na = FALSE)

filteredBOTH.tbtreatment <- filteredBOTH.tbtreatment %>% 
  mutate(regimen_tb_group = 
           as.factor(case_when(regimen_tb %in% c("Standard", "Rifampicin, Pyrazinamide, Isoniazid, Ethambutol", "Rifampicin", "RH, Pyrazinamide, Ethambutol", "Rifampicin, Pyrazinamide, Isoniazid, Ethambutol") ~ "Standard",
                                        regimen_tb %in% c("Pyrazinamide, Isoniazid, Ethambutol, Amikacin, Moxifloxacin, Cycloserine, PAS, Linezolid, Imipenem",
                                                          "RH, Pyrazinamide, Ethambutol, Amikacin, Moxifloxacin",
                                                          "Rifampicin, Pyrazinamide, Isoniazid, Ethambutol, Moxifloxacin, Linezolid",
                                                          "Rifampicin, Pyrazinamide, Isoniazid, Ethambutol, Moxifloxacin",
                                                          "Rifampicin, Pyrazinamide, Isoniazid, Ethambutol, Streptomycin, Moxifloxacin",
                                                          "Rifampicin, Isoniazid, Ethambutol, Levofloxacin",
                                                          "Rifampicin, Isoniazid, Ethambutol, Moxifloxacin",
                                                          "Rifampicin, Pyrazinamide, Isoniazid") ~ "HRZE (Rif, pyraz, ison, ethamb) plus at least one quinolone",
                               regimen == "Others" ~ "Others",
                               regimen == "Rifampicin, Pyrazinamide, Isoniazid" ~ "2HRZE 4HR",
                               str_detect(regimen_tb, "Rifabutin") ~ "Rifabutin based",
                                        TRUE ~ regimen_tb)))

tabyl(filteredBOTH.tbtreatment$regimen_tb_group)

#### TB outcome ----------------------------------------------------------------

filteredBOTH.tboutcome <- filteredBOTH.tbtreatment %>% 
  mutate(outcome_tb = case_when(outcome_tb %in% c("Treatment completed","Cured")  ~ "Sucessfull",
                                outcome_tb == "Treatment failed" ~ "Failed",
                                outcome_tb == "Died" ~ outcome_tb,
                                TRUE ~ NA
                             ))

#### Site of TB ----------------------------------------------------------------

filteredBOTH.tbsite <- filteredBOTH.tboutcome %>% 
  mutate(site_tb = as.factor(case_when(disease_tbc == "TBC" ~ "Pulmonary",
                             type_tb_shcs == "TEX" ~ "Extrapulmonary",
                             TRUE ~ NA)))

#### Last changes before saving ------------------------------------------------

final <- filteredBOTH.tbsite %>% 
  mutate(who_stage = as.factor(who_stage),
         id = as.character(id),
         born = born) %>% 
  filter(fup_time >=0) %>% 
  rename(resistant_tb = resistance_tb_any)

#### lab data -----------------------------------------------------------

lab_both <- lab.filtered %>% 
  arrange(id, labdate) %>% 
  left_join(filteredBOTH.tbsite %>% dplyr::select(id, art_start_date, disease_tb, date_tb, presenting_tb, sex, age_at_art_start, cd4_baseline, region, who_stage, regimen), by = "id") %>% 
  mutate(time_diff = as.numeric(labdate - art_start_date, units = "days")) %>% 
  filter(time_diff >= -60) 
  
lab_cd4 <- lab_both %>% 
  dplyr::select(-rna) %>% 
  filter(!is.na(cd4)) %>% 
  group_by(id) %>% 
  arrange(labdate) %>%
  mutate(timepoint = row_number()) %>% 
  ungroup() %>% 
  rename(date_cd4 = labdate) 


lab_rna <- lab_both %>% 
  dplyr::select(-cd4) %>% 
  filter(!is.na(rna)) %>% 
  group_by(id) %>% 
  arrange(labdate) %>%
  mutate(timepoint = row_number()) %>% 
  ungroup() %>% 
  rename(date_cd4 = labdate)

#### Selection study population ------------------------------------------------

## ART ##

art_ch <- final %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date,
         eligibility_art ==1,
         recent_tb == 0) %>% 
  mutate(incident_tb = as.factor(case_incident_2m)) %>% 
  dplyr::select(-virus_type, -type_tb_shcs, -disease_tbc, -risk, -art_start_cd4,-eligibility_art, -case_incident_2m, -moddate, -enddate, - time_diff_ART, -time_diff_STOP, -resistance_tb, -labdate_cd4, -labdate_rna, -current_art, -regdate, -cdc_group, -last_fup_date, -stop, -exitdate, -resistance_tb_mdr) %>% 
  rename(gender = sex)

saveRDS(art_ch, "data_clean/ch/art_ch.rds")

## ART (only not-presenting) ##

art_ch_noTB <- art_ch %>% 
  filter(presenting_tb == 0)

saveRDS(art_ch_noTB, "data_clean/ch/art_noTB_ch.rds")

## TB ##

tb_ch <- final %>% 
  filter(date_tb <= as.Date("2022-12-31") & date_tb >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) %>% 
  mutate(incident_tb = as.factor(case_incident_2m)) %>% 
  dplyr::select(-virus_type, -type_tb_shcs, -disease_tbc, -risk, -art_start_cd4, -eligibility_art, -case_incident_2m, -moddate, -enddate, -time_diff_ART, -time_diff_STOP, -resistance_tb, -labdate_cd4, -labdate_rna, -current_art, -regdate, -cdc_group)

saveRDS(tb_ch, "data_clean/ch/tb_ch.rds") 

## Lab data ##

lab_ch <- lab_both %>% 
  filter(id %in% art_ch$id)

saveRDS(lab_ch, "data_clean/ch/lab_ch.rds")

cd4_ch <- lab_cd4 %>% 
  filter(id %in% art_ch$id)

saveRDS(cd4_ch, "data_clean/ch/cd4_ch.rds")

rna_ch <- lab_rna %>% 
  filter(id %in% art_ch$id)

saveRDS(rna_ch, "data_clean/ch/rna_ch.rds")

### overview TB

#tb_overview <- tb_ch %>%
 # group_by(resistance_tb, regimen_tb_group) %>%
  #summarise(count = n()) %>%
  #mutate(proportion = round(count / sum(count),2))

