#### Library -------------------------------------------------------------------

library(dplyr)
library(lubridate)
library(janitor)

#### Cleaning the separate data files and preparing to then join them ----------

## ART treatment ##
tblART <- read.csv("data_raw/RSA/tblART.csv") %>% 
  mutate(across(c(art_sd, art_ed), ~na_if(.x, "")),
         across(c(art_sd, art_ed), ~as.Date(.x, format = "%Y-%m-%d"))) %>% 
  select(patient, art_id, art_combination, art_sd, art_ed)

## Base table ##

tblBAS <- read.csv("data_raw/RSA/tblBAS_incl_naive_imp2.csv") %>% 
  mutate(across(c(birth_d, enrol_d, recart_d), ~na_if(.x, "")),
         across(c(birth_d, enrol_d, recart_d), ~as.Date(.x, format = "%Y-%m-%d")),
         age_at_art_start = interval(birth_d, recart_d) %/% years(1),
         sex = as.factor(case_when(sex == 1 ~ "Male",
                         sex == 2 ~ "Female",
                         TRUE ~ NA))) %>%  
  rename(born = birth_d, 
         art_start_date = recart_d) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  select(patient, enrol_d, born, sex, art_start_date, age_at_art_start, naive_y, naive_y_imp)

tblNAIVE <- tblBAS %>% 
  filter(naive_y_imp %in% c(1,NA)) %>% 
  select(patient, art_start_date, enrol_d) %>% 
  mutate(diff = enrol_d - art_start_date) %>% 
  filter(diff <= 30 | is.na(diff))

## CD4 Lab ##

tblLAB_CD4 <- read.csv("data_raw/RSA/tblLAB_CD4.csv") %>% 
  filter(cd4_u == 1) %>% 
  rename(cd4 = cd4_v) %>% 
  mutate(date_cd4 = as.Date(cd4_d, format = "%Y-%m-%d")) %>% 
  select(patient, date_cd4, cd4)

## RNA Lab ##

tblLAB_RNA <- read.csv("data_raw/RSA/tblLAB_RNA.csv") %>% 
  filter(!is.na(rna_v)) %>% 
  rowwise() %>% 
  mutate(rna = case_when(
    rna_v < 0 ~ round(runif(1, 0, abs(rna_v)), 0),
    TRUE ~ rna_v
  )) %>% 
  ungroup() %>% 
  mutate(date_rna = as.Date(rna_d, format = "%Y-%m-%d")) %>% 
  select(patient, date_rna, rna)

## Exit and last fup date ##

tblLTFU <- read.csv("data_raw/RSA/tblLTFU.csv") %>% 
  select(patient, drop_d, death_d) %>% 
  rename(last_fup_date = drop_d,
         exitdate = death_d) %>% 
  mutate(across(c(last_fup_date, exitdate), ~na_if(.x, "")),
         across(c(last_fup_date, exitdate), ~as.Date(.x, format = "%Y-%m-%d")))

last_visit <- read.csv("data_raw/RSA/tblVIS.csv") %>% 
  select(patient, vis_d) %>% 
  mutate(vis_d = as.Date(vis_d, format = "%Y-%m-%d")) %>% 
  group_by(patient) %>% 
  arrange(desc(vis_d)) %>% 
  slice(1) %>%
  ungroup()

tblLTFU2 <- tblLTFU %>% 
  left_join(last_visit, by = "patient") %>% 
  filter(vis_d > last_fup_date)

#' the variable last_fup_date is OK
#' i was checking if there where some visits after the last fup-date which would indicate
#' that the variable doesn't really capture the last fup-date.

## TB resistance ##

tbltb_res <- read.csv("data_raw/RSA/tbltb_res.csv") %>%
  group_by(patient) %>%
  mutate(
    resistance_tb_any = as.factor(case_when(
      any(drug_res == 1) ~ 1,
      TRUE ~ 0)),
    resistance_tb_rif = as.factor(case_when(
      any(drug_res == 1 & tb_drug == "RIF") ~ 1, 
      TRUE ~ 0)),
    resistance_tb_inh = as.factor(case_when(
      any(drug_res == 1 & grepl("^INH", tb_drug)) ~ 1, 
      TRUE ~ 0)),
    resistance_tb_mdr = as.factor(case_when(
      resistance_tb_rif == 1 & resistance_tb_inh == 1 ~ 1,
      TRUE ~ 0))) %>% 
  ungroup() %>%
  distinct(patient, .keep_all = TRUE) %>%
  select(patient, resistance_tb_any, resistance_tb_mdr)

## Med Table ##

#' read.csv("data_raw/RSA/")
#' Chido will send a tblMED which includes all meds given, we can use this to 
#' the TB treatment in tblTB and also the ART treatment maybe

## TB table ##

#' Chido will send us a new TB table

tblTB <- read.csv("data_raw/RSA/tblTB.csv") %>% 
  group_by(patient) %>%
  filter(reg_dmy == min(reg_dmy),
         startsWith(patient, "K"),
         cat %in% c(1,95,99,NA)) %>%
  ungroup() %>% 
  rename(date_tb = reg_dmy, 
         regimen_tb = regimen,
         site_tb = class) %>% 
  select(patient, date_tb, regimen_tb, site_tb, tb_outcome) %>% 
  mutate(across(c(date_tb, site_tb), ~na_if(.x, "")),
         date_tb = as.Date(date_tb, format = "%Y-%m-%d"),
         disease_tb = as.factor(1),
         site_tb = as.factor(case_when(site_tb %in% c(1,3) ~ "Pulmonary",
                             site_tb == 2 ~ "Extrapulmonary",
                             TRUE ~ NA)),
         regimen_tb = as.factor(case_when(regimen_tb == 1 ~ "2HRZE 4HR",
                                          regimen_tb == 2 ~ "2HRZES 1HRZE 5HRE",
                                          regimen_tb == 3 ~ "2HRZ 4HR",
                                          regimen_tb == 4 ~ "Others",
                                          TRUE ~ NA)),
         regimen_tb_group = as.factor(case_when(regimen_tb %in% c("2HRZE 4HR", "2HRZ 4HR") ~ "Standard",
                                                TRUE ~ regimen_tb)),
         outcome_tb = case_when(tb_outcome == 1 ~ "Completed",
                                tb_outcome == 2 ~ "Cured",
                                tb_outcome == 3 ~ "Failed",
                                tb_outcome == 4 ~ "Interrupted",
                                tb_outcome == 5 ~ "Defaulted",
                                tb_outcome == 6 ~ "Treatment ongoing",
                                tb_outcome == 7 ~ "Died",
                                TRUE ~ NA))


## WHO stage ##

tblVIS <- read.csv("data_raw/RSA/tblVIS.csv") %>% 
  filter(patient %in% tblNAIVE$patient) %>% 
  select(patient, vis_d, who_stage) %>% 
  mutate(vis_d = as.Date(vis_d, format = "%Y-%m-%d"),
         who_stage = as.factor(who_stage)) %>% 
  left_join(tblBAS %>% select(patient, art_start_date), by = "patient") %>% 
  group_by(patient) %>%
  mutate(difference = abs(vis_d - art_start_date),
         who_stage = case_when(who_stage %in% c(2,3) ~ "2/3",
                               who_stage == 9 ~ NA,
                               TRUE ~ who_stage)) %>%
  filter(difference <= 180) %>%
  arrange(patient, difference) %>% #per patient, then first the non-NA values, then arranged per time difference
  slice(1) %>%
  ungroup() %>%
  select(patient, who_stage)

#' hier muss ich anders slicen, die ersten nicht-NA beobachtungen nehmen

## CD4 baseline values ##

baseline_cd4 <- tblLAB_CD4 %>% 
  left_join(tblBAS %>% select(patient, art_start_date), 
            by = "patient") %>% 
  mutate(cd4_baseline = as.numeric(ifelse(date_cd4 >= (art_start_date - 180) & date_cd4 <= (art_start_date + 30), 
                               cd4, NA))) %>%
  filter(!is.na(cd4_baseline)) %>% 
  arrange(abs(art_start_date - date_cd4)) %>% 
  group_by(patient) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  ungroup() %>% 
  mutate(cd4_group = as.factor(
    case_when(
      cd4_baseline >= 0 & cd4_baseline <= 99 ~ "0-99",
      cd4_baseline >= 100 & cd4_baseline <= 349 ~ "100-349",
      cd4_baseline >= 350 ~ "350+",
      TRUE ~ NA
    )
  )) %>% 
  select(patient, cd4_baseline, cd4_group)

## RNA baseline values ##

baseline_rna <- tblLAB_RNA %>% 
    left_join(tblBAS %>% select(patient, art_start_date), 
              relationship = "many-to-many", 
              by = "patient") %>% 
    mutate(rna_baseline = as.numeric(ifelse(date_rna >= (art_start_date - 180) & date_rna <= (art_start_date + 30), 
                                 rna, NA))) %>% 
    filter(!is.na(rna_baseline)) %>% 
    arrange(abs(art_start_date - date_rna)) %>% 
    group_by(patient) %>% 
    distinct(patient, .keep_all = TRUE) %>% 
    ungroup() %>% 
    mutate(
        rna_group = as.factor(
          case_when(
          rna_baseline >= 0 & rna_baseline <= 999 ~ "0-999",
          rna_baseline >= 1000 & rna_baseline <= 9999 ~ "1000-9999",
          rna_baseline >= 10000 ~ "10000+",
          TRUE ~ NA
        ))) %>% 
    select(patient, rna_baseline, rna_group)

## CD4 @ TB diagnosis ##
  
tb_cd4 <- tblLAB_CD4 %>% 
    left_join(tblTB %>% select(patient, date_tb), 
              relationship = "many-to-many", 
              by = "patient") %>% 
    mutate(tb_diag_cd4 = ifelse(date_cd4 >= (date_tb - 120) & date_cd4 <= (date_tb + 120), 
                                 cd4, NA)) %>% 
    filter(!is.na(tb_diag_cd4)) %>% 
    arrange(abs(date_tb - date_cd4)) %>% 
    group_by(patient) %>% 
    distinct(patient, .keep_all = TRUE) %>% 
    ungroup() %>% 
    select(patient, tb_diag_cd4)
  
## RNA @ TB diagnosis ## 
  
tb_rna <- tblLAB_RNA %>% 
    left_join(tblTB %>% select(patient, date_tb), 
              relationship = "many-to-many", 
              by = "patient") %>% 
    mutate(tb_diag_rna = as.numeric(ifelse(date_rna >= (date_tb - 120) & date_rna <= (date_tb + 120), 
                                rna, NA))) %>% 
    filter(!is.na(tb_diag_rna)) %>% 
    arrange(abs(date_tb - date_rna)) %>% 
    group_by(patient) %>% 
    distinct(patient, .keep_all = TRUE) %>% 
    ungroup() %>% 
    select(patient, tb_diag_rna)
  
## ART treatment ##
  
id_to_drug <- tibble(
  art_combination = c("J05AG03", "J05AF07", "J05AF09", "J05AF05", "J05AF04", "J05AG01", "J05AR10", "J05AF06", "J05AJ03", "J05AF01",
                      "J05AF02", "J05A", "J05AE03", "J05AE08", "J05AJ01", "J05AG02", "J05AG04", "J05AE10", "J05AX07", "J05AE02",
                      "J05AF03", "J05AE04", "J05AR11", "J05AR06", "J05AR07", "J05AR02", "J05AR27", "J05AR01", "J05AR05", "J05AR03",
                      "J05AR13", "J05AR12", "J05AR23", "J05AR04"),
  art_id = c(
    "J05AG03", "J05AF07", "J05AF09", "J05AF05", "J05AF04", "J05AG01", "J05AR10", "J05AF06", "J05AJ03", "J05AF01",
    "J05AF02", "J05A", "J05AE03", "J05AE08", "J05AJ01", "J05AG02", "J05AG04", "J05AE10", "J05AX07", "J05AE02",
    "J05AF03", "J05AE04", "J05AR11", "J05AR06", "J05AR07", "J05AR02", "J05AR27", "J05AR01", "J05AR05", "J05AR03",
    "J05AR13", "J05AR12", "J05AR23", "J05AR04"),
  drug = c(
    "EFV", "TDF", "FTC", "3TC", "d4T", "NVP", "LPV/r", "ABC", "DTG", "AZT", "ddI", "ART unspecified", "/r",
    "ATV", "RAL", "DLV", "ETV", "DRV", "T20", "IDV", "ddC", "NFV", "3TC, TDF, EFV", "FTC, TDF, EFV",
    "d4T, 3TC, NVP", "3TC, ABC", "3TC, TDF, DTG", "AZT, 3TC", "AZT, 3TC, NVP", "TDF, FTC", "3TC, ABC, DTG", "3TC, TDF",
    "ATV/r", "AZT, 3TC, ABC"))

treatment_art <- tblART %>%
  left_join(tblBAS %>% select(patient, art_start_date), 
            relationship = "many-to-many",
            by = "patient") %>% 
  filter(abs(art_sd - art_start_date) <= 120) %>% 
  left_join(id_to_drug %>% select(drug, art_id), by = c("art_id")) %>% 
  left_join(id_to_drug %>% select(drug, art_combination), by = c("art_combination")) %>% 
  group_by(patient) %>%
  mutate(abs_diff = abs(art_sd - art_start_date)) %>%
  arrange(abs_diff, .by_group = TRUE) %>%
  summarize(
    treatment = coalesce(
      first(na.omit(drug.y)),
      paste(unique(na.omit(drug.x)), collapse = ", "),
      NA
    )
  )

treatment_art <- treatment_art %>% 
  mutate(
    treatment = case_when(
      treatment == "3TC, TDF, DTG" ~ "TDF + 3TC/FTC + DTG",
      treatment %in% c("3TC, TDF, EFV", 
                       "FTC, TDF, EFV",
                       "TDF, FTC, EFV",
                       "TDF, EFV, 3TC",
                       "NVP, 3TC, TDF",
                       "EFV, 3TC, TDF",
                       "3TC, NVP, TDF",
                       "3TC, TDF, FTC, EFV",
                       "3TC, TDF, TFC, NVP, EFV",
                       "3TC, TDF, NVP",
                       "3TC, TDF, NVP, EFV",
                       "ddl, 3TC, TDF, EFV",
                       "d4T, 3TC, TDF, FTC, EFV",
                       "d4T, 3TC, TDF, EFV, LPV/r",
                       "d4T, 3TC, TDF, EFV",
                       "3TC, TDF, FTC, NVP, EFV",
                       "3TC, TDF, EFV, DTG",
                       "AZT, 3TC, TDF, NVP",
                       "AZT, TDF, NVP, 3TC") ~ "TDF + 3TC/FTC + EFV/NVP",
      treatment == "AZT, 3TC, EFV, LPV/r" ~ "AZT + 3TC/FTC + LPV/r",
      TRUE ~ treatment,
    ),
    regimen = case_when(treatment == "TDF + 3TC/FTC + EFV/NVP" ~ "NNRTI-based",
                        treatment %in% c("TDF + 3TC/FTC + DTG", "AZT + 3TC/FTC + LPV/r") ~ "INSTI-based",
                        TRUE ~ NA)
  )

#### Joining them together -----------------------------------------------------
  
df <- tblBAS %>% 
    left_join(tblTB, by = "patient") %>% 
    left_join(tblLTFU, by = "patient") %>% 
    left_join(baseline_cd4, by = "patient") %>% 
    left_join(baseline_rna, by = "patient") %>% 
    left_join(tb_cd4, by = "patient") %>% 
    left_join(tb_rna, by = "patient") %>% 
    left_join(tbltb_res, by = "patient") %>% 
    left_join(treatment_art, by = "patient") %>% 
    left_join(tblVIS, by = "patient") %>% 
    mutate(cohort = as.factor("RSA"),
           region = as.factor("Sub-Saharan Africa"),
           incident_tb = as.factor(case_when(
             date_tb >= art_start_date + 60 ~ 1,
             TRUE ~ 0)),
           prevalent_tb = as.factor(case_when(
             date_tb < art_start_date + 60 & date_tb > art_start_date - 60 ~ 1,
             TRUE ~ 0)),
           recent_tb = as.factor(case_when(
             date_tb < art_start_date - 60 & date_tb > art_start_date - 360 ~ 1,
             TRUE ~ 0)),
           presenting_tb = as.factor(case_when(
             prevalent_tb == 1 | recent_tb == 1 ~ 1,
             TRUE ~ 0)),
           last_persontime = as.numeric(case_when(!is.na(exitdate) ~ exitdate - art_start_date,
                                       TRUE ~ last_fup_date - art_start_date))) %>% 
  select(-enrol_d) %>% 
  distinct(patient, .keep_all = TRUE) %>%
  rename(id = patient) %>% 
  filter(born < exitdate | is.na(exitdate))

#### Defining study populations ------------------------------------------------

#' Using tblBAS as a basis, using the inclusion criteria defined in the SAP
#' As mentioned by Chido, for the ART analysis, it is important to check, 
#' if the patient is ART-naive when entering the cohort. If not, this patient also 
#' has to be excluded. 
#' This can be made sure if patient is ART-naive upon enrollment and double checked
#' in comparing Date of first antiretroviral treatment initiation with enrollement
#' in Khayelisha cohort. 

## TB ##
  
df_tb <- df %>%
  filter(between(date_tb, as.Date("2010-01-01"), as.Date("2022-12-31")),
         sex != 9,
         !is.na(born) & !is.na(art_start_date),
         age_at_art_start >= 16,
         !is.na(date_tb))

saveRDS(df_tb, "data_clean/rsa/tb_rsa.rds")

## ART ##  

df_art <- df %>% 
  filter(id %in% tblNAIVE$patient) %>% 
  filter(between(art_start_date, as.Date("2010-01-01"), as.Date("2022-12-31")),
         sex != 9,
         !is.na(born) & !is.na(art_start_date),
         age_at_art_start >= 16) %>% 
  rename(gender = sex)
  
saveRDS(df_art, "data_clean/rsa/art_rsa.rds")

#flextable::flextable(tabyl(df_art$treatment))

#### Lab data ------------------------------------------------------------------

## cd4 ##

lab_cd4 <- tblLAB_CD4 %>% 
  filter(patient %in% df_art$id) %>% 
  left_join(df %>% select(id, art_start_date, disease_tb, date_tb, presenting_tb), 
            by = c("patient" = "id")) %>% 
  arrange(patient, date_cd4) %>% 
  group_by(patient) %>% 
  mutate(timepoint = row_number(),
         time_diff = as.numeric(date_cd4 - art_start_date, units = "days")) %>% 
  ungroup() %>% 
  rename(id = patient)

saveRDS(lab_cd4, "data_clean/rsa/cd4_rsa.rds")
  
## rna ##

lab_rna <- tblLAB_RNA %>%
  filter(patient %in% df_art$id) %>% 
  left_join(df %>% select(id, art_start_date, disease_tb, date_tb, presenting_tb), by = c("patient" = "id")) %>% 
  arrange(patient, date_rna) %>%
  group_by(patient) %>% 
  mutate(timepoint = row_number(),
    time_diff = as.numeric(date_rna - art_start_date, units = "days")) %>% 
  ungroup() %>% 
  rename(id = patient)
  
saveRDS(lab_rna, "data_clean/rsa/rna_rsa.rds")

