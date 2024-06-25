#### Library -------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(janitor)
library(stringr)
library(flextable)
library(data.table)

source("utils/functions.R")

#### Cleaning the separate data files and preparing to then join them ----------

## Base table ------------------------------------------------------------------

#' This is the base of the whole database. Here, each patient has a row with their characteristics
#' This is the first broad definition of the study population based on the art start date

tblBAS.sa <- read.csv("data_raw/RSAnew/tblBAS.csv")  
tblBAS.sa2 <- tblBAS.sa %>%   
  mutate(RECART_D = as.character(RECART_D),
         across(c(BIRTH_D, ENROL_D, RECART_D), ~na_if(.x, "")),
         across(c(BIRTH_D, ENROL_D, RECART_D), ~as.Date(.x, format = "%Y-%m-%d")),
         age_at_art_start = interval(BIRTH_D, RECART_D) %/% years(1),
         sex = as.factor(case_when(SEX == 1 ~ "Male",
                                   SEX == 2 ~ "Female",
                                   TRUE ~ NA))) %>% 
  rename(born = BIRTH_D, 
         art_start_date = RECART_D,
         patient = PATIENT) %>%
  filter(between(art_start_date, as.Date("2017-01-01"), as.Date("2022-12-31")),
         !is.na(art_start_date),
         !is.na(born)) %>% 
  dplyr::select(patient, born, sex, art_start_date, age_at_art_start)

## ART treatment ---------------------------------------------------------------

#' I use this dataset first to assign the ART treatment but also to exclude patients that already received ART treatment
#' before the art_start_data which would mean they are not ART naive, which is an inclusion criteria for our study. 
#' This means, I will use this dataset now to further narrow down the study population.

tblART.sa <- read.csv("data_raw/RSAnew/tblART.csv") 
tblART.sa2 <- tblART.sa %>% 
  filter(patient %in% tblBAS.sa2$patient) %>% 
  mutate(across(c(art_sd, art_ed), ~as.Date(.x, format = "%Y-%m-%d"))) %>% 
  filter(between(art_sd, as.Date("2017-01-01"), as.Date("2022-12-31"))) %>% 
  dplyr::select(patient, art_id, art_sd, art_ed) %>% 
  left_join(tblBAS.sa2 %>% dplyr::select(patient, art_start_date), by = "patient") %>% 
  mutate(diff = art_sd - art_start_date)

notARTnaive <- tblART.sa2 %>% 
  filter(diff < 0)

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

tblART.sa3 <- tblART.sa2 %>% 
  # Step 1: Filter out patients that are not ART naive
  filter(patient %nin% notARTnaive$patient, diff < 180) %>%
  # Step 3: For each patient, keep only the row with the smallest diff in case of duplicates in the art_id column
  group_by(patient, art_id) %>%
  slice_min(order_by = diff, n = 1) %>%
  ungroup() %>%
  # Step 4: Join with id_to_drug to get drug information
  left_join(id_to_drug %>% dplyr::select(drug, art_id), by = c("art_id")) %>%
  group_by(patient) %>%
  arrange(diff, .by_group = TRUE) %>%
  # Step 5: Summarize treatment with conditional logic for diff = 0
  summarize(
    treatment = {
      drugs <- unique(na.omit(drug[diff == 0]))
      if (length(drugs) == 0) {
        drugs <- unique(na.omit(drug))
      }
      drugs <- setdiff(drugs, "ART unspecified")
      if (length(drugs) == 0) {
        drugs <- unique(na.omit(drug[diff == 0]))
        if (length(drugs) == 0) {
          drugs <- unique(na.omit(drug))
        }
      }
      paste(drugs, collapse = ", ")
    }
  ) %>%
  ungroup() %>% 
  mutate(
    treatment = case_when(
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "EFV")) ~ "TDF + 3TC/FTC + EFV",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "LPV")) ~ "TDF + 3TC/FTC + LPV/r",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "NVP")) ~ "TDF + 3TC/FTC + NVP",
      str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "LPV") ~ "AZT + 3TC/FTC + LPV/r",
      str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "AZT + 3TC/FTC + DTG",
      str_detect(treatment, "TDF") & str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "TDF + AZT + 3TC/FTC + DTG",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & str_detect(treatment, "DTG") ~ "TDF + 3TC/FTC + DTG",
      str_detect(treatment, "TDF") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "LPV") | str_detect(treatment, "ATV") | str_detect(treatment, "DTG")) ~ "TDF + 3TC/FTC + LPV/r or ATV/r or DTG",
      str_detect(treatment, "AZT") & (str_detect(treatment, "3TC") | str_detect(treatment, "FTC")) & (str_detect(treatment, "LPV") | str_detect(treatment, "ATV") | str_detect(treatment, "DTG")) ~ "AZT + 3TC/FTC + LPV/r or ATV/r or DTG",
      TRUE ~ str_replace_all(treatment, ",", " ++ ")  # Replace spaces with plus signs
    ),
    regimen = case_when(treatment %in% c("TDF + 3TC/FTC + EFV", "TDF + 3TC/FTC + NVP", "AZT + 3TC/FTC + DTG", "TDF + AZT + 3TC/FTC + DTG") ~ "NNRTI-based",
                        treatment %in% c("AZT + 3TC/FTC + LPV/r", "TDF + 3TC/FTC + LPV/r", "TDF + 3TC/FTC + DTG") ~ "InSTI-based",
                        treatment == "ART unspecified" ~ "Other",
                        !is.na(treatment) ~ "Other",
                        TRUE ~ NA)) %>% 
  mutate(treatment = case_when(
    str_detect(treatment, "\\+\\+") ~ "Other",  # Escape the + characters
    treatment == "ART unspecified" ~ "Other",
    TRUE ~ treatment
  ))

## Base table 2 ----------------------------------------------------------------

#' I exclude patients that were given ART treatment before art start date
#' I do this because this should always be the reference dataset which represents 
#' the study population.

tblBAS.sa3 <- tblBAS.sa2 %>% 
  filter(patient %nin% notARTnaive$patient)
  
## CD4 Lab ---------------------------------------------------------------------

tblLAB_CD4.sa <- read.csv("data_raw/RSAnew/tblLAB_CD4.csv") 
tblLAB_CD4.sa2 <- tblLAB_CD4.sa %>% 
  filter(CD4_U == 1) %>% 
  rename(cd4 = CD4_V,
         patient = PATIENT) %>% 
  mutate(date_cd4 = as.Date(CD4_D, format = "%Y-%m-%d")) %>% 
  dplyr::select(patient, date_cd4, cd4)

## RNA Lab ---------------------------------------------------------------------

tblLAB_RNA.sa <- read.csv("data_raw/RSAnew/tblLAB_RNA.csv") 
tblLAB_RNA.sa2 <- tblLAB_RNA.sa %>% 
  mutate(date_rna = as.Date(RNA_D, format = "%Y-%m-%d")) %>% 
  rename(rna = RNA_V,
         patient = PATIENT) %>% 
  dplyr::select(patient, date_rna, rna)

## TB table --------------------------------------------------------------------

tblTB.sa <- read.csv("data_raw/RSANEW/tblTB.csv") 
tblTB.sa2 <- tblTB.sa %>% 
  rename(date_tb = REG_D, 
         regimen_tb = REGIMEN,
         site_tb = CLASS,
         patient = PATIENT,
         resistant_tb= RESISTANT) %>% 
  mutate(across(c(date_tb), ~na_if(.x, "")),
         date_tb = as.Date(date_tb, format = "%Y-%m-%d"),
         disease_tb = 1,
         site_tb = as.factor(case_when(site_tb %in% c(1,3) ~ "Pulmonary",
                                       site_tb == 2 ~ "Extrapulmonary",
                                       TRUE ~ NA)),
         regimen_tb = as.factor(case_when(regimen_tb == 1 ~ "2HRZE 4HR",
                                          regimen_tb == 2 ~ "2HRZES 1HRZE 5HRE",
                                          regimen_tb == 3 ~ "2HRZ 4HR",
                                          regimen_tb == 4 ~ "Others",
                                          TRUE ~ NA)),
         outcome_tb = case_when(TB_OUTCOME %in% c(1,2) ~ "Sucessfull",
                                TB_OUTCOME == 3 ~ "Failed",
                                TB_OUTCOME == 4 ~ "Interrupted",
                                TB_OUTCOME == 5 ~ "Lost to follow-up",
                                TB_OUTCOME == 6 ~ "Treatment ongoing",
                                TB_OUTCOME == 7 ~ "Died",
                                TRUE ~ NA)) %>% 
  group_by(patient) %>% 
  arrange(date_tb)

tblTB.sa3 <- tblTB.sa2 %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  ungroup() %>% 
  mutate(
  regimen_tb = as.factor(regimen_tb),
  regimen_tb_group = as.factor(case_when(regimen_tb == "2HRZE 4HR" ~ "Standard",
                                         TRUE ~ regimen_tb))) %>%
dplyr::select(patient, date_tb, regimen_tb, disease_tb, regimen_tb_group, site_tb, outcome_tb, resistant_tb)

## below I define prevalent ect. this may be problematic when I only consider the first occurence of TB?

tblTB.sa.count <- tblTB.sa2 %>% 
  summarize(count_TB = n()) 

## WHO stage table -------------------------------------------------------------

tblVIS.sa <- read.csv("data_raw/RSAnew/tblVIS.csv")

tblVIS.sa2 <- tblVIS.sa %>% 
  filter(!is.na(WHO_STAGE)) %>% 
  dplyr::select(PATIENT, VIS_D, WHO_STAGE) %>% 
  filter(PATIENT %in% tblBAS.sa$patient) %>% 
  rename(patient = PATIENT,
         who_stage = WHO_STAGE) %>% 
  mutate(vis_d = as.Date(VIS_D, format = "%Y-%m-%d")) %>% 
  left_join(tblBAS.sa %>% dplyr::select(patient, art_start_date), by = "patient")  %>% 
  mutate(difference = abs(vis_d - art_start_date),
         who_stage = as.factor(
           case_when(
             who_stage == 9 ~ NA,
             TRUE ~ as.character(who_stage))),
         who_stage = as.factor(case_when(who_stage %in% c(1,2)  ~ "1/2", 
                                         who_stage %in% c(3,4) ~ "3/4",
                                         TRUE ~ NA))) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  select(patient, who_stage)

## Exit and last fup date ------------------------------------------------------

#' this is the official ltfu table which includes dates of drop out but not date of
#' last visit for everyone

tblLTFU.sa <- read.csv("data_raw/RSAnew/tblLTFU.csv") 
tblLTFU.sa2 <- tblLTFU.sa %>% 
  dplyr::select(PATIENT, DROP_D, DEATH_D, DROP_Y, DEATH_Y, L_ALIVE_D) %>% 
  rename(drop_out = DROP_D,
         exitdate = DEATH_D,
         patient = PATIENT,
         ltfu = DROP_Y, 
         last_info = L_ALIVE_D) %>% 
  mutate(across(c(drop_out, exitdate, last_info), ~na_if(.x, "")),
         across(c(drop_out, exitdate, last_info), ~as.Date(.x, format = "%Y-%m-%d")))

## CD4 baseline values ---------------------------------------------------------

baseline_cd4.sa <- tblLAB_CD4.sa2 %>% 
  filter(patient %in% tblBAS.sa3$patient) %>% 
  left_join(tblBAS.sa3 %>% dplyr::select(patient, art_start_date), 
            by = "patient") %>% 
  mutate(cd4_baseline = as.numeric(ifelse(date_cd4 >= (art_start_date - 180) & date_cd4 <= (art_start_date + 30) & cd4 <= 2000, 
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
  dplyr::select(patient, cd4_baseline, cd4_group)

## RNA baseline values ##

baseline_rna.sa <- tblLAB_RNA.sa2 %>% 
  left_join(tblBAS.sa3 %>% dplyr::select(patient, art_start_date), 
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
  dplyr::select(patient, rna_baseline, rna_group)

## CD4 @ TB diagnosis ----------------------------------------------------------

tb_cd4.sa <- tblLAB_CD4.sa2 %>% 
  left_join(tblTB.sa2 %>% dplyr::select(patient, date_tb), 
            by = "patient") %>% 
  mutate(tb_diag_cd4 = ifelse(date_cd4 >= (date_tb - 60) & date_cd4 <= (date_tb + 60), 
                              cd4, NA)) %>% 
  filter(!is.na(tb_diag_cd4)) %>% 
  arrange(abs(date_tb - date_cd4)) %>% 
  group_by(patient) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  ungroup() %>% 
  dplyr::select(patient, tb_diag_cd4)

## RNA @ TB diagnosis ----------------------------------------------------------

tb_rna.sa <- tblLAB_RNA.sa2 %>% 
  left_join(tblTB.sa2 %>% dplyr::select(patient, date_tb), 
            by = "patient") %>% 
  mutate(tb_diag_rna = as.numeric(ifelse(date_rna >= (date_tb - 60) & date_rna <= (date_tb + 60), 
                                         rna, NA))) %>% 
  filter(!is.na(tb_diag_rna)) %>% 
  arrange(abs(date_tb - date_rna)) %>% 
  group_by(patient) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  ungroup() %>% 
  dplyr::select(patient, tb_diag_rna)

#### Joining them together -----------------------------------------------------

df.sa <- tblBAS.sa3 %>% 
  left_join(tblTB.sa3, by = "patient") %>% 
  left_join(tblLTFU.sa2, by = "patient") %>% 
  left_join(baseline_cd4.sa, by = "patient") %>% 
  left_join(baseline_rna.sa, by = "patient") %>% 
  left_join(tb_cd4.sa, by = "patient") %>% 
  left_join(tb_rna.sa, by = "patient") %>% 
  left_join(tblVIS.sa2, by = "patient") %>% 
  left_join(tblART.sa3, by = "patient") %>% 
  mutate(cohort = as.factor("RSA"),
         region = as.factor("Sub-Saharan Africa"),
         who_stage = as.factor(who_stage),
         ltfu = as.numeric(ltfu),
         resistant_tb = as.numeric(resistant_tb),
         tb_diag_cd4 = as.numeric(tb_diag_cd4),
         DEATH_Y = as.numeric(DEATH_Y),
         
         incident_tb = as.factor(case_when(
           date_tb >= art_start_date + 60 ~ 1,
           TRUE ~ 0)),
         presenting_tb = as.factor(case_when(
           date_tb <= art_start_date + 60 & date_tb >= art_start_date - 60 ~ 1,
           TRUE ~ 0)),
         recent_tb = as.factor(case_when(
           date_tb < art_start_date - 60 & date_tb > art_start_date - 360 ~ 1,
           TRUE ~ 0)),
         fup_time = as.numeric(case_when(
           ltfu == 1 ~ drop_out - art_start_date,
           DEATH_Y == 1 ~ exitdate - art_start_date,
           TRUE ~ pmin(as.Date("2022-12-31") - art_start_date, last_info - art_start_date))),
         disease_tb = ifelse(is.na(date_tb),0,1)) %>% 
  rename(id = patient,
         death = DEATH_Y) %>% 
  filter(born < exitdate | is.na(exitdate))

## Data check fup_time ---------------------------------------------------------
#' for those that have a fup_time of NA, I will check if they have a last visit date
#' if they do, I will use this as the fup_time. 
#' I decide to just filter them out after all, they do not provide any information
#' and most of them had their last visit way before art start

# tblART.na <- df.sa %>% 
#   filter(is.na(fup_time)) 
# 
# tblVIS_check <- tblVIS.sa %>% 
#   filter(PATIENT %in% tblART.na$id) %>%
#   group_by(PATIENT) %>%
#   arrange(VIS_D) %>%
#   distinct(PATIENT, .keep_all = TRUE) %>% 
#   select(PATIENT, VIS_D) %>% 
#   mutate(VIS_D = as.Date(VIS_D, format = "%Y-%m-%d")) %>% 
#   rename(id = PATIENT)
# 
# df.sa2 <- df.sa %>%
#   left_join(tblVIS_check, by = "id") %>%
#   mutate(fup_time = case_when(
#     is.na(fup_time) & !is.na(VIS_D) ~ as.numeric(difftime(VIS_D, art_start_date, units = "days")),
#     TRUE ~ fup_time)) 
#   dplyr::select(-VIS_D)

## Defining study populations --------------------------------------------------

## ART ##  

df_art.sa <- df.sa %>% 
  filter(between(art_start_date, as.Date("2017-01-01"), as.Date("2022-12-31")),
         !is.na(sex),
         !is.na(born) & !is.na(art_start_date),
         !is.na(fup_time) & fup_time >= 0,
         age_at_art_start >= 16,
         age_at_art_start < 100,
         recent_tb == 0,
         (date_tb >= art_start_date - 60 | is.na(date_tb))) %>% 
  mutate(who_stage = as.factor(case_when(presenting_tb == 1 ~ "3/4",
                               TRUE ~ who_stage))) %>% 
  rename(gender = sex) %>% 
  dplyr::select(-ENROL_D, -drop_out, -exitdate, -last_info)

saveRDS(df_art.sa, "data_clean/rsa/art_rsa.rds")

art_rsa_noTB <- df_art.sa %>% 
  filter(presenting_tb == 0)

saveRDS(art_rsa_noTB, "data_clean/rsa/art_noTB_rsa.rds")

## Data summary ----------------------------------------------------------------


df_art.sa2 <- df_art.sa %>% 
  mutate(presenting_tb = ifelse(presenting_tb ==0, "No", "Yes"))

df_art.sa2$treatment <- factor(df_art.sa2$treatment, levels = c("TDF + 3TC/FTC + EFV",
                                                     "TDF + 3TC/FTC + DTG",
                                                     "AZT + 3TC/FTC + LPV/r",
                                                     "TDF + 3TC/FTC + NVP",
                                                     "TDF + 3TC/FTC + LPV/r",
                                                     "TDF + 3TC/FTC + LPV/r or ATV/r or DTG",
                                                     "AZT + 3TC/FTC + DTG",
                                                     "AZT + 3TC/FTC + LPV/r or ATV/r or DTG",
                                                     "Other"))

tbl_overall <- df_art.sa2 %>% 
  tbl_summary(
    include = c(
      `gender`, 
      `age_at_art_start`, 
      `cd4_baseline`, 
      `rna_baseline`,
      `who_stage`,
      `fup_time`,
      `treatment`),
    digits = list(everything() ~ c(0)))  

tbl1B <- df_art.sa2 %>% 
  tbl_summary(by = presenting_tb,
    include = c(
      `gender`, 
      `age_at_art_start`, 
      `cd4_baseline`, 
      `rna_baseline`,
      `who_stage`,
      `fup_time`,
      `treatment`),
    digits = list(everything() ~ c(0))) 

tbl_merge <- tbl_merge(
  tbls = list(tbl_overall, tbl1B),
  tab_spanner = c("Overall", "presenting TB"))
tbl_merge
gt_tbl <- as_gt(tbl_merge)
#gtsave(gt_tbl, "results/table1.html")

## Formulations for Methods section --------------------------------------------

#' Person-time was calculated from ART initiation (baseline) to TB diagnosis, 
#' death or drop-out or last follow-up information or end of study period. 
#' Whatever came first. 

## Lab data --------------------------------------------------------------------

#### Lab data ------------------------------------------------------------------

## cd4 ##

lab_cd4 <- tblLAB_CD4.sa2 %>% 
  filter(patient %in% df_art.sa$id) %>% 
  left_join(df_art.sa %>% dplyr::select(id, art_start_date, disease_tb, date_tb, presenting_tb, gender, age_at_art_start, cd4_baseline, who_stage, regimen), 
            by = c("patient" = "id")) %>% 
  group_by(patient) %>% 
  arrange(date_cd4) %>%
  mutate(time_diff = as.numeric(date_cd4 - art_start_date, units = "days")) %>% 
  filter(time_diff >= -60) %>% 
  mutate(timepoint = row_number()) %>% 
  ungroup() %>% 
  rename(id = patient) %>% 
  filter(id %in% df_art.sa$id)

saveRDS(lab_cd4, "data_clean/rsa/cd4_rsa.rds")
  
lab_cd4 %>% 
  group_by(id) %>%
  summarise(time_point = max(timepoint)) %>% 
  ggplot2::ggplot(aes(x = time_point)) +
  geom_histogram() +
  labs(x = "Number of CD4 measurements",
       y = "Number of patients")

## rna ##

lab_rna <- tblLAB_RNA.sa2 %>%
  filter(patient %in% df_art.sa$id) %>% 
  left_join(df_art.sa %>% dplyr::select(id, art_start_date, disease_tb, date_tb, presenting_tb), by = c("patient" = "id")) %>% 
  group_by(patient) %>% 
  arrange(date_rna) %>%
  mutate(time_diff = as.numeric(date_rna - art_start_date, units = "days")) %>% 
  filter(time_diff >= -60) %>% 
  mutate(timepoint = row_number()) %>% 
  ungroup() %>% 
  rename(id = patient)

lab_rna %>% 
  group_by(id) %>%
  summarise(time_point = max(timepoint)) %>% 
  ggplot2::ggplot(aes(x = time_point)) +
  geom_histogram() +
  labs(x = "Number of RNA measurements",
       y = "Number of patients")

saveRDS(lab_rna, "data_clean/rsa/rna_rsa.rds")


