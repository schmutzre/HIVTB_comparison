#### Library -------------------------------------------------------------------

library(tidyverse)
library(lubridate)

swiss <- readRDS("data_clean/art_ch.rds")

#### Cleaning the separate datafile --------------------------------------------

## ART treatment ##

tblART <- read.csv("data_raw/RSA/tblART.csv") %>% 
  select(patient, art_id, art_combination, art_sd, art_ed)
# ask Lukas

## Base table ##

tblBAS <- read.csv("data_raw/RSA/tblBAS.csv") %>% 
  mutate(across(c(birth_d, enrol_d, recart_d), ~na_if(.x, "")),
         across(c(birth_d, enrol_d, recart_d,), ~as.Date(.x, format = "%Y-%m-%d")),
         age_at_art_start = interval(birth_d, recart_d) %/% years(1),
         sex = as.factor(case_when(sex == 1 ~ "Male",
                         sex == 2 ~ "Female",
                         TRUE ~ NA))) %>% 
  filter(!is.na(recart_d)) %>% 
  rename(born = birth_d,
         risk = mode, 
         art_start_date = recart_d,
         risk = mode) %>% 
  select(patient, born, sex, art_start_date, age_at_art_start, risk, naive_y)

## CD4 Lab ##

tblLAB_CD4 <- read.csv("data_raw/RSA/tblLAB_CD4.csv") 
#Lukas fragen wie verschiedene Masseinheiten überführen

tblLAB_RNA <- read.csv("data_raw/RSA/tblLAB_RNA.csv")%>% 
  mutate(rna = case_when(rna_v < 0 ~ rna_l,
                         TRUE ~ rna_v)) %>% 
  rename(date_rna = rna_d) %>% 
  select(patient, date_rna, rna)

## Exit and last fup date ##

tblLTFU <- read.csv("data_raw/RSA/tblLTFU.csv") %>% 
  select(patient, drop_d, death_d) %>% 
  rename(last_fup_date = drop_d,
         exitdate = death_d) %>% 
  mutate(across(c(last_fup_date, exitdate), ~na_if(.x, "")),
         across(c(last_fup_date, exitdate), ~as.Date(.x, format = "%Y-%m-%d")))

## TB resistance ##

tbltb_res <- read.csv("data_raw/RSA/tbltb_res.csv") %>%
  group_by(patient) %>%
  mutate(resistance_tb = as.factor(case_when(
    any(drug_res == 1) ~ 1,
    TRUE ~ 0
  ))) %>%
  ungroup() %>%
  distinct(patient, .keep_all = TRUE) %>%
  select(patient, resistance_tb)

## TB table ##

tblTB <- read.csv("data_raw/RSA/tblTB.csv") %>% 
  group_by(patient) %>%
  filter(reg_dmy == min(reg_dmy),
         startsWith(patient, "K"),
         cat %in% c(1,95,99,NA)) %>%
  ungroup() %>% 
  rename(date_tb = reg_dmy, 
         regimen_tb = regimen,
         site_tb = class) %>% 
  select(patient, date_tb, regimen_tb, site_tb) %>% 
  mutate(across(c(date_tb, site_tb), ~na_if(.x, "")),
         date_tb = as.Date(date_tb, format = "%Y-%m-%d"),
         disease_tb = as.factor(1),
         site_tb = as.factor(case_when(site_tb == 1 ~ "Pulmonary",
                             site_tb == 2 ~ "Extrapulmonary",
                             TRUE ~ NA)),
         regimen_tb = as.factor(case_when(regimen_tb == 1 ~ "2HRZE 4HR",
                                          regimen_tb == 2 ~ "2HRZES 1HRZE 5HRE",
                                          regimen_tb == 3 ~ "2HRZ 4HR",
                                          regimen_tb == 4 ~ "Others",
                                          TRUE ~ NA)))

## WHO stage ##

tblVIS <- read.csv("data_raw/RSA/tblVIS.csv") %>% 
  select(patient, vis_d, who_stage) %>% 
  left_join(tblBAS, by = "patient") %>%
  mutate(vis_d = as.Date(vis_d, format = "%Y-%m-%d")) %>% 
  group_by(patient) %>%
  mutate(difference = abs(vis_d - art_start_date)) %>%
  filter(difference <= 60) %>%
  arrange(patient, difference) %>%
  slice(1) %>%
  ungroup() %>%
  select(-difference) # use vis_d to select who_stage at ART start

### Code should run but takes to long for it not to have any values in it

#### Combining ART & TB dataset ------------------------------------------------

#' Similar to the SHCS dataset, I will first unify both study populations and then
#' in the last step, slice the dataset and produce two seperate datsets for the
#' ART study population and the TB study population.

df <- tblBAS %>% 
  left_join(tblTB, by = "patient") %>% 
  left_join(tblLTFU, by = "patient") %>% 
  mutate(cohort = "RSA",
         region = "Sub-Saharan Africa",
         prevalent_tb = case_when(
           date_tb < art_start_date + 60 & date_tb > art_start_date - 60 ~ 1,
           TRUE ~ 0),
         recent_tb = case_when(
           date_tb < art_start_date - 60 & date_tb > art_start_date - 360 ~ 1,
           TRUE ~ 0),
         presenting_tb = case_when(
           prevalent_tb == 1 | recent_tb == 1 ~ 1,
           TRUE ~ 0),
         last_persontime = case_when(!is.na(exitdate) ~ exitdate - art_start_date,
                                     TRUE ~ last_fup_date - art_start_date))
  











df_tb <- df %>%
  filter(between(date_tb, as.Date("2010-01-01"), as.Date("2022-12-31")),
         sex != 9,
         !is.na(born) & !is.na(art_start_date),
         age_at_art_start >= 16,
         !is.na(date_tb))

df_art <- df %>% 
  filter(between(art_start_date, as.Date("2010-01-01"), as.Date("2022-12-31")),
         sex != 9,
         !is.na(born) & !is.na(art_start_date),
         age_at_art_start >= 16)
         
#### Defining study population -------------------------------------------------

#' Using tblBAS as a basis, using the inclusion criteria defined in the SAP
#' As mentioned by Chido, it is important to check, if the patient is ART-naive 
#' when entering the cohort. If not, this patient also has to be excluded. 
#' This can be made sure if patient is ART-naive upon enrollment and double checked
#' in comparing Date of first antiretroviral treatment initiation with enrollement
#' in Khayelisha cohort. 
#' Afterwards, I'll use the Date of first antiretroviral treatment initiation as
#' ART_start_date.

#### Create base file ----------------------------------------------------------

#### Study population
  
df_naive <- tblBAS %>% 
  select(patient, birth_d, enrol_d, sex, mode, naive_y, recart_d) %>% 
  mutate(across(c(birth_d, enrol_d, recart_d), ~na_if(.x, ""))) %>%  # Replace empty strings with NA
  mutate(across(c(birth_d, enrol_d, recart_d), ~as.Date(.x, format = "%Y-%m-%d"))) %>% 
  filter(naive_y %in% c(1,9),
         recart_d >= enrol_d - 30) %>% 
  select(-c(enrol_d, naive_y, recart_d))
  

