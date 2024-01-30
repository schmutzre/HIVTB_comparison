#### Library ------------------------------------------------------------------- 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  haven,
  janitor,
  gridExtra,
  stringr,tidyverse
)

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

#' The followng df are the study populations. Once patients starting ART between 2010-2022 and once patients having a TB diagnosis 
#' between 2010-2022. These dataframes are mostly overlapping. 

#### Filter studypopulation (Union of HIV- and TB-dataset) ---------------------

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
    starts_with("tbd_drug_resist_"), starts_with("tbd_antim_resist_"), regdate
  ) %>%
  mutate(art_start_date = case_when(!is.na(haart_start_date) ~ haart_start_date,
                                    TRUE ~ art_start_date)) %>% 
  dplyr::select(-haart_start_date) %>% 
  # Step 3: Creating new variables
  mutate(
    cohort = as.factor("CH"),
    age_at_art_start = year(art_start_date) - born,
    prevalent_tb = as_factor(case_when(
      date_tb < art_start_date + 60 & date_tb > art_start_date - 60 ~ 1,
      TRUE ~ 0)),
    recent_tb = as_factor(case_when(
      date_tb < art_start_date - 60 & date_tb > art_start_date - 360 ~ 1,
      TRUE ~ 0)),
    presenting_tb = prevalent_tb) %>% 
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
      region %in% c("Southern Europe", "Western Europe", "Northern Europe", "Eastern Europe") ~ "Europe",
      region %in% c("Eastern Africa", "Middle Africa", "Southern Africa", "Western Africa") ~ "Sub-Saharan Africa",
      region %in% c("Eastern Asia", "South-Eastern Asia", "Southern Asia", "Western Asia") ~ "Asia",
      TRUE ~ region
    ))
  )

#### Type of infection ---------------------------------------------------------

filteredBOTH <- as_factor(filteredBOTH, levels="labels") %>% 
  mutate(risk = as.factor(case_when(risk == "Unknown" ~ NA,
                          TRUE ~ risk)))

#### Birth country / nationality -----------------------------------------------

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
                                      TRUE ~ 0)) %>% 
  dplyr::select(-region, -tbd_pat_birth) %>% 
  mutate(region = as.factor(region_born)) %>% 
  dplyr::select(-region_born)

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
                dplyr::select(id,cd4_baseline, labdate), by = "id") %>% 
    left_join(baseline_rna %>% 
                dplyr::select(id,rna_baseline, labdate), by = "id") %>% 
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
      cdc_group == "C" ~ 4
    ))) %>% 
    ungroup() 
  
filteredBOTH.who <- filteredBOTH.lab %>% 
    left_join(who_stages %>% dplyr::select(id, who_stage, cdc_group), by = "id") 
  
#### CD4 and VL @ TB-date ------------------------------------------------------

# Create tb_cd4 dataset
tb_cd4 <- filteredBOTH.who %>%
    left_join(lab, by = "id") %>%
    mutate(
      tb_diag_cd4 = ifelse(labdate >= (date_tb - 120) & labdate <= (date_tb + 120), cd4, NA)
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
      tb_diag_rna = ifelse(labdate >= (date_tb - 120) & labdate <= (date_tb + 120), rna, NA)
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
  left_join(admin %>% dplyr::select(id, last_fup_date), by = "id") %>% 
  mutate(last_persontime = case_when(
    !is.na(exitdate) ~ as.numeric(difftime(exitdate, art_start_date, units = "days")),
    TRUE ~ as.numeric(difftime(last_fup_date, art_start_date, units = "days"))
  ))
#' last persontime is defined as persontime until death or loss to follow up. For any incidence models (eg. TB or viral suppression)
#' the person years have to be re-calculated as I haven't incorporated persontime until incidence of any sort.

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
      TRUE ~ "Other"
    )) %>% 
  dplyr::select(-(num_art:art_start_date.y)) %>% 
  rename(art_start_date = art_start_date.x)

flextable::flextable(tabyl(filteredBOTH.regimen$treatment))

filteredBOTH.regimen <- filteredBOTH.regimen %>% 
  mutate(current_art = case_when(current_art == "" ~ NA,
                                 TRUE ~ current_art))

## this has to be done afterwards -->
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
  dplyr::select(id, date_tb, Rifampicin_R:Others_R) %>%
  pivot_longer(cols = Rifampicin_R:Others_R, names_to = "Drug", values_to = "value") %>% 
  filter(value == 1) %>%
  group_by(id) %>%
  summarise(
    resistance_tb = paste(gsub("_R$", "", Drug), collapse = ", ")
  ) %>%
  ungroup()

filteredBOTH.resist <- filteredBOTH.regimen %>% 
  left_join(tb_resi , by = "id") %>% 
  dplyr::select(-c(tbd_drug_resist___1:tbd_drug_resist_others)) %>% 
  mutate(resistance_tb_any = as.factor(case_when(!is.na(resistance_tb) ~ 1,
                   TRUE ~ 0)),
         resistance_tb_mdr = as.factor(case_when(
           (str_detect(resistance_tb, "Rifampicin")| str_detect(resistance_tb, "Rifabutin")) & str_detect(resistance_tb, "Isoniazid") ~ 1,
           TRUE ~ 0
         )))

#### TB treatment --------------------------------------------------------------

test <- complete_ch %>% 
  dplyr::select(tbd_antim_resist___1:tbd_antim_resist_others)
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
  mutate(disease_tb = as.factor(disease_tb),
         who_stage = as.factor(who_stage),
         id = as.character(id))

#### lab data -----------------------------------------------------------

lab_both <- lab.filtered %>% 
  arrange(id, labdate) %>% 
  group_by(id) %>% 
  mutate(timepoint = row_number()) %>% 
  ungroup() %>% 
  dplyr::select(-cd4date) %>%
  left_join(filteredBOTH.tbsite %>% dplyr::select(id, art_start_date, disease_tb, date_tb, presenting_tb), by = "id") %>% 
  mutate(time_diff = as.numeric(labdate - art_start_date, units = "days"),
         pre_2016 = as.factor(case_when(art_start_date <= as.Date("2016-12-31") ~ 1,
                                        TRUE ~ 0)))

lab_cd4 <- lab_both %>% 
  dplyr::select(-rna, -timepoint) %>% 
  filter(!is.na(cd4)) %>% 
  mutate(timepoint = row_number()) %>% 
  rename(date_cd4 = labdate)

## check cd4 over 2,000 ##
#' I already prevented the baseline cd4 to be over 2000, now i look at the 
#' non-baseline cases individually

ids_over2k <- lab_cd4 %>%
  filter(cd4 > 2000) %>%
  distinct(id)

check <- lab_cd4 %>%
  filter(id %in% ids_over2k$id)

plot_group <- function(ids, data) {
  group_data <- data %>% filter(id %in% ids$id)
  ggplot(group_data, aes(x = time_diff, y = cd4)) +
    geom_point() +
    geom_line() +
    facet_wrap(~id, scales = "free") +
    theme_minimal()
}

group1_idsCH <- ids_over2k[1:20,]
group2_idsCH <- ids_over2k[21:40,]
group3_idsCH <- ids_over2k[41:57,]

plot_group(group1_idsCH, lab_cd4)
plot_group(group2_idsCH, lab_cd4)
plot_group(group3_idsCH, lab_cd4)

jumpy <- 
  c(17948, 19363, 19581, 19598, 26816, 27172, 27222, 27277, 27325, 27427, 28529, 
  28542, 31892, 32297, 32571, 43015, 46711, 52425, 52504, 60720, 60867, 90307, 
  90370, 90566, 91355, 57464, 42704, 42876, 34025)

lab_cd4 <- lab_cd4 %>% 
  filter(cd4 <= 2000 | id %nin% jumpy)

lab_rna <- lab_both %>% 
  dplyr::select(-cd4, -timepoint) %>% 
  filter(!is.na(rna)) %>% 
  mutate(timepoint = row_number()) %>% 
  rename(date_rna = labdate)

#### Selection study population ------------------------------------------------

## ART ##

art_ch <- final %>% 
  filter(art_start_date <= as.Date("2022-12-31") & art_start_date >= as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date,
         eligibility_art ==1,
         recent_tb == 0,
         (date_tb >= art_start_date - 60 | is.na(date_tb))) %>% 
  mutate(incident_tb = as.factor(case_incident_2m),
         pre_2016 = as.factor(case_when(art_start_date <= as.Date("2016-12-31") ~ 1,
                            TRUE ~ 0))) %>% 
  dplyr::select(-virus_type, -type_tb_shcs, -disease_tbc, -risk, -art_start_cd4,-eligibility_art, -case_incident_2m, -moddate, -enddate, - time_diff_ART, -time_diff_STOP, -resistance_tb, -labdate_cd4, -labdate_rna, -current_art, -regdate, -cdc_group) %>% 
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

