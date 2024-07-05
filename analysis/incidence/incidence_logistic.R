#### libraries -----------------------------------------------------------------

library(tidyverse)
library(jtools)
library(tableone)
library(mice)
library(jtools)
library(brglm2)
library(gtsummary)
library(purrr)

#### data preparation ---------------------------------------------------------

data <- readRDS("data_clean/art_noTB.rds")
custom_breaks <- c(16, 34, 44, 100)

df <- data %>% 
  rename("sex" = gender, 
         "age" = age_at_art_start,
         "bcd4" = cd4_baseline,
         "brna" = rna_baseline,
         "who" = who_stage,
         "regio" = region,
         "pres_tb" = presenting_tb) %>%
  mutate(bcd4_tr = sqrt(bcd4),
         brna_tr = log10(brna + 1),
         regimen = as.factor(regimen),
         regio = relevel(regio, ref = "Europe/Northern America"),
         agegroup = cut(age, breaks = custom_breaks, include.lowest = TRUE),
         agegroup = as.factor(agegroup)) %>% 
  dplyr::select(sex, agegroup, regio, bcd4_tr, brna_tr, who, regimen, incident_tb, cohort)

df_ch <- df %>% filter(cohort == "CH") %>%
  dplyr::select(-cohort) 

df_rsa <- df %>% filter(cohort == "RSA") %>% 
  dplyr::select(-cohort)

#### complete record analysis --------------------------------------------------

#' The complete record analysis is valid under the assumption of probability 
#' of complete record not associated with outcome (incidence of TB) given covariates.

### CH ###

m1_ch <- glm(incident_tb ~ sex + 
               agegroup + 
               who +
               regio +  bcd4_tr + brna_tr, 
                 family = binomial(link = "logit"), 
                 data = df_ch, method = "brglmFit")


tbl_regCH.complete <- m1_ch |>
  tbl_regression(exponentiate = TRUE,
                 add_estimate_to_reference_rows = TRUE,
                 label = list(sex = "Sex",
                              agegroup = "Age group",
                              who = "WHO stage",
                              regio = "Region",
                              bcd4_tr = "CD4 count (sqrt)",
                              brna_tr = "RNA (log10)")) %>% 
  bold_labels()


### RSA ###

m1_rsa <- glm(incident_tb ~ sex + agegroup + bcd4_tr + who, 
                  family = binomial(link = "logit"), 
                  data = df_rsa, method = "brglmFit")

tbl_regRSA.complete <- m1_rsa |>
  tbl_regression(exponentiate = TRUE,
                 add_estimate_to_reference_rows = TRUE,
                 label = list(sex = "Sex",
                              agegroup = "Age group",
                              who = "WHO stage",
                              bcd4_tr = "CD4 count (sqrt)")) %>% 
  bold_labels()

#### imputed record analysis ---------------------------------------------------

### CH ###

ch.imputed <- readRDS("data_clean/imputed/ch_imp.rds")

fit_ch_imp <- with(data = ch.imputed, 
                   exp = glm(incident_tb ~ 
                               sex + age + who + regio + bcd4_tr + brna_tr,
                             family = binomial(link = "logit"), method = "brglmFit"))

tbl_regCH.imputed <- fit_ch_imp |>
  tbl_regression(exponentiate = TRUE,
                 add_estimate_to_reference_rows = TRUE,
                 label = list(sex = "Sex",
                              age = "Age group",
                              who = "WHO stage",
                              regio = "Region",
                              bcd4_tr = "CD4 count (sqrt)",
                              brna_tr = "RNA (log10)")) %>% 
  bold_labels()
  
### RSA ###

rsa.imputed <- readRDS("data_clean/imputed/rsa_imp.rds")

fit_rsa_imp <- with(data = rsa.imputed,
                    exp = glm(incident_tb ~ sex + age + bcd4_tr + who, 
                              family=binomial(link = "logit"), method = "brglmFit"))

tbl_regRSA.imputed <- fit_rsa_imp |>
  tbl_regression(exponentiate = TRUE,
                 add_estimate_to_reference_rows = TRUE,
                 label = list(sex = "Sex",
                              age = "Age group",
                              who = "WHO stage",
                              bcd4_tr = "CD4 count (sqrt)")) %>% 
  bold_labels()

#### combine results ------------------------------------------------------------

tbl.log.ch <- tbl_merge(
  tbls = list(tbl_regCH.complete, tbl_regCH.imputed),
  tab_spanner = c("**Complete**", "**Imputed**")) 

tbl.log.rsa <- tbl_merge(
  tbls = list(tbl_regRSA.complete, tbl_regRSA.imputed),
  tab_spanner = c("**Complete**", "**Imputed**")) 

saveRDS(tbl.log.ch, file = "results/incidenceTB/tbl_log_ch.rds")
saveRDS(tbl.log.rsa, file = "results/incidenceTB/tbl_log_rsa.rds")
