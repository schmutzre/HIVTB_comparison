#### libraries -----------------------------------------------------------------

library(tidyverse)
library(jtools)
library(tableone)
library(mice)
library(jtools)
library(brglm2)
library(gtsummary)

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

#### investigate systematic missingness ----------------------------------------

data.to.impute_ch <- df_ch 

data.to.impute_rsa <- df_rsa %>%
  dplyr::select(-regio, -brna_tr)

syst_miss <- function(df_ti, return) {
  
  #' the function prints the whole systematic missingness analysis and adds the incomplete column to
  #' the respective dataset for the glm
  
  data.to.impute <- df_ti
  
  data.to.impute$incomplete <- complete.cases(data.to.impute)
  
  vars <- names(data.to.impute)
  
  cts <- vars[sapply(data.to.impute, is.numeric)]
  
  cat <- vars[(!(vars %in% cts))]
  
  tablemiss <- CreateTableOne(vars = vars, 
                              strata="incomplete", 
                              data = data.to.impute, 
                              factorVars = cat, 
                              includeNA = F)
  
  if (return == "data") {
    return(data.to.impute)
  } else if (return == "table") {
    return(tablemiss)
  } else if (return == "md") {
    return(md.pattern(data.to.impute, rotate.names = T))
  }
  else {
    return(list(dataframe1 = data.to.impute, dataframe2 = tablemiss))
  }
  
}

### CH ###

data.to.impute_ch <- syst_miss(data.to.impute_ch, return = "data")
syst_miss(data.to.impute_ch, return = "table")
syst_miss(data.to.impute_ch, return = "md")

miss.mod_ch <- glm(incomplete ~ 
                     sex +
                     agegroup + 
                     regio + 
                     bcd4_tr + 
                     brna_tr + 
                     #who + 
                     #incident_tb +
                     regimen, 
                   family= binomial(link = "logit"), 
                   data = data.to.impute_ch,
                   control = list(maxit = 50))

summ(miss.mod_ch, exp = T)

data.to.impute_ch$incomplete <- NULL

### RSA ### 

data.to.impute_rsa <- syst_miss(data.to.impute_rsa, return = "data")
syst_miss(data.to.impute_rsa, return = "table")
syst_miss(data.to.impute_rsa, return = "md") 

miss.mod_rsa <- glm(incomplete ~ 
                      sex +
                      agegroup + 
                      bcd4_tr +
                      incident_tb +
                      regimen + 
                      incident_tb, 
                    family= binomial(link = "logit"), 
                    data = data.to.impute_rsa, 
                    control = list(maxit = 50))

summ(miss.mod_rsa)

data.to.impute_rsa$incomplete <- NULL

#### imputation ----------------------------------------------------------------

#' here i include variables that are included in the scientific model 
#' (see complete record analysis) and variables that were significantly related 
#' to missingness in the previous analysis.

K <- 5 # number of imputed datasets

### CH ###

pmat_ch <- matrix(
  c(0,1,1,1,1,1,1,1, # sex
    1,0,1,1,1,1,1,1, # agegroup
    1,1,0,1,1,1,1,1, # regio
    1,1,1,0,1,1,1,1, # bcd4_tr
    1,1,1,1,0,1,1,1, # brna_tr
    1,1,1,1,1,0,1,1, # who
    1,1,1,1,1,1,0,1, # regimen
    1,1,1,1,1,1,1,0), #incident_tb
  nrow = 8,
  byrow = TRUE,
  dimnames = list(c("sex", "agegroup", "regio", "bcd4_tr", "brna_tr", "who", "regimen", "incident_tb"),
                  c("sex", "agegroup", "regio", "bcd4_tr", "brna_tr", "who", "regimen", "incident_tb"))
)

# shell imputation (for post processing)

ini_ch <- mice(data = data.to.impute_ch, maxit = 0)
post_ch <- ini_ch$post
post_ch["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post_ch["brna_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 10))"

# imputation

ch.imp <- mice(data = data.to.impute_ch,
               m = K,
               method = c("",
                          "",
                          "polyreg",  # regio
                          "norm",  # bcd4_tr
                          "norm", #brna_tr 
                          "logreg", # who
                          "polyreg", # regimen
                          ""), 
               predictorMatrix = pmat_ch,
               seed = 1569,
               maxit = 20,
               post = post_ch)

# check imputations
summary(complete(ch.imp, 0)) # with missingness, index 0
summary(complete(ch.imp, 1)) # imputed index 1:K
plot(ch.imp)

### RSA ###

pmat_rsa <- matrix(
  c(0,1,1,1,1,0, # sex
    1,0,1,1,1,0, # agegroup
    1,1,0,1,1,0, # bcd4_tr
    1,1,1,0,1,0, # who
    1,1,1,1,0,0, # regimen
    1,1,1,1,1,0), # incident_tb
  nrow = 6, 
  byrow = TRUE,
  dimnames = list(c("sex", "agegroup", "bcd4_tr", "who", "regimen", "incident_tb"),
                  c("sex", "agegroup", "bcd4_tr", "who", "regimen", "incident_tb")))

# shell imputation (for post processing) #

ini_rsa <- mice(data = data.to.impute_rsa, 
                maxit = 0)
post_rsa <- ini_rsa$post
post_rsa["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post_ch["bcd4_tr"] <- paste(post_ch["bcd4_tr"], "bcd4 <- imp[[j]][, i]^2; imp[[j]][, i] <- cut(bcd4, breaks = c(-Inf, 99, 349, Inf), labels = c('0-99', '100-349', '350+'))", sep = "; ")

# imputation #

rsa.imp <- mice(data = data.to.impute_rsa,
                m = K,
                method = c("", #sex
                           "", #agegroup
                           "norm",  # bcd4_tr
                           "logreg", #who
                           "polyreg", #regimen
                           ""), #incident_tb
                predictorMatrix = pmat_rsa,
                seed = 1569,
                maxit = 20,
                post = post_rsa)

# check imputations
summary(complete(rsa.imp, 0)) # with missingness, index 0
summary(complete(rsa.imp, 1)) # imputed index 1:K

#### imputed record analysis ---------------------------------------------------

### CH ###

fit_ch_imp <- with(data = ch.imp, 
                   exp = glm(incident_tb ~ 
                               sex + agegroup + who + regio + bcd4_tr + brna_tr,
                             family = binomial(link = "logit"), method = "brglmFit"))

tbl_regCH.imputed <- fit_ch_imp |>
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

fit_rsa_imp <- with(data = rsa.imp,
                    exp = glm(incident_tb ~ sex + agegroup + bcd4_tr + who, 
                              family=binomial(link = "logit"), method = "brglmFit"))

tbl_regRSA.imputed <- fit_rsa_imp |>
  tbl_regression(exponentiate = TRUE,
                 add_estimate_to_reference_rows = TRUE,
                 label = list(sex = "Sex",
                              agegroup = "Age group",
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
