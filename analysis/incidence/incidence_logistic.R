#### libraries -----------------------------------------------------------------

library(tidyverse)
library(jtools)
library(tableone)
library(mice)

#### data preparation ---------------------------------------------------------

data <- readRDS("data_clean/art_noTB.rds")

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
         regio = relevel(regio, ref = "Europe/Northern America")) %>% 
  dplyr::select(sex, age, regio, bcd4_tr, brna_tr, who, regimen, incident_tb, cohort)

df_ch <- df %>% filter(cohort == "CH") %>%
  dplyr::select(-cohort) 

df_rsa <- df %>% filter(cohort == "RSA") %>% 
  dplyr::select(-cohort)

#### complete record analysis --------------------------------------------------

#' The complete record analysis is valid under the assumption of probability 
#' of complete record not associated with outcome (incidence of TB) given covariates.

### CH ###

m1_ch <- glm(incident_tb ~ sex + age + who + regio +  bcd4_tr + brna_tr, 
                 family = "binomial", 
                 data = df_ch)

summ(m1_ch, exp = T)

### RSA ###

m1_rsa <- glm(incident_tb ~ sex + age + bcd4_tr, 
                  family = "binomial", 
                  data = df_rsa)

summ(m1_rsa, exp = T)

#### investigate systematic missingness ----------------------------------------

data.to.impute_ch <- df_ch 

data.to.impute_rsa <- df_rsa %>%
  dplyr::select(-regio, - brna_tr, -who)

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
                     age + 
                     regio + 
                     bcd4_tr + 
                     brna_tr + 
                     #who + 
                     #incident_tb +
                     regimen, 
                   family="binomial", 
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
                      age + 
                      bcd4_tr +
                      incident_tb +
                      #regimen + # does not converge when i include it
                      incident_tb, 
                    family="binomial", 
                    data = data.to.impute_rsa, 
                    control = list(maxit = 50))

summ(miss.mod_rsa)

data.to.impute_rsa$incomplete <- NULL

#### imputation ----------------------------------------------------------------

#' here i include variables that are included in the scientific model 
#' (see complete record analysis) and variables that were significantly related 
#' to missingness in the previous analsysis.

K <- 20 # number of imputed datasets

### CH ###

pmat_ch <- matrix(
  c(0,1,1,1,1,1,1,0, # sex
    1,0,1,1,1,1,1,0, # age
    1,1,0,1,1,1,1,0, # regio
    1,1,1,0,1,1,1,0, # bcd4_tr
    1,1,1,1,0,1,1,0, # brna_tr
    1,1,1,1,1,0,1,0, # who
    1,1,1,1,1,1,0,0, # regimen
    1,1,1,1,1,1,1,0), #incident_tb
  nrow = 8,
  byrow = TRUE)

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
  c(0,1,1,1,1, # sex
    1,0,1,1,1, # age
    1,1,0,1,1, # bcd4_tr
    1,1,1,0,1, # regimen
    1,1,1,1,0), # incident_tb
  nrow = 5, 
  byrow = TRUE)

# shell imputation (for post processing) #

ini_rsa <- mice(data = data.to.impute_rsa, 
                maxit = 0)
post_rsa <- ini_rsa$post
post_rsa["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"

# imputation #

rsa.imp <- mice(data = data.to.impute_rsa,
                m = K,
                method = c("", #sex
                           "", #age
                           "norm",  # bcd4_tr
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
                               sex + age + who + regio + bcd4_tr + brna_tr,
                             family = "binomial"))

fit_ch_imp |>
  pool() |>
  summary(conf.int = TRUE, exponentiate = TRUE) |>
  tibble::column_to_rownames("term") |>
  round(3)

### RSA ###

fit_rsa_imp <- with(data = rsa.imp,
                    exp = glm(incident_tb ~ sex + age + bcd4_tr, 
                              family="binomial"))

fit_rsa_imp |>
  pool() |>
  summary(conf.int = TRUE, exponentiate = TRUE) |>
  tibble::column_to_rownames("term") |>
  round(3)
