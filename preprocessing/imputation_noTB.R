#### libraries -----------------------------------------------------------------

library(tidyverse)
library(mice)
library(tableone)
library(jtools)
library(miceadds)

#### data preparation ---------------------------------------------------------

data.no.tb <- readRDS("data_clean/art_noTB.rds")
custom_breaks <- c(16, 34, 44, 100)

data.no.tb <- data.no.tb %>% 
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
         agegroup = as.factor(agegroup),
         incident_tb = as.numeric(as.character(incident_tb)),
                fup_days = case_when(
                  incident_tb == 1 ~ as.numeric(difftime(date_tb, art_start_date, units = "days")),
                  TRUE ~ fup_time),
                fup_years = fup_days / 360,
                cohort = fct_relevel(cohort, "RSA"),
                agegroup = cut(age, breaks = custom_breaks, labels = c("16-34","35-44", "45+"), include.lowest = TRUE),
         event_type = as.factor(case_when(
             incident_tb == 1 ~ 1,
             death == 1 ~ 2, # only for cases that are incident = 0 (thats how case_when works)
             TRUE ~ 0)),
         sex = fct_drop(sex),
                cd4_group = fct_drop(cd4_group),
                cd4_group = fct_relevel(cd4_group, "350+")) %>%  # Loss to fup is not considered a competing risk) 
  dplyr::select(sex, age, regio, bcd4_tr, bcd4, brna_tr, brna, who, regimen, incident_tb, id, cohort, fup_years, event_type, agegroup) 

data.to.impute.ch <- data.no.tb %>% filter(cohort == "CH") 

data.to.impute.rsa <- data.no.tb %>% filter(cohort == "RSA") 

#### imputation ----------------------------------------------------------------

#### investigate systematic missingness ----------------------------------------

#' syst_miss <- function(df_ti, return) {
#'   
#'   #' the function prints the whole systematic missingness analysis and adds the incomplete column to
#'   #' the respective dataset for the glm
#'   
#'   data.to.impute <- df_ti
#'   
#'   data.to.impute$incomplete <- complete.cases(data.to.impute %>% dplyr::select(-id))
#'   
#'   vars <- names(data.to.impute %>%  select(-id))
#'   
#'   cts <- vars[sapply(data.to.impute, is.numeric)]
#'   
#'   cat <- vars[(!(vars %in% cts))]
#'   
#'   tablemiss <- CreateTableOne(vars = vars, 
#'                               strata="incomplete", 
#'                               data = data.to.impute, 
#'                               factorVars = cat, 
#'                               includeNA = F)
#'   
#'   if (return == "data") {
#'     return(data.to.impute)
#'   } else if (return == "table") {
#'     return(tablemiss)
#'   } else if (return == "md") {
#'     return(md.pattern(data.to.impute, rotate.names = T))
#'   }
#'   else {
#'     return(list(dataframe1 = data.to.impute, dataframe2 = tablemiss))
#'   }
#'   
#' }

data.no.tb$incomplete <- complete.cases(data.no.tb)
data.to.impute.ch$incomplete <- complete.cases(data.to.impute.ch)
data.to.impute.rsa$incomplete <- complete.cases(data.to.impute.rsa)

miss.mod.ch <- glm(incomplete ~ 
                     sex +
                     age + 
                     regio + 
                     bcd4_tr + 
                     brna_tr + 
                     #who + 
                     #incident_tb +
                     regimen, 
                   family= binomial(link = "logit"), 
                   data = data.to.impute.ch,
                   control = list(maxit = 50))

summ(miss.mod.ch, exp = T)

data.to.impute.ch$incomplete <- NULL

### RSA ### 

#data.to.impute.rsa <- syst_miss(data.to.impute.rsa, return = "data")
#syst_miss(data.to.impute.rsa, return = "table")
#syst_miss(data.to.impute.rsa, return = "md") 

miss.mod.rsa <- glm(incomplete ~ 
                      sex +
                      age + 
                      bcd4_tr +
                      incident_tb +
                      regimen + 
                      incident_tb, 
                    family= binomial(link = "logit"), 
                    data = data.to.impute.rsa, 
                    control = list(maxit = 50))

summ(miss.mod.rsa)

data.to.impute.rsa$incomplete <- NULL
data.no.tb$incomplete <- NULL

#### imputation ----------------------------------------------------------------

#' here i include variables that are included in the scientific model 
#' (see complete record analysis) and variables that were significantly related 
#' to missingness in the previous analysis.

K <- 20 # number of imputed datasets

### both ###

ini <- mice(data.no.tb, maxit=0) # dry run without iterations to get the predictor matrix
pred <- ini$predictorMatrix
pred[, c('id', 'fup_years', 'event_type', 'agegroup', 'bcd4', 'brna')] <- 0
post <- ini$post
post["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post["brna_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 10))"

both.imp <- mice(data = data.no.tb,
               m = K,
               method = c("", #sex
                          "", #age
                          "polyreg",  # regio
                          "norm",  # bcd4_tr
                          "~ I(bcd4_tr^2)", #bcd4
                          "norm", #brna_tr 
                          "~ I(10^brna_tr)", #brna
                          "logreg", # who
                          "polyreg", # regimen
                          "", #indicent_tb
                          "", #id
                          "", #cohort
                          "", #fup_years
                          "", #event_type
                          ""), #agegroup
               predictorMatrix = pred,
               seed = 1569,
               maxit = 20, #should be 20
               post = post)

summary(complete(both.imp, 0)) # with missingness, index 0
summary(complete(both.imp, 1)) # imputed index 1:K

saveRDS(both.imp, file = "data_clean/imputed/both_imp.rds")

### CH ###

ini.ch <- mice(data.to.impute.ch, maxit=0) # dry run without iterations to get the predictor matrix
pred.ch <- ini.ch$predictorMatrix
pred.ch[, c('id', 'cohort', 'fup_years', 'event_type', 'agegroup', 'bcd4', 'brna')] <- 0
post.ch <- ini.ch$post
post.ch["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post.ch["brna_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 10))"

ch.imp <- mice(data = data.to.impute.ch,
               m = K,
               method = c("", #sex
                          "", #age
                          "polyreg",  # regio
                          "norm",  # bcd4_tr
                          "~ I(bcd4_tr^2)", #bcd4
                          "norm", #brna_tr 
                          "~ I(10^brna_tr)", #brna
                          "logreg", # who
                          "polyreg", # regimen
                          "", #indicent_tb
                          "", #id
                          "", #cohort
                          "", #fup_years
                          "", #event_type
                          ""), #agegroup
               predictorMatrix = pred.ch,
               seed = 1569,
               maxit = 20, #should be 20
               post = post.ch)

summary(complete(ch.imp, 0)) # with missingness, index 0
summary(complete(ch.imp, 1)) # imputed index 1:K
plot(ch.imp)
saveRDS(ch.imp, file = "data_clean/imputed/ch_imp.rds")


### RSA ###

ini.rsa <- mice(data.to.impute.rsa, maxit=0) # dry run without iterations to get the predictor matrix
pred.rsa <- ini.rsa$predictorMatrix
pred.rsa[, c('id', 'cohort', 'fup_years', 'event_type', 'agegroup', 'bcd4', 'brna')] <- 0
post.rsa <- ini.rsa$post
post.rsa["bcd4_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post.rsa["brna_tr"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 10))"

# imputation #

rsa.imp <- mice(data = data.to.impute.rsa,
                m = K,
                method = c("", #sex
                           "", #age
                           "polyreg",  # regio
                           "norm",  # bcd4_tr
                           "~ I(bcd4_tr^2)", #bcd4
                           "norm", #brna_tr 
                           "~ I(10^brna_tr)", #brna
                           "logreg", # who
                           "polyreg", # regimen
                           "", #indicent_tb
                           "", #id
                           "", #cohort
                           "", #fup_years
                           "", #event_type
                           ""), #agegroup
                predictorMatrix = pred.rsa,
                seed = 1569,
                maxit = 20,
                post = post.rsa)

saveRDS(rsa.imp, file = "data_clean/imputed/rsa_imp.rds")

# check imputations
summary(complete(rsa.imp, 0)) # with missingness, index 0
summary(complete(rsa.imp, 1)) # imputed index 1:K

imp.dta <- mice::complete(both.imp, "all")
imp.dta.long <- mice::complete(both.imp, "long", include = T)
saveRDS(imp.dta, file = "data_clean/imputed/both_imp.complete.rds")
saveRDS(imp.dta.long, file = "data_clean/imputed/both_imp.complete.long.rds")

ch.imp.dta <- mice::complete(ch.imp, "all")
ch.imp.dta.long <- mice::complete(ch.imp, "long", include = T)
saveRDS(ch.imp.dta, file = "data_clean/imputed/ch_imp.complete.rds")
saveRDS(ch.imp.dta.long, file = "data_clean/imputed/ch_imp.complete.long.rds")

rsa.imp.dta <- mice::complete(rsa.imp, "all")
rsa.imp.dta.long <- mice::complete(rsa.imp, "long", include = T)
saveRDS(rsa.imp.dta, file = "data_clean/imputed/rsa_imp.complete.rds")
saveRDS(rsa.imp.dta.long, file = "data_clean/imputed/rsa_imp.complete.long.rds")

