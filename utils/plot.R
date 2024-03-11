library(tidyr)
library(tidyverse)

#### Plotting CD4 immun recovery (interaction) ---------------------------------

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- chol(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

pred_trend <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model$gam)
  
  # Generate new data for predictions
  pred_data <- with(data, expand.grid(time_trans = seq(min(time_trans), max(time_trans), length = 200),
                                     cohort = levels(cohort),
                                     presenting_tb = levels(presenting_tb))) %>% 
    mutate(interaction_term = interaction(cohort, presenting_tb))
  
  prediction_cd4 <- predict(model$gam, type = "response", se.fit = TRUE, newdata = pred_data)
  fit <- prediction_cd4$fit
  se.fit <-  prediction_cd4$se.fit
  fit_org <- fit^2
  se.fit_org <- abs(2 * fit) * se.fit #delta method 
  
  N <- 10000
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  Cg <- predict(model$gam, pred_data, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  pred <- transform(cbind(data.frame(prediction_cd4), pred_data),
                    fit_org = fit_org,
                    uprP = fit_org + (2 * se.fit_org),
                    lwrP = fit_org - (2 * se.fit_org),
                    uprS = fit_org + (crit * se.fit_org),
                    lwrS = fit_org - (crit * se.fit_org)) %>% 
    mutate(presenting_tb = case_when(presenting_tb == 0 ~ "Without prevalent TB",
                                     TRUE ~ "With prevalent TB"),
           cohort = fct_relevel(cohort, "RSA"),
           time_org = time_trans - 60,
           lwrS = ifelse(lwrS < 0, 0, lwrS))

  return(pred)
  
}

#### Plotting CD4 immun recovery (seperate) ------------------------------------

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- chol(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

pred_trend_uni <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model$gam)
  
  # Generate new data for predictions
  pred_data <- with(data, expand.grid(time_diff = seq(min(time_diff), max(time_diff), length = 200),
                                      cohort = levels(cohort)))
  
  prediction_cd4 <- predict(model$gam, type = "response", se.fit = TRUE, newdata = pred_data)
  fit <- prediction_cd4$fit
  se.fit <-  prediction_cd4$se.fit
  fit_org <- fit^2
  se.fit_org <- abs(2 * fit) * se.fit #delta method
  
  N <- 10000
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  Cg <- predict(model$gam, pred_data, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  pred <- transform(cbind(data.frame(prediction_cd4), pred_data),
                    fit_org = fit_org,
                    uprP = fit_org  + (1.96 * se.fit_org),
                    lwrP = fit_org  - (1.96 * se.fit_org),
                    uprS = fit_org  + (crit * se.fit_org),
                    lwrS = fit_org  - (crit * se.fit_org)) %>% 
    mutate(cohort = fct_relevel(cohort, "RSA"),
           lwrS = if_else(lwrS < 0, 0, lwrS))
  
  return(pred)
  
}
