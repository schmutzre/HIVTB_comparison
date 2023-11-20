library(tidyr)
library(tidyverse)

rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- chol(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

plot_trend <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, expand.grid(time_diff = seq(min(time_diff), max(time_diff), length = 200),
                                 cohort = levels(cohort),
                                 KEEP.OUT.ATTRS = FALSE))
  
  # Compute predicted values and standard errors
  pred <- predict(model, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  
  # Simulate the conditional distribution of the model's random effects
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  # Compute linear predictor matrix
  Cg <- predict(model, newd, type = "lpmatrix")
  
  # Simulate deviance residuals
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  
  # Compute the maximum absolute scaled deviance (MASD) for each simulation
  masd <- apply(absDev, 2L, max)
  
  # Compute the 95th percentile of the MASD distribution
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  # Add the computed intervals to the prediction data frame
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit)) %>% 
    mutate(cohort = fct_relevel(cohort, "RSA"))
  
  pred$lwrS_adj <- pmax(pred$lwrS, 0) # lower end of Confidence band shouldnt go below 0 (not possible)
  
  # Plot the trend with ggplot2
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
    geom_line(aes(y = fit, color = cohort)) +
    theme_classic()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) +
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2")) +
    coord_cartesian(ylim = c(0,800))
  
  return(trend.plot)
}

plot_trend.raw.rna <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, expand.grid(time_diff = seq(min(time_diff), max(time_diff), length = 200),
                                 cohort = levels(cohort),
                                 KEEP.OUT.ATTRS = FALSE))
  
  # Compute predicted values and standard errors
  pred <- predict(model, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  
  # Simulate the conditional distribution of the model's random effects
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  # Compute linear predictor matrix
  Cg <- predict(model, newd, type = "lpmatrix")
  
  # Simulate deviance residuals
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  
  # Compute the maximum absolute scaled deviance (MASD) for each simulation
  masd <- apply(absDev, 2L, max)
  
  # Compute the 95th percentile of the MASD distribution
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  # Add the computed intervals to the prediction data frame
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit)) %>% 
    mutate(cohort = fct_relevel(cohort, "RSA"))
  
  pred$lwrS_adj <- pmax(pred$lwrS, 0) # lower end of Confidence band shouldnt go below 0 (not possible)
  
  # Plot the trend with ggplot2
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_line(data = data, aes(group = factor(id), y = log10(rna+1)), color = "lightgrey", alpha = .2) +
    geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS, fill = cohort), alpha = 0.2) +
    geom_line(aes(y = fit, color = cohort)) +
    theme_classic()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = log10(400)) +
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2"))
  
  
  return(trend.plot)
}


plot_trend.raw.cd4 <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, expand.grid(time_diff = seq(min(time_diff), max(time_diff), length = 200),
                                 cohort = levels(cohort),
                                 KEEP.OUT.ATTRS = FALSE))
  
  # Compute predicted values and standard errors
  pred <- predict(model, newd, se.fit = TRUE)
  se.fit <- pred$se.fit
  
  # Simulate the conditional distribution of the model's random effects
  BUdiff <- rmvn(N, mu = rep(0, nrow(Vb)), sig = Vb)
  
  # Compute linear predictor matrix
  Cg <- predict(model, newd, type = "lpmatrix")
  
  # Simulate deviance residuals
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  
  # Compute the maximum absolute scaled deviance (MASD) for each simulation
  masd <- apply(absDev, 2L, max)
  
  # Compute the 95th percentile of the MASD distribution
  crit <- quantile(masd, prob = 0.95, type = 8)
  
  # Add the computed intervals to the prediction data frame
  pred <- transform(cbind(data.frame(pred), newd),
                    uprP = fit + (2 * se.fit),
                    lwrP = fit - (2 * se.fit),
                    uprS = fit + (crit * se.fit),
                    lwrS = fit - (crit * se.fit)) %>% 
    mutate(cohort = fct_relevel(cohort, "RSA"))
  
  # Plot the trend with ggplot2
plot_name <- deparse(substitute(model))
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_line(data = data, aes(group = factor(id), y = cd4), color = "lightgrey", alpha = .1) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
    geom_line(aes(y = fit, color = cohort)) +
    theme_classic()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = 350) +
    geom_vline(xintercept = 0, linetype = "dotted")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none") +
    scale_color_manual(values = wes_palette("Moonrise2")) +
    scale_fill_manual(values = wes_palette("Moonrise2")) 
  
  return(trend.plot)
}

plot_pred.cd4 <- function(model, rec_data, i) {
  beta <- coef(model)
  V <- vcov(model)
  
  # simulate parameters
  n_sim <- 10000
  Cv <- chol(V)
  set.seed(1)
  nus <- rnorm(n_sim * length(beta))
  beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = n_sim)
  
  # make linear predictions
  X <- expand.grid(time_diff = seq(min(rec_data$time_diff), max(rec_data$time_diff), 10), id = 1)
  covar_sim <- predict(model, newdata = X, type = "lpmatrix", exclude = "s(id)") # exclude random effects
  linpred_sim <- covar_sim %*% beta_sims
  
  # simulate predictions (sample from normal distribution)
  pred <- matrix(rnorm(n = prod(dim(linpred_sim)), mean = linpred_sim, sd = sqrt(summary(model)$scale)), 
                 nrow = nrow(linpred_sim), ncol = ncol(linpred_sim))
  
  # create data frame
  pred_df <- t(pred) %>%
    data.frame() %>%
    set_names(sort(X$time_diff)) %>%
    gather() %>%
    rename(time_diff = key, 
           trans_pred = value) %>%
    mutate(time_diff = as.numeric(gsub("X", "", time_diff))) %>%
    group_by(time_diff) %>%
    mutate(draw = 1:n()) %>%
    ungroup()
  
  # summarize
  pred_df_sum <- pred_df %>%
    group_by(time_diff) %>%
    summarize(mean_pred = mean(trans_pred),
              lower_pred = quantile(trans_pred, .025),
              upper_pred = quantile(trans_pred, .975)) %>%
    ungroup() 
  
  # plot
  p <- ggplot(mapping = aes(time_diff)) +
    geom_line(data = rec_data, mapping = aes(y = cd4_trans, group = factor(id)), color = "grey", alpha = .2) +
    geom_line(data = pred_df_sum, mapping = aes(y = mean_pred), color = wes_palette("Moonrise2")[i], linewidth = 1.5) + 
    geom_line(data = pred_df_sum, mapping = aes(y = lower_pred), color = wes_palette("Moonrise2")[i], linetype = "dashed") +
    geom_line(data = pred_df_sum, mapping = aes(y = upper_pred), color = wes_palette("Moonrise2")[i], linetype = "dashed") +
    geom_vline(mapping = aes(xintercept = 0), linetype = "dotted") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x = "Days since ART", y = expression(sqrt(CD4))) +
    geom_hline(mapping = aes(yintercept = sqrt(350)), color = "black")  
  
  return(p)
}

# To use the function:
# p <- plot_pred.int(m.cd4, rec.cd4)
# print(p)

plot_pred.rna <- function(model, rec_data) {
  beta <- coef(model)
  V <- vcov(model)
  
  # simulate parameters
  n_sim <- 10000
  Cv <- chol(V)
  set.seed(1)
  nus <- rnorm(n_sim * length(beta))
  beta_sims <- beta + t(Cv) %*% matrix(nus, nrow = length(beta), ncol = n_sim)
  
  # make linear predictions
  X <- expand.grid(time_diff = seq(min(rec_data$time_diff), max(rec_data$time_diff), 10), id = 1)
  covar_sim <- predict(model, newdata = X, type = "lpmatrix", exclude = "s(id)") # exclude random effects
  linpred_sim <- covar_sim %*% beta_sims
  
  # simulate predictions (sample from normal distribution)
  pred <- matrix(rnorm(n = prod(dim(linpred_sim)), mean = linpred_sim, sd = sqrt(summary(model)$scale)), 
                 nrow = nrow(linpred_sim), ncol = ncol(linpred_sim))
  
  # create data frame
  pred_df <- t(pred) %>%
    data.frame() %>%
    set_names(sort(X$time_diff)) %>%
    gather() %>%
    rename(time_diff = key, 
           trans_pred = value) %>%
    mutate(time_diff = as.numeric(gsub("X", "", time_diff))) %>%
    group_by(time_diff) %>%
    mutate(draw = 1:n()) %>%
    ungroup()
  
  # summarize
  pred_df_sum <- pred_df %>%
    group_by(time_diff) %>%
    summarize(mean_pred = mean(trans_pred),
              lower_pred = max(0,quantile(trans_pred, .025)),
              upper_pred = quantile(trans_pred, .975)) %>%
    ungroup() 
  
  # plot
  p <- ggplot(mapping = aes(time_diff)) +
    geom_line(data = rec_data, mapping = aes(y = trans.rna, group = factor(id)), color = "grey", alpha = .2) +
    geom_line(data = pred_df_sum, mapping = aes(y = mean_pred), color = "blue") + 
    geom_line(data = pred_df_sum, mapping = aes(y = lower_pred), color = "blue", linetype = "dashed") +
    geom_line(data = pred_df_sum, mapping = aes(y = upper_pred), color = "blue", linetype = "dashed") +
    geom_vline(mapping = aes(xintercept = 0), linetype = "dotted") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x = "Days since ART", y = "log10(viral-load)") +
    geom_hline(mapping = aes(yintercept = log10(400)), color = "black")  
  
  return(p)
}



