rmvn <- function(n, mu, sig) { ## MVN random deviates
  L <- chol(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}

plot_trend <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, data.frame(time_diff = seq(min(time_diff), max(time_diff), length = 200), id = 1))
  
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
                    lwrS = fit - (crit * se.fit))
  
  # Plot the trend with ggplot2
  plot_name <- deparse(substitute(model))
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red") +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "blue") +
    geom_line(aes(y = fit), color = "red") +
    theme_bw()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(trend.plot)
}

plot_trend.raw.rna <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, data.frame(time_diff = seq(min(time_diff), max(time_diff), length = 200), id = 1))
  
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
                    lwrS = fit - (crit * se.fit))
  
  pred$lwrS_adj <- pmax(pred$lwrS, 0) # lower end of Confidence band shouldnt go below 0 (not possible)
  
  # Plot the trend with ggplot2
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_line(data = data, aes(group = factor(id), y = log10(rna+1)), color = "grey", alpha = .1) +
    geom_ribbon(aes(ymin = lwrS_adj, ymax = uprS), alpha = 0.2, fill = "red") +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "blue") +
    geom_line(aes(y = fit), color = "red") +
    theme_bw()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = log10(400)) +
    geom_vline(xintercept = 0, linetype = "dashed")+
    labs(x = "", y = "log10(viral-load)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  
  return(trend.plot)
}


plot_trend.rawcd4 <- function(model, data, N = 10000) {
  
  # Compute Vb, the covariance matrix of the random effects
  Vb <- vcov(model, unconditional = TRUE)
  
  # Generate new data for predictions
  newd <- with(data, data.frame(time_diff = seq(min(time_diff), max(time_diff), length = 200), id = 1))
  
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
                    lwrS = fit - (crit * se.fit))
  
  # Plot the trend with ggplot2
plot_name <- deparse(substitute(model))
  trend.plot <- ggplot(pred, aes(x = time_diff)) +
    geom_line(data = data, aes(group = factor(id), y = sqrt(cd4)), color = "grey", alpha = .1) +
    geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2, fill = "red") +
    geom_ribbon(aes(ymin = lwrP, ymax = uprP), alpha = 0.2, fill = "blue") +
    geom_line(aes(y = fit), color = "red") +
    theme_bw()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0)) +
    geom_hline(yintercept = sqrt(350)) +
    geom_vline(xintercept = 0, linetype = "dashed")+
    labs(x = "", y = "sqrt(CD4 count)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(trend.plot)
}




