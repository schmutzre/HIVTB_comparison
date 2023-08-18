#### Library ####

library(tidyverse)
library(mgcv)

source("utils/plot.R")
  
#### Data ####

rec_long_ch <- readRDS("data_clean/art_ch.long.rds") 

rec <- rec_long_ch %>% 
  filter(time_diff > -30 & time_diff < 365)

rec.cd4 <- rec %>% 
  filter(!is.na(cd4)) %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1),
    trans.cd4 = sqrt(cd4)
  )

rec.rna <- rec %>% 
  filter(!is.na(rna)) %>% 
  mutate(
    trans.rna = case_when(
      rna == 0 ~ log10(runif(1, min = 1, max = 50)),
      TRUE ~ log10(rna)
    )
  )

#### Models ####

m.cd4 <- gam(trans.cd4 ~ s(time_diff, k = 5) + s(id, bs="re"), data = rec.cd4, method = "REML")
m.rna <- gam(trans.rna ~ s(time_diff, k = 7) + s(id, bs="re"), data = rec.rna, method = "REML")

#### Plots  ####

trend.cd4 <- plot_trend(m.cd4, rec.cd4) +
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "", y = "sqrt(CD4 count)") +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 30))

trend.rna <- plot_trend(m.rna, rec.cd4)+
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  labs(x = "", y = "log10(viral-load)") +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 7.5)) 

print(trend.cd4)
print(trend.rna)

#### Prediction intervals ####

# Replace these lines with your own data
# Example:
x <- rec.cd4$time_diff
y <- rec.cd4$trans.cd4

beta <- coef(m.cd4)
Vb <- vcov(m.cd4)

## simulate replicate beta vectors from posterior...
Cv <- chol(Vb)
n.rep=100
nb <- length(beta)
br <- t(Cv) %*% matrix(rnorm(n.rep*nb),nb,n.rep) + beta

## turn these into replicate linear predictors...
xp <- seq(from = -30, to = 360, length.out = 1000)
Xp <- predict(m.cd4,newdata=data.frame(time_diff=xp, id = 1),type="lpmatrix")
fv <- Xp%*%br ## ... finally, replicate expected value vectors

## and estimated scale...

plot(rep(xp,n.rep),fv,pch=".", ylim = c(0,40)) ## plotting replicates
points(x,y,pch=19,cex=.5) ## and original data

## compute 95% prediction interval...
PI <- apply(fv,1,quantile,prob=c(.025,0.975))
## and plot it...
lines(xp,PI[1,],col=2,lwd=2);lines(xp,PI[2,],col=2,lwd=2)

## Confidence interval for comparison...
pred <- predict(m.cd4,newdata=data.frame(time_diff=xp, id = 1),se=TRUE)
lines(xp,pred$fit,col=3,lwd=2)
u.ci <- pred$fit + 2*pred$se.fit
l.ci <- pred$fit - 2*pred$se.fit
lines(xp,u.ci,col=3,lwd=2);lines(xp,l.ci,col=3,lwd=2)

