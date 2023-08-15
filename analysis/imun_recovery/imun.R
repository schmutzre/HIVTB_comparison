##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  lmerTest,
  mgcv
)

# Set seed for reproducibility
set.seed(123)

#### data import ----

rec_long_ch <- readRDS("data_clean/art_ch.long.rds") %>% 
  filter(!(is.na(cd4) & is.na(rna)))

rec_sa <- ##

# Join this dataframe with the original dataframe

rec <- rec_long_ch %>% 
  filter(time_diff > -30 & time_diff < 365)

rec.cd4 <- rec %>% 
  filter(!is.na(cd4))

rec.rna <- rec %>% 
  filter(!is.na(rna))

#### assumptions ----
ggplot(rec.cd4, aes(x= sqrt(cd4))) +
  geom_histogram() +
  theme_bw()

ggplot(rec.rna, aes(x=log10(rna+1))) +
  geom_histogram() +
  theme_bw()

# Q-Q plot for cd4
ggplot(rec, aes(sample = sqrt(cd4))) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Q-Q plot for cd4")

#### raw plot CD4 ----

cd4_yearRAW <- rec.cd4 %>%
  ggplot(aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = sqrt(cd4)), color = "grey", alpha = .2) +
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "", title = "raw data + GAM") +
  scale_x_continuous(limits = c(-30,360)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,50)) +
  theme(plot.title = element_text(hjust = 0.5))

print(cd4_yearRAW)

#### raw plot HIV-RNA ----

rna_yearRAW <- rec.rna %>%
  ggplot(aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = log10(rna+1)), color = "grey", alpha = .2) +
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "Days since ART start", title = "") +
  scale_x_continuous(limits = c(-30,360)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,10)) +
  theme(plot.title = element_text(hjust = 0.5))

print(rna_yearRAW)

#### spline model ----
## CD4 
model.cd4_spline <- gam(sqrt(cd4) ~ s(time_diff) + s(id, bs="re"), data = rec.cd4, method = "REML")

summary(model.cd4_spline)

residuals.cd4 = resid(model.cd4_spline)
hist(residuals.cd4, main = "Histogram of Residuals", xlab = "Residuals")

## RNA 

model.rna_spline <- gam(log10(rna+1) ~ s(time_diff) + s(id, bs="re"), data = rec.rna, method = "REML")

summary(model.rna_spline)

residuals.rna = resid(model.rna_spline)
hist(residuals.rna, main = "Histogram of Residuals", xlab = "Residuals")

### Plotting ----

#CD4
newdata.cd4 <- data.frame(
  time_diff = seq(min(rec.cd4$time_diff), max(rec.cd4$time_diff), length.out=100),
  id = rep(rec.cd4$id)[1], 100)

preds.cd4 <- predict(model.cd4_spline, newdata=newdata.cd4, se=TRUE)
newdata.cd4$fit <- preds$fit
newdata.cd4$lowerCI <- preds.cd4$fit - 1.96 * preds.cd4$se.fit
newdata.cd4$upperCI <- preds.cd4$fit + 1.96 * preds.cd4$se.fit

trend.cd4 <- ggplot(rec.cd4, aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = sqrt(cd4)), color = "grey", alpha = .1) +
  geom_line(data=newdata.cd4, aes(x = time_diff, y=fit), color="red") +
  geom_ribbon(data=newdata.cd4, aes(x = time_diff, ymin=lowerCI, ymax=upperCI), alpha=0.2, fill="red")+
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "", title = "raw data + GAM") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 50)) +
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(hjust = 0.5))

trend.cd42 <- ggplot(newdata.cd4, aes(x = time_diff)) +
  geom_line(aes(y=fit), color="red") +
  geom_ribbon(aes(x = time_diff, ymin=lowerCI, ymax=upperCI), alpha=0.2, fill="red")+
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "", title = "GAM", y = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-30, 360), ylim = c(15, 25)) +
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(hjust = 0.5))

print(trend.cd42)

## RNA

plot(model.rna_spline, se=TRUE)

newdata.rna <- data.frame(
  time_diff = seq(min(rec.rna$time_diff), max(rec.rna$time_diff), length.out=100),
  id = rep(rec.rna$id)[1], 100)

preds.rna <- predict(model.rna_spline, newdata=newdata.rna, se=TRUE)
newdata.rna$fit <- preds.rna$fit
newdata.rna$lowerCI <- preds.rna$fit - 1.96 * preds.rna$se.fit
newdata.rna$upperCI <- preds.rna$fit + 1.96 * preds.rna$se.fit

trend.rna <- ggplot(rec.rna, aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = log10(rna)), color = "grey", alpha = .1) +
  geom_line(data=newdata.rna, aes(x = time_diff, y=fit), color="red") +
  geom_ribbon(data=newdata.rna, aes(x = time_diff, ymin=lowerCI, ymax=upperCI), alpha=0.2, fill="red") +
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "", title = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 10)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.title = element_text(hjust = 0.5))

trend.rna2 <- ggplot(newdata.rna, aes(x = time_diff)) +
  geom_line(aes(y=fit), color="red") +
  geom_ribbon(aes(x = time_diff, ymin=lowerCI, ymax=upperCI), alpha=0.2, fill="red")+
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "", title = "", y = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 6)) +
  scale_y_continuous(expand = c(0,0))+
  theme(plot.title = element_text(hjust = 0.5))

print(trend.rna2)

trends <- grid.arrange(trend.cd4, trend.cd42, trend.rna, trend.rna2, ncol = 2)

ggsave(plot = trends, filename = "results/slopes.png")
