##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  lmerTest,
  mgcv,
  gridExtra
)

source("utils/plot.R")

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
  filter(!is.na(cd4)) %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1),
    trans.cd4 = sqrt(cd4)
)

#### data exploration ####

summary_stats.cd4 <- rec.cd4 %>% 
  group_by(group_30days) %>% 
  summarise(
    median_trans.cd4 = median(trans.cd4, na.rm = TRUE),
    Q1 = quantile(trans.cd4, 0.25, na.rm = TRUE),
    Q3 = quantile(trans.cd4, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(trans.cd4, 0.05, na.rm = TRUE),
    Qhigh = quantile(trans.cd4, 0.95, na.rm = TRUE),
    mean = mean(trans.cd4, na.rm = TRUE),
    se = sd(trans.cd4, na.rm = TRUE) / sqrt(n()),
    sd = sd(trans.cd4, na.rm = TRUE)
  )

IQR.cd4 <- ggplot(summary_stats.cd4, aes(x = time_diff_midpoint, y = median_trans.cd4)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw() +
  geom_hline(yintercept = sqrt(350)) 

IQR.cd4

se.cd4 <- ggplot(summary_stats.cd4, aes(x = time_diff_midpoint, y = mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean - sd, 
        ymax = mean + sd),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw()+   
  geom_hline(yintercept = sqrt(350)) 
  
rec.rna <- rec %>% 
  filter(!is.na(rna)) %>% 
  mutate(
    trans.rna = case_when(
      rna == 0 ~ log10(runif(1, min = 1, max = 50)),
      TRUE ~ log10(rna)
    )
  )

## checking distribution of measurements
rec.rna <- rec.rna %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1)
  )

num_unique_ids <- rec.rna %>% 
  summarise(n_distinct(id))

num_unique_ids2 <- rec.cd4 %>% 
  summarise(n_distinct(id))

# Calculate Q1 and Q3 for each group
summary_stats.rna <- rec.rna %>% 
  group_by(group_30days) %>% 
  summarise(
    median_trans.rna = median(trans.rna, na.rm = TRUE),
    Q1 = quantile(trans.rna, 0.25, na.rm = TRUE),
    Q3 = quantile(trans.rna, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(trans.rna, 0.05, na.rm = TRUE),
    Qhigh = quantile(trans.rna, 0.95, na.rm = TRUE),
    Mean = mean(trans.rna, na.rm = TRUE),
    se = sd(trans.rna, na.rm = TRUE) / sqrt(n()),
    sd = sd(trans.rna, na.rm = TRUE)
  )

# Plot the data with IQR error bars
IQR.rna <- ggplot(summary_stats.rna, aes(x = time_diff_midpoint, y = median_trans.rna)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

se.rna <- ggplot(summary_stats.rna, aes(x = time_diff_midpoint, y = Mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Mean - sd, 
        ymax = Mean + sd),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

IQR.rna

ggplot(rec.rna, aes(x = time_diff)) +
  geom_histogram(binwidth = 7) +
  theme_bw() +
  labs(x = "Time Difference", y = "Frequency")

distance <- grid.arrange(IQR.cd4, IQR.rna, ncol = 2)

se <- grid.arrange(se.cd4, se.rna, ncol = 2)

ggsave(plot = distance, file = "results/slopes/IQR.png")

#### assumptions ####
ggplot(rec.cd4, aes(x= sqrt(cd4))) +
  geom_histogram() +
  theme_bw()

ggplot(rec.rna, aes(x=log10(rna+1))) +
  geom_histogram(binwidth = 0.5) +
  theme_bw()

# Q-Q plot for cd4
ggplot(rec, aes(sample = sqrt(cd4))) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Q-Q plot for cd4")

#### raw plot CD4 ####

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

#### raw plot HIV-RNA ####

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

#### model ####
## CD4 
m.cd4 <- gam(sqrt(cd4) ~ s(time_diff, k = 5) + s(id, bs="re"), data = rec.cd4, method = "REML")

plot.gam(m.cd4)

summary(m.cd4)

residuals.cd4 = resid(m.cd4)
hist(residuals.cd4, main = "Histogram of Residuals", xlab = "Residuals")

ggplot(rec.cd4, aes(x = time_diff, y = residuals.cd4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

## RNA 

m.rna <- gam(trans.rna ~ s(time_diff, k = 7) + s(id, bs="re"), data = rec.rna, method = "REML")

gam.check(m.rna)

summary(m.rna)

plot.gam(m.rna, 
         seWithMean = TRUE,
         unconditional=TRUE)

ggplot(rec.rna, aes(x = time_diff)) +
  geom_histogram(binwidth = 14) +
  theme_bw() +
  labs(x = "Time Difference", y = "Frequency")

residuals.rna = resid(m.rna)
hist(residuals.rna, main = "Histogram of Residuals", xlab = "Residuals")

rec.rna$residuals_rna <- resid(m.rna)

ggplot(rec.rna, aes(x = time_diff, y = residuals_rna)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

#### Plotting ----

# I calculate the simultaneous CI instead of point-wise https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

### Fit only

#CD4

trend.cd4 <- plot_trend(m.cd4, rec.cd4) +
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  labs(x = "", y = expression(sqrt(CD4))) +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 30))

trend.cd4

#HIV-RNA

trend.rna <- plot_trend(m.rna, rec.cd4)+
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  labs(x = "", y = "log10(viral-load)") +
  coord_cartesian(xlim = c(-30, 360), ylim = c(0, 7.5)) 

trend.rna

trends <- grid.arrange(trend.cd4, trend.rna, ncol = 2)

ggsave(plot = trends, filename = "results/slopes/trends.jpeg")

### raw + fitted

trend.cd4.raw <- plot_trend.rawcd4(m.cd4, rec.cd4) 

trend.rna.raw <- plot_trend.raw.rna(m.rna, rec.rna)

trends.raw <- grid.arrange(trend.cd4.raw, trend.rna.raw, ncol = 2)

ggsave(plot = trends.raw, filename = "results/slopes/trends.raw.png")

### all 

trends_complete <- grid.arrange(trend.cd4.raw, trend.cd4, trend.rna.raw, trend.rna, ncol =2)
ggsave(plot = trends_complete, filename = "results/slopes/trends.both.png")

trends.raw

#### Prediction interval #### 
#https://mikl.dk/post/2019-prediction-intervals-for-gam/

## CD4

prediction_cd4 <- plot_pred.cd4(m.cd4, rec.cd4)

prediction_cd4

prediction_rna <- plot_pred.rna(m.rna, rec.rna)

prediction_rna

##
conf.pred <- grid.arrange(trend.cd4.raw, prediction_cd4, trend.rna.raw, prediction_rna, ncol = 2)

ggsave(plot = conf.pred, filename = "results/slopes/confvspred.png")
