##### set-up -------------------------------------------------------------------

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

#### data ---------------------------------------------------------------

rna_ch <- readRDS("data_clean/ch/rna_ch.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"))
cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"))
rna_rsa <- readRDS("data_clean/rsa/rna_rsa.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))
cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"))

cd4 <- rbind(cd4_ch %>% select(-timepoint),
             cd4_rsa %>% select(-timepoint)) %>% 
  filter(time_diff > -60 & time_diff < 365) %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1),
    cd4_trans = sqrt(cd4))

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

rna <- rbind(rna_ch %>% select(-timepoint),
             rna_rsa %>% select(-timepoint)) %>% 
  filter(time_diff > -60 & time_diff < 365) %>% 
  mutate(group_30days = as.integer((time_diff %/% 30) + 1),
         rna_trans = log10(rna + 1))

rna_npres <- rna %>% 
  filter(presenting_tb == 0)

rna_pres <- rna %>% 
  filter(presenting_tb == 1)

#### models --------------------------------------------------------------------

### cd4 ###

m.cd4_npres <- gam(cd4 ~ s(time_diff, k = 4, by = cohort), 
             data = cd4_npres, 
             method = "REML")

m.cd4_pres <- gam(cd4 ~ s(time_diff, k = 4, by = cohort), 
                   data = cd4_pres, 
                   method = "REML")

### rna ###

m.rna_npres <- gam(rna_trans ~ s(time_diff, k = 4, by = cohort), 
             data = rna_npres, 
             method = "REML")

m.rna_pres <- gam(rna_trans ~ s(time_diff, k = 4, by = cohort), 
             data = rna_pres, 
             method = "REML")

#### data exploration ----------------------------------------------------------

# cd4 #

summary_stats_cd4 <- cd4 %>% 
  group_by(group_30days) %>% 
  summarise(
    md_trans_cd4 = median(cd4_trans, na.rm = TRUE),
    Q1 = quantile(cd4_trans, 0.25, na.rm = TRUE),
    Q3 = quantile(cd4_trans, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(cd4_trans, 0.05, na.rm = TRUE),
    Qhigh = quantile(cd4_trans, 0.95, na.rm = TRUE),
    mean = mean(cd4_trans, na.rm = TRUE),
    se = sd(cd4_trans, na.rm = TRUE) / sqrt(n()),
    sd = sd(cd4_trans, na.rm = TRUE)
  )

iqr_cd4 <- ggplot(summary_stats_cd4, aes(x = time_diff_midpoint, y = md_trans_cd4)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw() +
  geom_hline(yintercept = sqrt(350)) 

iqr_cd4

sd_cd4 <- ggplot(summary_stats_cd4, aes(x = time_diff_midpoint, y = mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = mean - sd, 
        ymax = mean + sd),
    width = 10
  ) +
  labs(x = 'Time Difference (days)', y = 'Square root of CD4') +
  theme_bw()+   
  geom_hline(yintercept = sqrt(350)) 
  
sd_cd4

# rna #

summary_stats_rna <- rna %>% 
  group_by(group_30days) %>% 
  summarise(
    md_trans_rna = median(rna_trans, na.rm = TRUE),
    Q1 = quantile(rna_trans, 0.25, na.rm = TRUE),
    Q3 = quantile(rna_trans, 0.75, na.rm = TRUE),
    time_diff_midpoint = mean(time_diff, na.rm = TRUE),
    Qlow = quantile(rna_trans, 0.05, na.rm = TRUE),
    Qhigh = quantile(rna_trans, 0.95, na.rm = TRUE),
    Mean = mean(rna_trans, na.rm = TRUE),
    se = sd(rna_trans, na.rm = TRUE) / sqrt(n()),
    sd = sd(rna_trans, na.rm = TRUE))

iqr_rna <- ggplot(summary_stats_rna, aes(x = time_diff_midpoint, y = md_trans_rna)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Q1, ymax = Q3),
    width = 10) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

sd_rna <- ggplot(summary_stats_rna, aes(x = time_diff_midpoint, y = Mean)) +
  geom_point() +
  geom_errorbar(
    aes(ymin = Mean - sd, 
        ymax = Mean + sd),
    width = 10) +
  labs(x = 'Time Difference (days)', y = 'Log-transformed RNA') +
  theme_bw() +
  geom_hline(yintercept = log10(400)) 

iqr_rna

plot_se <- grid.arrange(sd_cd4, sd_rna, ncol = 2)

ggsave(plot = plot_se, file = "results/slopes/se.png")

#### model check ---------------------------------------------------------------

# cd4 #

plot.gam(m.cd4)
summary(m.cd4)
residuals.cd4 = resid(m.cd4)
hist(residuals.cd4, main = "Histogram of Residuals", xlab = "Residuals")

ggplot(rec.cd4, aes(x = time_diff, y = residuals.cd4)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

# rna #

plot.gam(m.rna)
summary(m.rna)
residuals.rna = resid(m.rna)
hist(residuals.rna, main = "Histogram of Residuals", xlab = "Residuals")

ggplot(rec.rna, aes(x = time_diff, y = residuals.rna)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  theme_bw() +
  labs(x = "Time Difference", y = "Residuals")

#### plots ---------------------------------------------------------------------

#### Confidence intervals + raw  ----

# I calculate the simultaneous CI instead of point-wise https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

## cd4 ##

## not presenting

trend_cd4_npres <- plot_trend.raw.cd4(m.cd4_npres, cd4_npres) +
  facet_wrap(~cohort) +
  ylab("not presenting with TB")+
  coord_cartesian(ylim = c(0,2000))
  
trend_cd4_npres

## presenting

trend_cd4_pres <- plot_trend.raw.cd4(m.cd4_pres, cd4_pres) +
  facet_wrap(~cohort) +
  theme(strip.text.x = element_blank(), strip.background = element_blank())+
  ylab("presenting with TB") +
  coord_cartesian(ylim = c(0,2000))

trend_cd4_pres

## combine

trend_cd4 <- ggarrange(trend_cd4_npres, 
                             trend_cd4_pres, 
                             ncol = 1) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "cd4")

trend_cd4

ggsave(plot = trend_cd4, filename = "results/slopes/cd4.png", 
       width = 16, height = 11, units = "cm")

# ------------------------- 

## rna ##

## not presenting

trend_rna_npres <- plot_trend.raw.rna(m.rna_npres, rna_npres) +
  facet_wrap(~cohort) +
  ylab("not presenting with TB")

trend_rna_npres

## presenting

trend_rna_pres <- plot_trend.raw.rna(m.rna_pres, rna_pres) +
  facet_wrap(~cohort) +
  theme(strip.text.x = element_blank(), strip.background = element_blank())+
  ylab("presenting with TB")

trend_rna_pres

## combine

trend_rna <- ggarrange(trend_rna_npres, 
                       trend_rna_pres, 
                       ncol = 1) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = "log10(rna)")

trend_rna

ggsave(plot = trend_rna, filename = "results/slopes/rna.png", 
       width = 16, height = 11, units = "cm")

#### prediction interval ----

### cd4

## CH non-presenting ##

cd4_npres_ch <- cd4_npres %>% 
  filter(cohort == "CH")

m.cd4_npres_ch <- gam(cd4_trans ~ s(time_diff, k = 3), 
             data = cd4_npres_ch, 
             method = "REML")

pred_cd4_npres_ch <- plot_pred.cd4(m.cd4_npres_ch, cd4_npres_ch, 2) 
pred_cd4_npres_ch

## CH presenting ##

cd4_pres_ch <- cd4_pres %>% 
  filter(cohort == "CH")

m.cd4_pres_ch <- gam(cd4_trans ~ s(time_diff, k = 3), 
                      data = cd4_pres_ch, 
                      method = "REML")

pred_cd4_pres_ch <- plot_pred.cd4(m.cd4_pres_ch, cd4_pres_ch, 2) 
pred_cd4_pres_ch

## RSA non-presenting ##

cd4_npres_rsa <- cd4_npres %>% 
  filter(cohort == "RSA")

m.cd4_npres_rsa <- gam(cd4_trans ~ s(time_diff, k = 3), 
                      data = cd4_npres_rsa, 
                      method = "REML")

pred_cd4_npres_rsa <- plot_pred.cd4(m.cd4_npres_rsa, cd4_npres_rsa, 2) 
pred_cd4_npres_rsa

## RSA presenting ##

cd4_pres_rsa <- cd4_pres %>% 
  filter(cohort == "RSA")

m.cd4_pres_rsa <- gam(cd4_trans ~ s(time_diff, k = 3), 
                     data = cd4_pres_rsa, 
                     method = "REML")

pred_cd4_pres_rsa <- plot_pred.cd4(m.cd4_pres_rsa, cd4_pres_rsa, 2) 
pred_cd4_pres_rsa

## combine them ##

pred_cd4 <- ggarrange(pred_cd4_npres_rsa, 
                      pred_cd4_npres_ch, 
                      pred_cd4_pres_rsa,
                      pred_cd4_pres_ch, 
                      ncol = 2) %>% 
  annotate_figure(bottom = "Days after ART start",
                  left = expression(sqrt(CD4)))