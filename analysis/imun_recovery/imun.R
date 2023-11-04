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
  filter(time_diff > -100 & time_diff < 365) %>% 
  mutate(
    group_30days = as.integer((time_diff %/% 30) + 1),
    cd4_trans = sqrt(cd4))

rna <- rbind(rna_ch %>% select(-timepoint),
             rna_rsa %>% select(-timepoint)) %>% 
  filter(time_diff > -100 & time_diff < 365) %>% 
  mutate(group_30days = as.integer((time_diff %/% 30) + 1),
         rna_trans = log10(rna + 1))

#### models --------------------------------------------------------------------

### cd4 ###

m.cd4 <- gam(cd4_trans ~ s(time_diff, k = 4, by = cohort), 
             data = cd4, 
             method = "REML")

lm_cd4 <- lmer(sqrt(cd4) ~ time_diff + (1 | id), data = cd4)

### rna ###

m.rna <- gam(rna_trans ~ s(time_diff, k = 5, by = cohort), 
             data = rna, 
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
  
se.cd4

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

### raw ###

# cd4 #

cd4_raw <- cd4 %>%
  ggplot(aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = sqrt(cd4)), color = "grey", alpha = .2) +
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "", title = "raw data + GAM") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5))

print(cd4_raw)

# rna #

rna_raw <- rna %>%
  ggplot(aes(x = time_diff)) +
  geom_line(aes(group = factor(id), y = log10(rna+1)), color = "grey", alpha = .2) +
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() +
  labs(x = "Days since ART start", title = "") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5))+
  facet_wrap(~cohort) 

print(rna_raw)

### fit ###

# I calculate the simultaneous CI instead of point-wise https://fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/

# cd4 #

trend_cd4 <- plot_trend(m.cd4, cd4) +
  geom_hline(yintercept = sqrt(350)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(x = "Days since ART start", y = expression(sqrt(CD4))) +
  coord_cartesian(xlim = c(-60, 360), ylim = c(0, 30)) +
  facet_wrap(~cohort) +
  theme_minimal()

trend_cd4

# rna #

trend_rna <- plot_trend(m.rna, rna) +
  geom_hline(yintercept = log10(400)) +
  geom_vline(xintercept = 0, linetype = "dotted")+
  labs(x = "Days since ART start", y = "log10(viral-load)") +
  coord_cartesian(xlim = c(-60, 360), ylim = c(0, 7.5)) +
  facet_wrap(~cohort) +
  theme_minimal()

trend_rna

### raw + fit ###

trend.cd4.raw <- plot_trend.raw.cd4(m.cd4, cd4) +
  facet_wrap(~cohort) 
  
trend.cd4.raw

ggsave(plot = trend.cd4.raw, filename = "results/slopes/cd4.png", 
       width = 16, height = 11, units = "cm")

trend.rna.raw <- plot_trend.raw.rna(m.rna, rna) +
  theme_classic() +
  facet_wrap(~cohort) 

trend.rna.raw
  
ggsave(plot = trend.rna.raw, filename = "results/slopes/rna.png", 
       width = 16, height = 11, units = "cm")

