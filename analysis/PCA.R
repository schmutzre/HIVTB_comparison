#### libraries ####
library(tidyverse)
library(PCAmixdata)
library(gridExtra)



#### data ####

#ch#
hiv.ch <- readRDS("data_clean/art_ch.rds")
tb.ch <- readRDS("data_clean/tb_ch.rds")
incident <- hiv.ch %>% 
  filter(case_incident_2m == 1)

incidenttb <- tb.ch %>% 
  filter(case_incident_2m == 1)

#### incidence vs non-incidence ####

hiv.ch <- hiv.ch %>% 
  select(sex, age_at_ART_start, cd4_baseline, rna_baseline, who_stage, regimen, case_incident_2m)

hiv.ch$case_incident_2m <- as.factor(hiv.ch$case_incident_2m)

split <- splitmix(hiv.ch)

X1 <- scale(split$X.quanti) 
X2 <- split$X.quali 

incident_TB <- X2$case_incident_2m
X2$case_incident_2m <- NULL

# Combine X1, X2, and incident_TB
combined_df <- data.frame(X1, X2, incident_TB)

# Reorder combined_df by incident_TB
combined_df <- combined_df[order(combined_df$incident_TB), ]

# Split them again
X1_reordered <- as.matrix(combined_df[, 1:ncol(X1)])
X2_reordered <- combined_df[, (ncol(X1) + 1):(ncol(X1) + ncol(X2))]
incident_TB <- combined_df$incident_TB

# Now, perform PCAmix
res.pcamix <- PCAmix(X.quanti=X1_reordered, X.quali=X2_reordered, rename.level=TRUE, graph=FALSE)

# And then plot
plot(res.pcamix, choice="ind", 
     axes=c(1, 2),
     coloring.ind=incident_TB, 
     label=FALSE,
     posleg="topleft", 
     main="Observations")

#### suppressed at tb date vs non-suppressed at tb date (for incident cases) ####

tb.ch <- tb.ch %>% 
  filter(case_incident_2m == 1) %>% 
  select(sex, disease_tbc, risk, age_at_ART_start, region_born, who_stage, cd4_baseline, rna_baseline, regimen, tb_diag_cd4, tb_diag_rna) %>% 
  mutate(sup.TB = as.factor(case_when(tb_diag_rna < 400 & tb_diag_cd4 > 350 ~ 'sup', 
                            tb_diag_rna > 400 & tb_diag_cd4 < 350 ~ 'non.sup',
                            tb_diag_rna < 400 ~ 'rna.sup',
                            tb_diag_cd4 > 350 ~ 'cd4.sup')),
         sup.TB.cd4 = as.factor(case_when(tb_diag_cd4 > 350 ~ 'sup',
                                          TRUE ~ 'non.sup')),
         sup.TB.rna = as.factor(case_when(tb_diag_rna < 400 ~ 'sup',
                                          TRUE ~ 'non.sup')),
         cd4_baseline = sqrt(cd4_baseline),
         rna_baseline = log10(rna_baseline+1)) %>% 
  select(-c(tb_diag_cd4, tb_diag_rna))

split.tb <- splitmix(tb.ch)

X1.tb <- scale(split.tb$X.quanti) 
X2.tb <- split.tb$X.quali 

sup.TB <- X2.tb$sup.TB
X2.tb$sup.TB <- NULL

res.pcamix.tb <- PCAmix(X.quanti=X1.tb, X.quali=X2.tb, rename.level=TRUE, graph=FALSE)

# And then plot
plot(res.pcamix.tb, choice="ind", 
     axes=c(1, 2),
     coloring.ind=sup.TB, 
     label=FALSE,
     posleg="topleft", 
     main="Observations")

#### plotting different comparisons ####
#cd4

p1 <- tb.ch %>%
  ggplot(aes(x = sex, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw()

p2 <- tb.ch %>%
  ggplot(aes(x = disease_tbc, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw() +
  theme(legend.position = "none")

p3 <- tb.ch %>%
  ggplot(aes(x = risk, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw()+
  theme(legend.position = "none")

p4 <- tb.ch %>%
  ggplot(aes(x = age_at_ART_start, fill = sup.TB.cd4)) +
  geom_histogram() +
  theme_bw()+
  theme(legend.position = "none")

p5 <- tb.ch %>%
  ggplot(aes(x = region_born, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw() +
  theme(legend.position = "none")

p6 <- tb.ch %>%
  ggplot(aes(x = who_stage, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw()+
  theme(legend.position = "none")

p7 <- tb.ch %>%
  ggplot(aes(x = regimen, fill = sup.TB.cd4)) +
  geom_bar() +
  theme_bw()+
  theme(legend.position = "none")

p8 <- tb.ch %>%
  ggplot(aes(x = cd4_baseline, fill = sup.TB.cd4)) +
  geom_histogram(binwidth = 1) +
  theme_bw()+
  theme(legend.position = "none")

p9 <- tb.ch %>%
  ggplot(aes(x = rna_baseline, fill = sup.TB.cd4)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw()+
  theme(legend.position = "none")

p10 <- grid.arrange(p1,p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)

ggsave(plot = p10, filename = "results/descriptive/cd4.sup.png", width = 20, height = 8)

#rna
p1.1 <- tb.ch %>%
  group_by(sex) %>% 
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = sex, y = prop, fill = sup.TB.rna)) +
  geom_bar(stat = "identity") +
  theme_bw()

p2.1 <- tb.ch %>%
  group_by(disease_tbc) %>% 
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = disease_tbc, fill = sup.TB.rna)) +
  geom_bar() +
  theme_bw() +
  theme(legend.position = "none")

p3.1 <- tb.ch %>%
  ggplot(aes(x = risk, fill = sup.TB.rna)) +
  geom_bar() +
  theme_bw()+
  theme(legend.position = "none")


p4.1 <- tb.ch %>%
  ggplot(aes(x = age_at_ART_start, fill = sup.TB.rna)) +
  geom_histogram() +
  theme_bw()+
  theme(legend.position = "none")


p5.1 <- tb.ch %>%
  ggplot(aes(x = region_born, fill = sup.TB.rna)) +
  geom_bar() +
  theme_bw() +
  theme(legend.position = "none")


p6.1 <- tb.ch %>%
  ggplot(aes(x = who_stage, fill = sup.TB.rna)) +
  geom_bar() +
  theme_bw()+
  theme(legend.position = "none")


p7.1 <- tb.ch %>%
  ggplot(aes(x = regimen, y = prop, fill = sup.TB.rna)) +
  geom_bar(stat = "identity") +
  theme_bw()+
  theme(legend.position = "none")

p7.1
p8.1 <- tb.ch %>%
  ggplot(aes(x = cd4_baseline, fill = sup.TB.rna)) +
  geom_histogram(binwidth = 1) +
  theme_bw()+
  theme(legend.position = "none")


p9.1 <- tb.ch %>%
  ggplot(aes(x = rna_baseline, fill = sup.TB.rna)) +
  geom_histogram(binwidth = 0.2) +
  theme_bw()+
  theme(legend.position = "none")

p10.1 <- grid.arrange(p1.1, p2.1, p3.1, p4.1, p5.1, p6.1, p7.1, p8.1, p9.1, ncol = 3)

ggsave(plot = p10.1, filename = "results/descriptive/rna.sup.png", width = 20, height = 8)
