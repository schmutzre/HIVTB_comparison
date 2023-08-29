#### libraries ####
library(tidyverse)
library(PCAmixdata)

#### data ####

#ch#
hiv.ch <- readRDS("data_clean/art_ch.rds")
tb.ch <- readRDS("data_clean/tb_ch.rds")
incident <- hiv.ch %>% 
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

#### suppressed at tb date vs non-suppressed at tb date ####


tb.ch <- tb.ch %>% 
  select(sex, disease_tbc, age_at_ART_start, cd4_baseline, rna_baseline, regimen, tb_diag_cd4, tb_diag_rna) %>% 
  mutate(sup.TB = as.factor(case_when(tb_diag_rna < 400 & tb_diag_cd4 > 350 ~ 'sup', 
                            tb_diag_rna > 400 & tb_diag_cd4 < 350 ~ 'non.sup',
                            tb_diag_rna < 400 ~ 'rna.sup',
                            tb_diag_cd4 > 350 ~ 'cd4.sup')),
         cd4_baseline = sqrt(cd4_baseline),
         rna_baseline = log10(rna_baseline)) %>% 
  select(-c(tb_diag_cd4, tb_diag_rna))

split.tb <- splitmix(tb.ch)

X1.tb <- scale(split.tb$X.quanti) 
X2.tb <- split.tb$X.quali 

sup.TB <- X2.tb$sup.TB
X2.tb$sup.TB <- NULL

res.pcamix.tb <- PCAmix(X.quanti=X1.tb, X.quali=X2.tb, rename.level=TRUE, graph=FALSE)

# And then plot
plot(res.pcamix.tb, choice="ind", 
     axes=c(2, 3),
     coloring.ind=sup.TB, 
     label=FALSE,
     posleg="topleft", 
     main="Observations")

