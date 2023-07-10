##### Libraries ----

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  dplyr, # for data wrangling
  ggplot2, # for plotting
  haven,
  lme4,
  lmerTest
)


#### data import ----

rec_ch <- readRDS("data_clean/art_lab_ch.rds")
rec_long_ch <- readRDS("data_clean/art_lab_ch_long.rds")

rec_sa <- ##
  
 # Create a dataframe with unique ids and a cohort for each
id_cohort_df <- data.frame(id = unique(rec_ch$id)) %>%
  mutate(cohort = ifelse(row_number() <= 2000, "CH", "SA"))

# Join this dataframe with the original dataframe
rec_test <- rec_long_ch %>%
  left_join(id_cohort_df, by = "id")
  
rec <- rec_test 

#### assumptions ----
ggplot(rec, aes(x= sqrt(cd4))) +
  geom_histogram() +
  theme_bw()

ggplot(rec, aes(x=log10(rna))) +
  geom_histogram() +
  theme_bw()

# Q-Q plot for cd4
ggplot(rec, aes(sample = sqrt(cd4))) +
  stat_qq() +
  stat_qq_line() +
  ggtitle("Q-Q plot for cd4")

# dependend variable doesn't have to be normally distributed, just the residuals. 
# Could still be a good idea to transform so the data isn't so skewed. 
# Have to check residuals if transformation was good


#### model ----

# For CD4

model_cd4 <- lmer(cd4 ~ time_diff + cohort + sex + (1|id), data = rec)
summary(model_cd4)

#TODO residuals plotten und auf normalität checken
#TODO spline model fitten mit mgcv package (gamm funktion benutzen)

# For RNA
model_rna <- lmer(value ~ days_since_art_start + cohort + (1|id), data = df_rna)
summary(model_rna)
