##### libraries ----------------------------------------------------------------

library(wesanderson)
library(gamm4)
library(jomo)
library(mitools)
library(purrr)
library(mice)
library(tidyr)
library(tidyverse)

set.seed(123)

source("utils/plot.R")

#### data preparation ----------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/cd4_ch.rds")%>% 
  mutate(id = as.factor(id),
         cohort = as.factor("CH"),
         cd4_trans = sqrt(cd4),
         time_trans = time_diff + 60,
         regimen = as.factor(regimen)) %>% 
  filter(time_diff >= -60 & time_diff < 420) %>% 
  dplyr::select(-timepoint, -cd4, -art_start_date, -date_tb, -timepoint,
                                          -pre_2016)

cd4_rsa <- readRDS("data_clean/rsa/cd4_rsa.rds") %>% 
  mutate(id = as.factor(id),
         cohort = as.factor("RSA"),
         cd4_trans = sqrt(cd4),
         time_trans = time_diff + 60,
         region = NA, 
         regimen = as.factor(regimen)) %>% 
  filter(time_diff >= -60 & time_diff < 420) %>% 
  dplyr::select(-timepoint, -cd4, -art_start_date, -date_tb, -timepoint,
         -pre_2016)

cd4 <- rbind(cd4_ch, cd4_rsa) %>% 
  mutate(time_trans2 = time_trans^2,
         time_trans3 = time_trans^3)

cd4_npres <- cd4 %>% 
  filter(presenting_tb == 0)

cd4_pres <- cd4 %>% 
  filter(presenting_tb == 1)

#### complete record analysis [interaction] ---------------------------------------------------

cd4$interaction_term <- interaction(cd4$cohort, cd4$presenting_tb)

model_cd4_count <- gamm4(cd4_trans ~ cohort + presenting_tb  +
                           s(time_trans, by = interaction_term, bs = "cs", k = 3),
                                            random = ~ (1 | id),
                                            data = cd4)
summary(model_cd4_count$gam)

pred <- pred_trend(model_cd4_count, data = cd4)

trend_cd4 <- ggplot(pred, aes(x = time_org)) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 420), breaks = seq(-60, 420, 60)) +
  scale_y_continuous(limits = c(0,650), expand = c(0,0))+
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA)) 
  
trend_cd4

ggsave(trend_cd4, filename = "results/cd4_predEff.png",
       width = 20, height = 12.7, unit = "cm")

#### complete record analysis [t^2, t^3] ---------------------------------------

model_simple_ch <- lme4::lmer(cd4_trans ~ poly(time_trans, 3) * presenting_tb + (1 |id), 
                           data = cd4_ch)

model_simple_rsa <- lme4::lmer(cd4_trans ~ poly(time_trans, 3) * presenting_tb +(1|id), 
                      data = cd4_rsa)

dtaaa_ch <- expand.grid(
  time_trans = seq(min(cd4$time_trans), max(cd4$time_trans), length.out = 20),
  presenting_tb = levels(cd4$presenting_tb), 
  id = sample(cd4_ch$id,1))

# Get predictions and intervals

preds_ch <- predictInterval(model_simple_ch, dtaaa_ch)


dtaaa_ch <- cbind(dtaaa_ch,preds_ch) %>% 
  mutate(fit_org = fit^2,
         lower_org = lwr^2,
         upper_org = upr^2,
         time = time_trans-60,
         cohort = as.factor("CH"))

dtaaa_rsa <- expand.grid(
  time_trans = seq(min(cd4$time_trans), max(cd4$time_trans), length.out = 20),
  presenting_tb = levels(cd4$presenting_tb), 
  id = sample(cd4_rsa$id,1))

preds_rsa <- predictInterval(model_simple_rsa, dtaaa_rsa)

dtaaa_rsa <- cbind(dtaaa_rsa,preds_rsa) %>% 
  mutate(fit_org = fit^2,
         lower_org = lwr^2,
         upper_org = upr^2,
         time = time_trans-60,
         cohort = as.factor("RSA"))

dtaaa_predicted <- rbind(dtaaa_ch, dtaaa_rsa) %>% 
  mutate(cohort = relevel(cohort, ref = "RSA"),
         presenting_tb = case_when(presenting_tb == 0 ~ "Without prevalent TB",
                                          TRUE ~ "With prevalent TB"))

dtaaa_predicted %>% 
  ggplot(aes(x = time, color = cohort)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lower_org, ymax = upper_org, color = cohort, fill = cohort), alpha = 0.2) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 900)) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 
  

preddd <- predict(model_simple, se.fit = TRUE, newdata = dtaaa)
predic
dtaaa <- dtaaa %>% 
  mutate(prediction = preddd)

dtaaa %>% 




summ(model_simple)

#### complete record analysis [separate models] --------------------------------

### non-presenting ###

model_cd4_count_npres <- gamm4(cd4_trans ~ cohort +
                                s(time_diff, by = cohort, bs = "ps", k = 4),
                              random = ~ (1 | id),
                              data = cd4_npres)

pred_npres <- pred_trend_uni(model_cd4_count_npres, 
                   data = cd4_npres)

ggplot(pred_npres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 800)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 

### presenting ###

model_cd4_count_pres <- gamm4(cd4_trans ~ cohort +
                           s(time_diff, by = cohort, bs = "ps", k = 4),
                         random = ~ (1 | id),
                         data = cd4_pres)

pred_pres <- pred_trend_uni(model_cd4_count_pres, 
                         data = cd4_pres)

ggplot(pred_pres, aes(x = time_diff)) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lwrP, ymax = uprP, fill = cohort), alpha = 0.2) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)")) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(face = 2, size = 10)) 



#### imputation ----------------------------------------------------------------

### Set-up ----

intervisit_Md <- function(data) {
  
  result <- data %>%
    arrange(id, date_cd4) %>% 
    group_by(id) %>% 
    mutate(diff = date_cd4 - lag(date_cd4)) %>%
    filter(!is.na(diff)) %>% 
    ungroup() %>% 
    summarise(median_diff = median(diff, na.rm = TRUE)) %>% 
    as.numeric()
  
  return(result)
  
}

fill_values_for_id <- function(patient_data, global_grid, max_value = 480, md) {
  
  # Exclude measurement variables 'time_diff' and 'cd4_value' from unique value extraction
  measurement_vars <- c("time_trans","time_diff", "cd4_trans", "date_cd4")
  unique_val_columns <- setdiff(names(patient_data), measurement_vars)
  
  # Dynamically create a tibble of unique values for each non-measurement column
  unique_vals <- patient_data %>%
    select(all_of(unique_val_columns)) %>%
    summarise(across(everything(), unique))
  
  results <- map_df(global_grid, function(grid_day) {
    # If it's the last grid point, use max_value as the upper limit
    upper_limit <- if (grid_day == max(global_grid)) {
      max_value
    } else {
      grid_day + md
    }
    
    interval_data <- patient_data %>%
      filter(time_trans >= grid_day & time_trans < upper_limit)
    
    if (nrow(interval_data) > 0) {
      selected_measurement <- sample_n(interval_data, 1)
      tibble(time_trans = selected_measurement$time_trans, cd4_trans = selected_measurement$cd4_trans)
    } else {
      tibble(time_trans = grid_day, cd4_trans = NA_real_)
    }
  })
  
  # Join the unique values with the results
  results <- cross_join(results, unique_vals)
  
  return(results)
}

### Switzerland ----

#' in Switzerland there are usually multiple cd4 measurements during the first year

## checking assumptions for applying linear imputation ##

# normality of cd4_trans

cd4_ch %>% 
  ggplot() +
  geom_histogram(aes(x = cd4_trans, y = after_stat(density)), color = "white")+
  theme_classic() +
  scale_y_continuous(expand = c(0,0))

# linearity of time on cd4_trans

cd4_ch %>% 
  ggplot(aes(x = time_diff, y = cd4_trans)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, color = "blue") + 
  geom_smooth(method = "lm", se = FALSE, color = "red") +  
  theme_classic() +
  labs(x = "Time Difference", 
       y = "CD4 Trans", 
       title = "Relationship between Time Difference and CD4 Trans")

## applying imputation ##

md_intvstCH <- intervisit_Md(cd4_ch) # 82 days

# Create a global grid from 0 to 480 days
gridCH <- seq(from = 0, to = 480, by = md_intvstCH)

gridCH_filled <- cd4_ch %>%
  group_by(id) %>%
  group_modify(~fill_values_for_id(., gridCH, md = md_intvstCH)) %>%
  ungroup() %>% 
  mutate(presenting_tb = as.numeric(presenting_tb) - 1)

Y_CH <- data.frame(gridCH_filled$cd4_trans)
X_CH <- data.frame(gridCH_filled$presenting_tb, 
                   gridCH_filled$time_trans, 
                   gridCH_filled$sex, 
                   gridCH_filled$age_at_art_start)
clus_CH <- data.frame(gridCH_filled$id)
# run jomo
class.imputed_CH <- jomo(Y = Y_CH, 
                         X = X_CH, 
                         clus = clus_CH, 
                         nburn = 1000, 
                         nbetween = 1000,
                         nimp = 20)

# split imputations, and exclude original data (imputation "0")

imp.list_CH <- imputationList(split(class.imputed_CH, 
                                    class.imputed_CH$Imputation)[-1]) 

#diese funktion mal ausprobieren

# ch_list <- lapply(imp.list_CH$imputations, function(dataset) {
#   mutate(dataset, cohort = as.factor("CH"),
#          time_diff = gridCH_filled.time_trans - 60)
# })

ch1 <- imp.list_CH$imputations$"1" %>% mutate(clus = as.character(clus),
                                              cohort = as.factor("CH"),
                                              time = gridCH_filled.time_trans - 60) %>% 
  select(clus, time, gridCH_filled.cd4_trans, gridCH_filled.presenting_tb, cohort) %>% 
  rename(id = clus,
         cd4_trans = gridCH_filled.cd4_trans)

ch2 <- imp.list_CH$imputations$"2" %>% mutate(cohort = as.factor("CH"),
                                              time_diff = gridCH_filled.time_trans - 60) %>% 
  select(clus, gridCH_filled.time_trans, gridCH_filled.cd4_trans, gridCH_filled.presenting_tb, cohort)

### South Africa ----

md_intvstRSA <- intervisit_Md(cd4_rsa) 
gridRSA <- seq(from = 0, by = md_intvstRSA, length.out = 2)
#' For South Africa we only get two intervals, which means that jomo does not 
#' make sense - just use MICE for this

gridRSA_filled <- cd4RSA %>%
  group_by(id) %>%
  group_modify(~fill_values_for_id(., gridRSA, md = md_intvstRSA)) %>%
  mutate(row_id = ifelse(row_number() == 1, "first", "second")) %>%
  ungroup() %>% 
  pivot_wider(
  names_from = row_id,
  values_from = c(time_trans, cd4_trans),
  names_sep = "_")

gridRSA_imputation <- gridRSA_filled %>% 
  select(presenting_tb, sex, age_at_art_start, regimen, cd4_trans_first, cd4_trans_second, time_trans_first, time_trans_second)

pmat_rsa <- matrix(
  c(0,1,1,1,1,1,1,1, # presenting_tb
    1,0,1,1,1,1,1,1, # sex
    1,1,0,1,1,1,1,1, # age_at_art_start
    1,1,1,0,1,1,1,1, # regimen
    1,1,1,1,0,1,1,1, # cd4_trans_first
    1,1,1,1,1,0,1,1, # cd4_trans_second
    1,1,1,1,1,1,0,1, # time_trans_second
    1,1,1,1,1,1,1,0), # time_trans_second
  nrow = 8, 
  byrow = TRUE)

# shell imputation (for post processing) #

ini_rsa <- mice(data = gridRSA_imputation, 
                maxit = 0)
post_rsa <- ini_rsa$post
post_rsa["cd4_trans_first"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post_rsa["cd4_trans_second"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"

# imputation #

K <- 20

rsa.imp <- mice(data = gridRSA_imputation,
                m = K,
                method = c("", # presenting_tb
                           "", # sex
                           "",  # age_at_art_start
                           "polyreg", # regimen
                           "norm", # cd4_first
                           "norm", # cd4_second
                           "", # time_first
                           ""), # time_second
                predictorMatrix = pmat_rsa,
                seed = 1569,
                maxit = 20,
                post = post_rsa)

# check imputations
summary(complete(rsa.imp, 0)) # with missingness, index 0
summary(complete(rsa.imp, 2)) # imputed index 1:K

# rsa1 <- complete(rsa.imp,2) %>% 
#   mutate(cohort = as.factor("RSA"),
#          id = as.character(row_number())) %>% 
#   pivot_longer(
#     cols = c(time_trans_first, time_trans_second, cd4_trans_first, cd4_trans_second),
#     names_to = c(".value"),
#     names_pattern = "(.*)_.*"
#   ) 
#   mutate(time = time_trans - 60) %>% 
#   select(id, time, cd4_trans, presenting_tb, cohort) 
#   
#              
# rsa2 <- complete(rsa.imp,2)

#### imputed record analysis ---------------------------------------------------

#' here I have to somehowe add imputation i from CH and RSA together to analyse
#' them together

#cd4_1 <- rbind(ch1, rsa1)

process_and_merge <- function(rsa_imp, ch_imp, i) {
  
  rsa_i <- complete(rsa_imp, i)
  
  # Check if required columns are present in RSA dataset
  required_cols_rsa <- c("time_trans_first", "time_trans_second", "cd4_trans_first", "cd4_trans_second")
  if (!all(required_cols_rsa %in% names(rsa_i))) {
    stop(paste("Required columns missing in RSA dataset for iteration", i))
  }
  
  rsa_i <- rsa_i %>%
    mutate(cohort = as.factor("RSA"),
           id = as.character(row_number())) %>%
    pivot_longer(
      cols = required_cols_rsa,
      names_to = c(".value"),
      names_pattern = "(.*)_.*"
    )
  
  # Check if cd4_trans column is present after pivot_longer
  if (!"cd4_trans" %in% names(rsa_i)) {
    stop(paste("cd4_trans column not found in RSA dataset after pivot for iteration", i))
  }
  
  rsa_i <- rsa_i %>% 
    select(id, time_trans, cd4_trans, presenting_tb, cohort)
  
  ch_i <- ch_imp$imputations[[as.character(i)]]
  
  # Check if required columns are present in CH dataset
  required_cols_ch <- c("clus", "gridCH_filled.time_trans", "gridCH_filled.cd4_trans")
  if (!all(required_cols_ch %in% names(ch_i))) {
    stop(paste("Required columns missing in CH dataset for iteration", i))
  }
  
  ch_i <- ch_i %>%
    mutate(clus = as.character(clus),
           cohort = as.factor("CH")) %>%
    select(clus, gridCH_filled.time_trans, gridCH_filled.cd4_trans, gridCH_filled.presenting_tb, cohort) %>% 
    rename(id = clus, cd4_trans = gridCH_filled.cd4_trans, time_trans = gridCH_filled.time_trans,
           presenting_tb = gridCH_filled.presenting_tb)
  
  # Check if cd4_trans column is present in CH dataset
  if (!"cd4_trans" %in% names(ch_i)) {
    stop(paste("cd4_trans column not found in CH dataset for iteration", i))
  }
  
  merged_data <- rbind(ch_i, rsa_i) %>% 
    mutate(presenting_tb = as.factor(presenting_tb))
  return(merged_data)
}


# Merge all datasets

all_merged_data <- lapply(1:20, function(i) process_and_merge(rsa.imp, imp.list_CH, i))

fit_and_predict_int <- function(data) {
  
  # Create interaction term
  data$interaction_term <- interaction(data$cohort, data$presenting_tb)
  
  # Fit model
  model_cd4_count <- gamm4(cd4_trans ~ cohort + presenting_tb  +
                             s(time_trans, by = interaction_term, bs = "cs", k = 3),
                           random = ~ (1 | id),
                           data = data)
  
  # Generate predictions
  pred <- pred_trend(model_cd4_count, data = data)
  
  return(pred)
}

fit_and_predict_uni <- function(data){
   
  cd4_npres <- data %>% 
    filter(presenting_tb == 0)
  
  cd4_pres <- data %>% 
    filter(presenting_tb == 1)
  
  model_cd4_count_npres <- gamm4(cd4_trans ~ cohort +
                                   s(time_diff, by = cohort, bs = "ps", k = 4),
                                 random = ~ (1 | id),
                                 data = cd4_npres)
  
  pred_npres <- pred_trend_uni(model_cd4_count_npres, 
                               data = cd4_npres)
  
  model_cd4_count_pres <- gamm4(cd4_trans ~ cohort +
                                  s(time_diff, by = cohort, bs = "ps", k = 4),
                                random = ~ (1 | id),
                                data = cd4_pres)
  
  pred_pres <- pred_trend_uni(model_cd4_count_pres, 
                              data = cd4_pres)
  
}

# Apply the function to each dataset and store predictions
all_predictions <- lapply(all_merged_data, fit_and_predict_int)

calculate_row_means_with_original_columns <- function(predictions_list) {
  # Check if the list is empty
  if (length(predictions_list) == 0) {
    return(NULL)
  }
  
  # Ensure all data frames have the same number of rows
  n_rows <- nrow(predictions_list[[1]])
  if (!all(sapply(predictions_list, nrow) == n_rows)) {
    stop("All data frames must have the same number of rows")
  }
  
  # Take the common columns from the first dataframe
  common_columns <- predictions_list[[1]][, c("time_org", "presenting_tb", "cohort", "interaction_term")]
  
  # Initialize data frame to store means
  means_df <- data.frame(mean_fit = numeric(n_rows), 
                         mean_lwr = numeric(n_rows), 
                         mean_upr = numeric(n_rows))
  
  # Loop through each row
  for (i in 1:n_rows) {
    row_values_fit <- sapply(predictions_list, function(x) x$fit_org[i], simplify = "array")
    row_values_lwr <- sapply(predictions_list, function(x) x$lwrS[i], simplify = "array")
    row_values_upr <- sapply(predictions_list, function(x) x$uprS[i], simplify = "array")
    
    means_df$mean_fit[i] <- mean(row_values_fit, na.rm = TRUE)
    means_df$mean_lwr[i] <- mean(row_values_lwr, na.rm = TRUE)
    means_df$mean_upr[i] <- mean(row_values_upr, na.rm = TRUE)
  }
  
  # Combine the common columns with the calculated means
  final_df <- cbind(common_columns, means_df)
  
  return(final_df)
}

# model fit

final_output <- calculate_row_means_with_original_columns(all_predictions)

final_output <- final_output %>% 
  mutate(presenting_tb = ifelse(presenting_tb == "Without prevalent TB", "Not presenting TB", "Presenting TB"))

trend_cd42 <- ggplot(final_output, aes(x = time_org)) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(limits = c(0,650), expand = c(0,0))+
  geom_line(aes(y = mean_fit, color = cohort)) +
  geom_ribbon(aes(ymin = mean_lwr, ymax = mean_upr, fill = cohort), alpha = 0.2) +
  labs(y = expression("CD4 count (cells/µl)"),
       x = NULL) +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title.position = "plot",
        plot.title = element_text(size = 20, hjust = 0),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank()) +
  labs(title = expression("2c | Modeled CD4 trajectories"), 
       caption = "Mean as lines, 95%-CI as ribbons")

trend_cd42

ggsave(plot = trend_cd4, filename = "results/immun/trend_cd4.png", 
       width = 16, height = 11, units = "cm")

ggsave(plot = trend_cd42, filename = "results/immun/trend_cd4_imp.png", 
       width = 16, height = 11, units = "cm")
