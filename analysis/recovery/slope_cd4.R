##### libraries ----------------------------------------------------------------

library(wesanderson)
library(gamm4)
library(jomo)
library(mitools)
library(purrr)
library(mice)
library(tidyr)
library(tidyverse)
library(ggeffects)
library(data.table)
library(parallel)

set.seed(123)

#source("utils/plot.R")

#### data preparation ----------------------------------------------------------

cd4_ch <- readRDS("data_clean/ch/lab_ch.rds") %>% 
  mutate(
         cohort = as.factor("CH"),
         cd4_trans = sqrt(cd4),
         rna_trans = log10(rna +1),
         time_trans = time_diff + 60,
         regimen = as.factor(regimen)) %>% 
  filter(time_diff >= -60 & time_diff < 420)  %>% 
  dplyr::select(id, cohort, labdate, presenting_tb, cd4_baseline, time_diff, cd4_trans, time_trans, age_at_art_start, sex, rna_trans)

cd4_ch <- as.data.table(cd4_ch)

cd4_rsa <- readRDS("data_clean/rsa/lab_both_rsa.rds") %>% 
  mutate(
         cohort = as.factor("RSA"),
         cd4_trans = sqrt(cd4),
         rna_trans = log10(rna+1),
         time_trans = time_diff + 60,
         region = NA, 
         regimen = as.factor(regimen)) %>% 
  filter(time_diff >= -60 & time_diff < 420) %>% 
  rename(sex = gender,
         labdate = date) %>% 
  dplyr::select(id, cohort, labdate, presenting_tb, cd4_baseline, time_diff, cd4_trans, time_trans, age_at_art_start, sex, rna_trans) 

cd4_rsa <- as.data.table(cd4_rsa)

# cd4 <- rbind(cd4_ch, cd4_rsa) %>% 
#   mutate(time_trans2 = time_trans^2,
#          time_trans3 = time_trans^3) %>% 
#   mutate(time_trans_scaled = scale(time_trans))

# cd4_npres <- cd4 %>% 
#   filter(presenting_tb == 0)
# 
# cd4_pres <- cd4 %>% 
#   filter(presenting_tb == 1)
# 
# #### complete record analysis [interaction] ---------------------------------------------------
# 
# cd4$interaction_term <- interaction(cd4$cohort, cd4$presenting_tb)
# 
# model_cd4_count <- gamm4(cd4_trans ~ cohort + presenting_tb  +
#                            s(time_trans, by = interaction_term, bs = "cs", k = 3),
#                                             random = ~ (1 | id),
#                                             data = cd4)
# 
# pred <- pred_trend(model_cd4_count, data = cd4)
# 
# trend_cd4 <- ggplot(pred, aes(x = time_org)) +
#   facet_wrap(~presenting_tb) +
#   scale_x_continuous(expand = c(0,0), limits = c(-60, 420), breaks = seq(-60, 420, 60)) +
#   scale_y_continuous(limits = c(0,650), expand = c(0,0))+
#   geom_line(aes(y = fit_org, color = cohort)) +
#   geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
#   labs(x = "Days after ART start",
#        y = expression("CD4 count (cells/µl)")) +
#   geom_hline(yintercept = 350, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) +
#   theme_classic(base_size = 10) +
#   theme(legend.position = "bottom", legend.title = element_blank(),
#         panel.spacing.x = unit(.66, "cm"),
#         plot.title.position = "plot",
#         plot.title = element_text(face = 2, size = 10),
#         panel.background = element_rect(fill='transparent'), #transparent panel bg
#         plot.background = element_rect(fill='transparent', color=NA)) 
#   
# trend_cd4
# 
# #ggsave(trend_cd4, filename = "results/cd4_inter.png",
#  #      width = 20, height = 12.7, unit = "cm")
# 
# #### complete record analysis [t^2, t^3] interaction term ----------------------
# 
# model_simple <- lme4::lmer(cd4_trans ~ poly(time_trans_scaled, 3) + cohort * presenting_tb + (1 |id), 
#                               data = cd4)
# 
# # Create parameter grid for predictions
# parameters <- expand.grid(
#   time_trans_scaled = seq(min(cd4$time_trans_scaled), max(cd4$time_trans_scaled), length.out = 20),
#   presenting_tb = levels(cd4$presenting_tb), 
#   cohort = levels(cd4$cohort),
#   id = sample(cd4$id, 1)  # Use cd4$id instead of cd4_ch$id to avoid errors
# )
# 
# # Get predictions and intervals
# predictions_simple <- merTools::predictInterval(model_simple, parameters)
# 
# # Calculate the original scale of time_trans
# parameters <- parameters %>%
#   mutate(time_trans = time_trans_scaled * attr(scale(cd4$time_trans), 'scaled:scale') + attr(scale(cd4$time_trans), 'scaled:center'))
# 
# # Combine predictions with parameters and re-scale for plotting
# parameters <- cbind(parameters, predictions_simple) %>%
#   mutate(fit_org = fit^2,
#          lower_org = lwr^2,
#          upper_org = upr^2,
#          time = time_trans - 60,
#          cohort = relevel(cohort, ref = "RSA"),
#          presenting_tb = case_when(presenting_tb == 0 ~ "Without prevalent TB",
#                                    TRUE ~ "With prevalent TB"))
# 
# # Plot the results
# parameters %>%
#   ggplot(aes(x = time, color = cohort)) +
#   geom_line(aes(y = fit_org, color = cohort)) +
#   geom_ribbon(aes(ymin = lower_org, ymax = upper_org, color = cohort, fill = cohort), alpha = 0.2) +
#   facet_wrap(~presenting_tb) +
#   scale_x_continuous(expand = c(0, 0), limits = c(-60, 420), breaks = seq(-60, 420, 60)) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 1000)) +
#   labs(x = "Days after ART start",
#        y = expression("CD4 count (cells/µl)")) +
#   geom_hline(yintercept = 350, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) +
#   theme_classic(base_size = 10) +
#   theme(legend.position = "top", legend.title = element_blank(),
#         panel.spacing.x = unit(.66, "cm"),
#         plot.title.position = "plot",
#         plot.title = element_text(face = 2, size = 10))

#### complete record analysis [t^2, t^3] separate models -----------------------

model_simple_ch <- lme4::lmer(cd4_trans ~ poly(time_trans, 2) * presenting_tb + (1 |id), 
                           data = cd4_ch)

model_simple_rsa <- lm(cd4_trans ~ poly(time_trans, 2) * presenting_tb, 
                      data = cd4_rsa)

pred_ch <- predict_response(model_simple_ch, terms = c("time_trans [all]", "presenting_tb")) %>% as.data.table()  %>% mutate(cohort = as.factor("CH"))
pred_rsa <- predict_response(model_simple_rsa, terms = c("time_trans [all]", "presenting_tb")) %>% as.data.table() %>%  mutate(cohort = as.factor("RSA"))

pred_both <- rbind(pred_ch, pred_rsa) %>% 
  mutate(fit_org = predicted^2,
         lower_org = conf.low^2,
         upper_org = conf.high^2,
         time = x-60,
         presenting_tb = group,
         cohort = relevel(cohort, ref = "RSA"),
         presenting_tb = case_when(presenting_tb == 0 ~ "Without prevalent TB",
                                   TRUE ~ "With prevalent TB"))

pred_comp <- pred_both %>% 
  ggplot(aes(x = time, color = cohort)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lower_org, ymax = upper_org, color = cohort, fill = cohort), alpha = 0.2) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 750)) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       title = "Original dataset") +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title = element_text(hjust = 0.5))
pred_comp
#### complete record analysis [separate models] --------------------------------

### Switzerland ###

# model_gam_ch <- gamm4(cd4_trans ~ presenting_tb +
#                               s(time_trans, by = presenting_tb, bs = "ps", k = 4),
#                               random = ~ (1 | id),
#                               data = cd4_ch)
# 
# pred_ch <- pred_trend_uni(model_gam_ch, 
#                    data = cd4_ch)
# 
# pred_ch <- pred_ch %>% mutate(cohort = as.factor("CH"))
# 
# ggplot(pred_ch, aes(x = time_org)) +
#   facet_wrap(~presenting_tb) +
#   scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 800)) +
#   geom_line(aes(y = fit_org)) +
#   geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
#   labs(x = "Days after ART start",
#        y = expression("CD4 count (cells/µl)")) +
#   geom_hline(yintercept = 350, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) +
#   theme_classic(base_size = 10) +
#   theme(legend.position = "top", legend.title = element_blank(),
#         panel.spacing.x = unit(.66, "cm"),
#         plot.title.position = "plot",
#         plot.title = element_text(face = 2, size = 10)) 
# 
# ### South Africa ###
# 
# model_gam_rsa <- gamm4(cd4_trans ~ presenting_tb +
#                         s(time_trans, by = presenting_tb, bs = "ps", k = 4),
#                       random = ~ (1 | id),
#                       data = cd4_rsa)
# 
# pred_rsa <- pred_trend_uni(model_gam_rsa, 
#                           data = cd4_rsa)
# 
# pred_rsa <- pred_rsa %>% mutate(cohort = as.factor("RSA"))
# 
# ggplot(pred_rsa, aes(x = time_org)) +
#   facet_wrap(~presenting_tb) +
#   scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 800)) +
#   geom_line(aes(y = fit_org)) +
#   geom_ribbon(aes(ymin = lwrS, ymax = uprS), alpha = 0.2) +
#   labs(x = "Days after ART start",
#        y = expression("CD4 count (cells/µl)")) +
#   geom_hline(yintercept = 350, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) +
#   theme_classic(base_size = 10) +
#   theme(legend.position = "top", legend.title = element_blank(),
#         panel.spacing.x = unit(.66, "cm"),
#         plot.title.position = "plot",
#         plot.title = element_text(face = 2, size = 10)) 
# 
# #combine
# 
# pred3 <- rbind(pred_ch, pred_rsa) %>% mutate(cohort = fct_relevel(cohort, "RSA"))
# 
# ggplot(pred3, aes(x = time_org)) +
#   facet_wrap(~presenting_tb) +
#   scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 800)) +
#   geom_line(aes(y = fit_org, color = cohort)) +
#   geom_ribbon(aes(ymin = lwrS, ymax = uprS, fill = cohort), alpha = 0.2) +
#   labs(x = "Days after ART start",
#        y = expression("CD4 count (cells/µl)")) +
#   geom_hline(yintercept = 350, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#   scale_color_manual(values = wes_palette("Moonrise2")) +
#   scale_fill_manual(values = wes_palette("Moonrise2")) +
#   theme_classic(base_size = 10) +
#   theme(legend.position = "top", legend.title = element_blank(),
#         panel.spacing.x = unit(.66, "cm"),
#         plot.title.position = "plot",
#         plot.title = element_text(face = 2, size = 10)) 

#### imputation ----------------------------------------------------------------

### Set-up ----

intervisit_Md <- function(data) {
  
  result <- data %>%
    arrange(id, labdate) %>% 
    group_by(id) %>% 
    mutate(diff = labdate - lag(labdate)) %>%
    filter(!is.na(diff)) %>% 
    ungroup() %>% 
    summarise(median_diff = median(diff, na.rm = TRUE)) %>% 
    as.numeric()
  
  return(result)
  
}

fill_values_for_id <- function(patient_data, global_grid, max_value = 480, md) {
  
  # Exclude measurement variables 'time_diff' and 'cd4_value' from unique value extraction
  measurement_vars <- c("time_trans", "time_diff", "cd4_trans", "labdate", "rna_trans")
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
      tibble(time_trans = selected_measurement$time_trans, 
             cd4_trans = selected_measurement$cd4_trans,
             rna_trans = selected_measurement$rna_trans)
    } else {
      random_time_trans <- runif(1, grid_day, upper_limit)
      tibble(time_trans = random_time_trans, cd4_trans = NA_real_)
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

md_intvstCH <- intervisit_Md(cd4_ch) # 85 days

# Create a global grid from 0 to 420 days
gridCH <- seq(from = 0, to = 420, by = md_intvstCH)

gridCH_filled <- cd4_ch %>%
  group_by(id) %>%
  group_modify(~fill_values_for_id(., gridCH, md = md_intvstCH)) %>%
  ungroup() %>% 
  mutate(presenting_tb = as.numeric(presenting_tb) - 1)

Y_CH <- data.frame(gridCH_filled[, c("cd4_trans")])
X_CH <- data.frame(gridCH_filled[, c("presenting_tb", "time_trans", "sex", "age_at_art_start", "cohort")])
clus_CH <- data.frame(gridCH_filled[, c("id")])

# run jomo
class.imputed_CH <- jomo(Y = Y_CH, 
                         X = X_CH, 
                         clus = clus_CH, 
                         nburn = 1000, 
                         nbetween = 1000,
                         nimp = 20)
######
# split imputations, and exclude original data (imputation "0")

imp.list_CH <- imputationList(split(class.imputed_CH, 
                                    class.imputed_CH$Imputation)[-1]) 

predictions.ch <- lapply(seq_along(imp.list_CH$imputations), function(i) {
  # Fit the model to the ith imputed dataset
  m <- lmer(cd4_trans ~ poly(time_trans, 2) * presenting_tb + (1 | clus), 
            data = imp.list_CH$imputations[[i]])
  
  # Generate predictions for the ith imputed dataset
  predict_response(m, terms = c("time_trans [all]", "presenting_tb"))
})

pooled_predictions.ch <- pool_predictions(predictions.ch) %>% as_tibble() %>% 
  mutate(fit_org = predicted^2,
                 lower_org = conf.low^2,
                 upper_org = conf.high^2,
                 time = x-60,
                 cohort = as.factor("CH"),
         presenting_tb = group,
         cohort = as.factor("CH"))

saveRDS(pooled_predictions.ch, "results/pooled_predictions_cd4_ch.rds")

pooled_predictions.ch %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = fit_org)) +
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

### South Africa ----

md_intvstRSA <- intervisit_Md(cd4_rsa) 
gridRSA <- seq(from = 0, by = md_intvstRSA, length.out = 2)
#' For South Africa we only get two intervals, which means that jomo does not 
#' make sense - just use MICE for this

gridRSA_filled <- cd4_rsa %>%
  group_by(id) %>%
  group_modify(~fill_values_for_id(., gridRSA, md = md_intvstRSA)) %>%
  mutate(row_id = ifelse(row_number() == 1, "first", "second")) %>%
  ungroup() %>% 
  pivot_wider(
  names_from = row_id,
  values_from = c(time_trans, cd4_trans, rna_trans),
  names_sep = "_")

gridRSA_imputation <- gridRSA_filled %>% 
  select(id, cohort, presenting_tb, sex, age_at_art_start, cd4_trans_first, cd4_trans_second, time_trans_first, time_trans_second, rna_trans_first, rna_trans_second)

# shell imputation (for post processing) #

ini.rsa <- mice(data = gridRSA_imputation, 
                maxit = 0)
pred.rsa <- ini.rsa$predictorMatrix 
pred.rsa[, c('id', 'cohort')] <- 0
post.rsa <- ini.rsa$post
post.rsa["cd4_trans_first"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post.rsa["cd4_trans_second"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 50))"
post.rsa["rna_trans_first"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 6))"
post.rsa["rna_trans_second"] <- "imp[[j]][, i] <- squeeze(imp[[j]][, i], c(0, 6))"

# imputation #

K <- 20

rsa.imp <- mice(data = gridRSA_imputation,
                m = K,
                method = c("", # id
                           "", # cohort
                           "", # presenting_tb
                           "", # sex
                           "",  # age_at_art_start
                           "norm", # cd4_first
                           "norm", # cd4_second
                           "", # time_first
                           "", # time_second
                           "norm", # rna first
                           "norm"), # rna second
                predictorMatrix = pred.rsa,
                seed = 1569,
                maxit = 20,
                post = post.rsa, printFlag = F)

# check imputations
summary(complete(rsa.imp, 0)) # with missingness, index 0
summary(complete(rsa.imp, 2)) # imputed index 1:K

#### imputed record analysis ---------------------------------------------------

predictions.rsa <- lapply(1:K, function(i) {
  
  data <- mice::complete(rsa.imp, action = i) %>% 
    pivot_longer(cols = c(cd4_trans_first, cd4_trans_second, time_trans_first, time_trans_second),
                 names_to = c(".value", "measurement"),
                 names_pattern = "(.*)_(.*)")
  
  # Fit the model to the ith imputed dataset
  m <- lm(cd4_trans ~ poly(time_trans, 2) * presenting_tb, 
            data = data)
  
  # Generate predictions for the ith imputed dataset
  predict_response(m, terms = c("time_trans [all]", "presenting_tb"))
})

pooled_predictions.rsa <- pool_predictions(predictions.rsa) %>% as_tibble() %>% 
  mutate(fit_org = predicted^2,
         lower_org = conf.low^2,
         upper_org = conf.high^2,
         time = x-60,
         cohort = as.factor("CH"),
         presenting_tb = group,
         cohort = as.factor("RSA"))

pooled_predictions.rsa %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = fit_org)) +
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

pooled_predictions.ch <- readRDS("results/pooled_predictions_cd4_ch.rds")
pooled_predictions.rsa <- readRDS("results/pooled_predictions.rsa.cd4.rds")

pooled_predictions <- rbind(pooled_predictions.ch, pooled_predictions.rsa) %>% 
  mutate(cohort = relevel(cohort, ref = "RSA"),
         presenting_tb = case_when(presenting_tb == 0 ~ "Without prevalent TB",
                                   TRUE ~ "With prevalent TB"))

pred_imp <- pooled_predictions %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = fit_org, color = cohort)) +
  geom_ribbon(aes(ymin = lower_org, ymax = upper_org, color = cohort, fill = cohort), alpha = 0.2) +
  facet_wrap(~presenting_tb) +
  scale_x_continuous(expand = c(0,0), limits = c(-60, 360), breaks = seq(-60, 360, 60)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 750)) +
  labs(x = "Days after ART start",
       y = expression("CD4 count (cells/µl)"),
       title = "Imputed dataset") +
  geom_hline(yintercept = 350, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_color_manual(values = wes_palette("Moonrise2")) +
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.spacing.x = unit(.66, "cm"),
        plot.title = element_text(hjust = 0.5))

plot <- ggpubr::ggarrange(pred_comp, pred_imp, ncol = 1,  common.legend = TRUE, legend = "bottom")
plot
ggsave(plot, filename = "results/cd4_slope.png",
       width = 16, height = 12, unit = "cm")
