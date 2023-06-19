#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven
)

source("preprocessing/prep_ch.R")

#### Data ####

studypopulation_tb_ch <- read_dta("data_raw/shcs_509_hivall.dta") %>% #112
  filter(disease_tb == 1,
         date_tb <=as.Date("2022-12-31"),
         date_tb >=as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16, #two patients excluded
         is.na(exitdate) | exitdate >= art_start_date) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born) #110 in studypopulation (flowchart)

#presenting with tb
ch_tb <- studypopulation_tb_ch %>% 
  filter(date_tb >= art_start_date-365 & date_tb <= art_start_date+60) #presenting with tb

write_dta(ch_tb, "data_clean/presentingTB_ch.dta") #55
write_dta(studypopulation_tb_ch, "data_clean/tb20102023_ch.dta") #110

##adding additional columns
process_data("data_clean/presentingTB_ch.dta", "presentingTB_ch")
process_data("data_clean/tb20102023_ch.dta", "tb20102023_ch")

#notpresenting with tb
all <- readRDS("data_clean/df_ch.rds")
tb_yes <- readRDS("data_clean/presentingTB_ch.rds")
ch_notb <- anti_join(all, tb_yes, by = "id")

write_rds(ch_notb, "data_clean/notpresentingTB_ch.rds")

#remove unnecassery files
file.remove("data_clean/presentingTB_ch.dta")
file.remove("data_clean/tb20102023_ch.dta")
