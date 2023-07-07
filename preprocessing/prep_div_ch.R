#### Library #### 

if(!require(pacman)) install.packages("pacman")

pacman:: p_load(
  tidyverse,
  haven
)

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
  filter(art_start_date <=as.Date("2022-12-31"),
         art_start_date >=as.Date("2010-01-01"),
         !is.na(sex),
         !is.na(born),
         year(art_start_date) - born >= 16,
         is.na(exitdate) | exitdate >= art_start_date) %>% 
  mutate(age_at_ART_start = year(art_start_date)-born) %>% 
  filter(date_tb >= art_start_date-365 & date_tb <= art_start_date+60) #presenting with tb

##adding additional columns
process_data(studypopulation_tb_ch, "tb20102022_ch")
process_data(ch_tb, "presentingTB_ch")

#notpresenting with tb
ch_notb <- readRDS("data_clean/art_ch.rds") %>% 
  filter(disease_tb == 0)

write_rds(ch_notb, "data_clean/notpresentingTB_ch.rds")

