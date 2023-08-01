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

labels <- as_factor(studypopulation_tb_ch, levels="labels")

labels(studypopulation_tb_ch$tbd_pat_birth)

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

####### adding lab data to tb20102022 ####

tb20102022_ch <- readRDS("data_clean/tb20102022_ch.rds")
lab3 <- read_dta("data_raw/lab.dta")

lab3 <- lab3 %>% 
  select(id:cd4date, cd4, rna) %>% 
  filter(id %in% tb20102022_ch$id)

# First, sort the lab dataframe by id and labdate
lab3 <- lab3 %>% arrange(id, labdate)

# Add a column to count the occurrence of each id
lab3 <- lab3 %>% 
  group_by(id) %>% 
  mutate(row = row_number())

# Now reshape lab dataframe using pivot_wider
lab_wide <- lab3 %>%
  pivot_wider(names_from = row, 
              values_from = c(labdate, cd4, rna), 
              names_sep = "_") 

reshape_lab_data_ordered <- function(lab_data){
  lab_data <- lab_data %>% 
    arrange(id, labdate) %>%
    group_by(id) %>%
    mutate(row = row_number())
  
  # Split the data into separate dataframes for each measurement
  labdate_wide <- lab_data %>%
    select(id, row, labdate) %>%
    pivot_wider(names_from = row, 
                values_from = labdate, 
                names_prefix = "labdate_")
  
  cd4_wide <- lab_data %>%
    select(id, row, cd4) %>%
    pivot_wider(names_from = row, 
                values_from = cd4, 
                names_prefix = "cd4_")
  
  rna_wide <- lab_data %>%
    select(id, row, rna) %>%
    pivot_wider(names_from = row, 
                values_from = rna, 
                names_prefix = "rna_")
  
  # Combine the separate measurements back into one dataframe
  combined_wide <- reduce(list(labdate_wide, cd4_wide, rna_wide), left_join, by = "id")
  
  return(combined_wide)
}

# Use the function to reshape lab data
lab_wide_ordered <- reshape_lab_data_ordered(lab3)

# Now join with art_ch
art_ch_lab_ordered <- left_join(tb20102022_ch, lab_wide_ordered, by = "id") 

# Separate labdate, cd4, and rna columns
labdate_df <- lab_wide_ordered %>% 
  select(id, starts_with("labdate")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "labdate")

cd4_df <- lab_wide_ordered %>% 
  select(id, starts_with("cd4")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "cd4")

rna_df <- lab_wide_ordered %>% 
  select(id, starts_with("rna")) %>% 
  pivot_longer(-id, names_to = "timepoint", values_to = "rna")

cd4_df$timepoint <- gsub("cd4_", "", cd4_df$timepoint)
labdate_df$timepoint <- gsub("labdate_", "", labdate_df$timepoint)
rna_df$timepoint <- gsub("rna_", "", rna_df$timepoint)

# Join the labdate, cd4, and rna columns into one dataframe
lab_long <- labdate_df %>%
  full_join(cd4_df, by = c("id", "timepoint")) %>%
  full_join(rna_df, by = c("id", "timepoint"))

# Add a time_diff column

lab_long <- lab_long %>% 
  left_join(art_ch_lab_ordered, by = "id") %>% 
  select(id:sex, date_tb, art_start_date, cd4_group, rna_group, who_stage) %>% 
  mutate(time_diff = as.numeric(labdate - art_start_date, units = "days"))

saveRDS(lab_long, "data_clean/tb20102022long.rds")

