# Now i import each dataset
# Set the path to the folder containing the datasets
folder_path <- here("data_raw")

# Get a list of all files in the folder with the .dta extension
file_list <- list.files(folder_path, pattern = ".dta", full.names = TRUE)

# Import each dataset and assign the name based on the file name
for (file in file_list) {
  dataset_name <- tools::file_path_sans_ext(basename(file))
  assign(dataset_name, haven::read_dta(file))
}
library(foreign)
test <- read.dta("data_raw/shcs_509_tbcases.dta")
