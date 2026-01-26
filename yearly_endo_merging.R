#11/10/25
#load in yearly nlcd data

library(dplyr)
library(readr)
library(stringr)
library(purrr)

#set path for the folder
path <- "C:/Users/malpa/OneDrive/Documents/EndoHerbQGIS/nlcd_zonal_hist_all_years"

#list out files within this folder
files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)

  
# Read and process each file
all_data <- map_df(files, function(file) {
  
  # Extract the year from the file name 
  file_year <- str_extract(basename(file), "\\d{4}")
  
  # Read in the csv
  df <- read_csv(file)
  
  # Filter for rows matching that year 
  if ("year" %in% names(df)) {
    df <- df %>% filter(year == as.numeric(file_year))
  }
  
  # add a column for the file’s year
  df <- df %>% mutate(file_year = as.numeric(file_year))
  
  return(df)
})

View(all_data)



write.csv(all_data, file = "C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_yearly_nlcd.csv")
          
View(all_data %>% filter(Sample_id == "AM_AGHY_3") %>% select(Sample_id, year, HISTO1))
          