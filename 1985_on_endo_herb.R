#9/25/25
#making yearly endo_georefs so that qgis batch processes annual data more quickly
endo_herb_georef_allyears <- read.csv("endo_herb_georef.csv")

#################################### attempt 2 ########################################
endo_herb_georef_allyears$year <- as.numeric(endo_herb_georef_allyears$year)

cutoff_year <- 1985

filtered_endo_herb <- subset(endo_herb_georef_allyears, year >= cutoff_year)

unique_years <- unique(filtered_endo_herb$year)

for (year_val in unique_years) {
  yearly_data <- subset(filtered_endo_herb, year == year_val)
  year <- paste0(year_val, ".csv")
  write.csv(yearly_data, year, row.names = FALSE, file = )
}

###################### combining the 1985-2024 files for qgis processing ##################
write.csv(filtered_endo_herb, file= "C:/Users/malpa/OneDrive/Documents/EndoHerbQGIS/endo_herb_georef_30yr.csv")
