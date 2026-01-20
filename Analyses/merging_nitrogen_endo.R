#Mallory Tucker and Josh Fowler
#4/29/24
#Merging average yearly deposited nitrogen from 2000-2020 to create an average amount of nitrogen deposition for each specimen

#Load in (updated) endoherb data
endo_herbnew <- read.csv(file = "C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_georef.csv")


library(tidyverse) # for data manipulation and ggplot
# library(slider) # add on to tidyverse to calculate sliding windows for climate data
library(fuzzyjoin) # add on to tidyverse to merge tables on nearest values
library(readxl)
library(lubridate)
library(ggmap)
library(prism) # to import prism raster files
library(terra)
library(INLA) # using this here to work with the shape boundary
library(raster)

#making a raster stack of NO3 data
# nitrogen_mean = average nitrogen deposited 2000-2020
setwd("C:/Users/malpa/Box/endo_urb/NO3/NO3")
n_rasters <- list.files(pattern=".tif")
raster_stack <- stack(n_rasters)
NO3 <- mean(raster_stack)
#writeRaster(NO3, filename = "mean_NO3", format = "GTiff")


#######################################Now extracting the NO3 data at our points##################################
endo_herb <- endo_herbnew %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" & !is.na(County)) 

coords_df <- endo_herb %>% 
  distinct(lat,lon) %>% 
  filter(!is.na(lat), !is.na(lon))
coords <- cbind(coords_df$lon, coords_df$lat)

NO3_mean <- extract(NO3, coords)
NO3_meanpts <- cbind(coords,NO3_mean)
colnames(NO3_meanpts) <- c("lon", "lat", "NO3_mean")

endo_herb <- merge(endo_herb, NO3_meanpts, by = c("lon","lat"))
write.csv(endo_herb, file = "C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_nit.csv")



############################################### NH4 ####################################################

#making a raster stack of NH4 data
setwd("C:/Users/malpa/Downloads/NH4")
nh4_rasters <- list.files(pattern=".tif")
raster_stacknh4 <- stack(nh4_rasters)
NH4 <- mean(raster_stacknh4)
#writeRaster(NH4, filename = "mean_NH4", format = "GTiff")


######################################Now extracting the NH4 data at our points##################################

NH4_mean <- extract(NH4, coords)
NH4_meanpts <- cbind(coords,NH4_mean)
colnames(NH4_meanpts) <- c("lon", "lat", "NH4_mean")

endo_herb <- merge(endo_herb, NH4_meanpts, by = c("lon","lat"))
write.csv(endo_herb, file = "C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_nit.csv")

######################################## Total inorganic Nitrogen (TIN) #################################################

#making a raster stack of Total N data
setwd("C:/Users/malpa/Downloads/TIN")
tin_rasters <- list.files(pattern=".tif")
raster_stack_tin <- stack(tin_rasters)
TIN <- mean(raster_stack_tin)
#writeRaster(TIN, filename = "mean_TIN", format = "GTiff")

######################################Now extracting the TIN data at our points##################################

TIN_mean <- extract(TIN, coords)
TIN_meanpts <- cbind(coords,TIN_mean)
colnames(TIN_meanpts) <- c("lon", "lat", "TIN_mean")

endo_herb <- merge(endo_herb, TIN_meanpts, by = c("lon","lat"))
write.csv(endo_herb, file = "C:/Users/malpa/OneDrive/Documents/Endophyte_herbarium_urbanization/endo_herb_nit.csv")
