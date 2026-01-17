# Assessing relationship between climatic and anthropogenic global change variables
# Date: Nov 24, 2025
# Author: Joshua Fowler

library(tidyverse)
library(readxl)
library(lubridate)

library(raster)
library(terra)
library(exactextractr)
library(sf)



quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

######################################################
##### Read in the data ###########
######################################################
endo_herb_data <- read_csv(file = "/Users/joshuacfowler/Documents/R_projects/Endophyte_herbarium_urbanization/Analyses/endo_herb_nit.csv")



# converting the lat long to same crs as the rasters are stored in
# define a crs
crs <- paste(
  "GEOGCRS[\"unknown\",\n    DATUM[\"North American Datum 1983\",\n        ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n            LENGTHUNIT[\"metre\",1]],\n        ID[\"EPSG\",6269]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433],\n        ID[\"EPSG\",8901]],\n    CS[ellipsoidal,2],\n        AXIS[\"longitude\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]],\n        AXIS[\"latitude\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]]"
)
# epsg6703km <- paste(
#   "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
#   "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
#   "+units=km +no_defs"
# )

locations<- endo_herb_data %>% 
  st_as_sf(coords = c("lon", "lat"), crs = crs, remove = FALSE)  %>% 
  st_transform(crs) %>% 
  dplyr::select(lon,lat, Sample_id, Institution_specimen_id, new_id, Spp_code, year, month, day)

####################################################
##### Getting climate data from PRISM ##############
####################################################
# the rprism packages isn't working because of an update to prism's file system. But I still want to use prism over chelsa because chelsa won't let me download anything past 2021
# This a folder where I have stored zipped prism data on an external drive. Note that this folder contains the unzipped prism data, and only the annual .tif files
prism_dir <- c("/Volumes/Expansion with MacOS - Data/Users/joshuacfowler/Documents/prism_climate_data/an/800m/")
ppt_files <- list.files(path = paste0(prism_dir,"ppt/" ), pattern = ".tif$")
tmean_files <- list.files(path = paste0(prism_dir,"tmean/" ), pattern = ".tif$")

# load the raster files as a stack 
# this takes alittle bit, especially if using the full set of years

ppt_stack <- terra::rast(stack(paste0(prism_dir,"ppt/", ppt_files))) 
tmean_stack <- terra::rast(stack(paste0(prism_dir,"tmean/", tmean_files)))



# extracting the monthly values at our coordinates
coords_df <- locations %>% 
  dplyr::select(lon,lat) %>% 
  distinct() 

coords <- SpatialPoints(cbind(coords_df$lon, coords_df$lat), proj4string = CRS(crs))


buffers_10 <- st_buffer(st_as_sf(coords), dist = 10000) # 10 km buffer
buffers_30 <- st_buffer(st_as_sf(coords), dist = 30000) # 30 km buffer

# extract each monthly measurement as the mean of values within each buffer
terra::gdalCache(3276) # need to expand the "gdal cache" when using the full raster stack.

ppt_annual_10km <- exactextractr::exact_extract(ppt_stack, buffers_10, fun = "mean")
tmean_annual_10km <- exactextractr::exact_extract(tmean_stack, buffers_10, fun = "mean")


ppt_annual_30km <- exactextractr::exact_extract(ppt_stack, buffers_30, fun = "mean")
tmean_annual_30km <- exactextractr::exact_extract(tmean_stack, buffers_30, fun = "mean")



colnames_key <- str_sub(colnames(ppt_annual_10km), -4, end = nchar(colnames(ppt_annual_10km)))
colnames(ppt_annual_10km) <- colnames(ppt_annual_30km) <- colnames(tmean_annual_10km) <- colnames(tmean_annual_30km) <- colnames_key

ppt_annual_10km$prism_variable <- ppt_annual_30km$prism_variable <- "ppt"; 

tmean_annual_10km$prism_variable <- tmean_annual_30km$prism_variable <- "tmean"; 


ppt_annual_10km$buffer <- tmean_annual_10km$buffer <-  "10km"
ppt_annual_30km$buffer <- tmean_annual_30km$buffer <-  "30km"

ppt_annual_10km$lon <- ppt_annual_30km$lon <- tmean_annual_10km$lon <- tmean_annual_30km$lon <-  coords_df$lon
ppt_annual_10km$lat <- ppt_annual_30km$lat <- tmean_annual_10km$lat <- tmean_annual_30km$lat <-  coords_df$lat




PRISM_yearly_df <- as_tibble(rbind(ppt_annual_10km, ppt_annual_30km, tmean_annual_10km, tmean_annual_30km)) %>%   
  na.omit() %>%
  pivot_longer(cols = c(-lon, -lat, -prism_variable, -buffer), names_to = "year") %>% 
  pivot_wider(names_from = prism_variable, values_from = value) 

write_csv(PRISM_yearly_df, "PRISM_yearly_df.csv")
PRISM_yearly_df <- read_csv("PRISM_yearly_df.csv")


##########################################################################
##### Getting nitrogen deposition data downloaded from USGS ##############
##########################################################################
# This a folder where I have stored USGS nitrogen deposition data on an external drive. Note that this folder contains the unzipped tif files with monthly values
TIN_dir <- c("/Volumes/Expansion with MacOS - Data/Users/joshuacfowler/Documents/endo_herbarium_urbanization/TIN/")
NO3_dir <- c("/Volumes/Expansion with MacOS - Data/Users/joshuacfowler/Documents/endo_herbarium_urbanization/NO3/")
NH4_dir <- c("/Volumes/Expansion with MacOS - Data/Users/joshuacfowler/Documents/endo_herbarium_urbanization/NH4/")

TIN_files <- list.files(path = TIN_dir, pattern = ".tif$")
NO3_files <- list.files(path = NO3_dir, pattern = ".tif$")
NH4_files <- list.files(path = NH4_dir, pattern = ".tif$")

# load the raster files as a stack 
# this takes alittle bit, especially if using the full set of years

TIN_stack <- terra::rast(stack(paste0(TIN_dir, TIN_files))) 
NO3_stack <- terra::rast(stack(paste0(NO3_dir, NO3_files))) 
NH4_stack <- terra::rast(stack(paste0(NH4_dir, NH4_files))) 

# calculating annual sums within each raster stack

years <- str_split_i(names(TIN_stack), "_",2)
months <- str_split_i(names(TIN_stack), "_",3)

time(TIN_stack) <- as.Date(my(paste0(months, "/", years)))
time(NO3_stack) <- as.Date(my(paste0(months, "/", years)))
time(NH4_stack) <- as.Date(my(paste0(months, "/", years)))

TIN_stack_year <- tapp(TIN_stack, "years", sum)

NO3_stack_year <- tapp(NO3_stack, "years", sum)
NH4_stack_year <- tapp(NH4_stack, "years", sum)


# extracting the yearly values at our coordinates
coords_df <- locations %>% 
  dplyr::select(lon,lat) %>% 
  distinct() 

coords <- SpatialPoints(cbind(coords_df$lon, coords_df$lat), proj4string = CRS(crs))


buffers_10 <- st_buffer(st_as_sf(coords), dist = 10000) # 10 km buffer
buffers_30 <- st_buffer(st_as_sf(coords), dist = 30000) # 30 km buffer

# extract each location's measurement as the mean of values within each buffer
terra::gdalCache(3276) # need to expand the "gdal cache" when using the full raster stack.

TIN_10km <- exactextractr::exact_extract(TIN_stack_year, buffers_10, fun = "mean")
NO3_10km <- exactextractr::exact_extract(NO3_stack_year, buffers_10, fun = "mean")
NH4_10km <- exactextractr::exact_extract(NH4_stack_year, buffers_10, fun = "mean")

TIN_30km <- exactextractr::exact_extract(TIN_stack_year, buffers_30, fun = "mean")
NO3_30km <- exactextractr::exact_extract(NO3_stack_year, buffers_30, fun = "mean")
NH4_30km <- exactextractr::exact_extract(NH4_stack_year, buffers_30, fun = "mean")




colnames_key <- str_sub(colnames(TIN_10km), -4, end = nchar(colnames(TIN_10km)))
colnames(TIN_10km) <- colnames(NO3_10km) <- colnames(NH4_10km) <- colnames(TIN_30km) <- colnames(NO3_30km) <- colnames(NH4_30km) <- colnames_key

TIN_10km$variable <- TIN_30km$variable <- "TIN"; 
NO3_10km$variable <- NO3_30km$variable <- "NO3"; 
NH4_10km$variable <- NH4_30km$variable <- "NH4"; 



TIN_10km$buffer <- NO3_10km$buffer <- NH4_10km$buffer <-  "10km"
TIN_30km$buffer <- NO3_30km$buffer <- NH4_30km$buffer <-  "30km"

TIN_10km$lon <- NO3_10km$lon <- NH4_10km$lon <- TIN_30km$lon <- NO3_30km$lon <- NH4_30km$lon <- coords_df$lon
TIN_10km$lat <- NO3_10km$lat <- NH4_10km$lat <- TIN_30km$lat <- NO3_30km$lat <- NH4_30km$lat <- coords_df$lat





nitrogen_yearly_df <- as_tibble(rbind(TIN_10km, NO3_10km, NH4_10km, TIN_30km, NO3_30km, NH4_30km)) %>%   
  na.omit() %>%
  pivot_longer(cols = c(-lon, -lat, -variable, -buffer), names_to = "year") %>% 
  pivot_wider(names_from = variable, values_from = value) 

write_csv(nitrogen_yearly_df, "nitrogen_yearly_df.csv")
# nitrogen_yearly_df.csv <- read_csv("nitrogen_yearly_df.csv")

nitrogen_mean_df <- nitrogen_yearly_df %>% 
  ungroup() %>% 
  group_by(buffer,lon,lat) %>% 
  dplyr::summarise(mean_TIN = mean(TIN),
                   mean_NO3 = mean(NO3),
                   mean_NH4 = mean(NH4))

write_csv(nitrogen_mean_df, "nitrogen_mean_df.csv")


