# Assessing relationship between climatic and anthropogenic global change variables
# Date: Nov 24, 2025
# Author: Joshua Fowler

library(tidyverse)
library(readxl)


library(raster)
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
##### Getting SPEI data from PRISM #################
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

ppt_monthly_10km <- exactextractr::exact_extract(ppt_stack, buffers_10, fun = "mean", append_cols = TRUE)
tmean_monthly_10km <- exactextractr::exact_extract(tmean_stack, buffers_10, fun = "mean", append_cols = TRUE)


ppt_monthly_30km <- exactextractr::exact_extract(ppt_stack, buffers_30, fun = "mean", append_cols = TRUE)
tmean_monthly_30km <- exactextractr::exact_extract(tmean_stack, buffers_30, fun = "mean", append_cols = TRUE)



colnames_key <- str_sub(colnames(ppt_monthly_10km), -4, end = nchar(colnames(ppt_monthly_10km)))
colnames(ppt_monthly_10km) <- colnames(ppt_monthly_30km) <- colnames(tmean_monthly_10km) <- colnames(tmean_monthly_30km) <- colnames_key

ppt_monthly_10km$prism_variable <- ppt_monthly_30km$prism_variable <- "ppt"; 

tmean_monthly_10km$prism_variable <- tmean_monthly_30km$prism_variable <- "tmean"; 


ppt_monthly_10km$buffer <- tmean_monthly_10km$buffer <-  "10km"
ppt_monthly_30km$buffer <- tmean_monthly_30km$buffer <-  "30km"

ppt_monthly_10km$lon <- ppt_monthly_30km$lon <- tmean_monthly_10km$lon <- tmean_monthly_30km$lon <-  coords_df$lon
ppt_monthly_10km$lat <- ppt_monthly_30km$lat <- tmean_monthly_10km$lat <- tmean_monthly_30km$lat <-  coords_df$lat




PRISM_yearly_df <- as_tibble(rbind(ppt_monthly_10km, ppt_monthly_30km, tmean_monthly_10km, tmean_monthly_30km)) %>%   
  na.omit() %>%
  pivot_longer(cols = c(-lon, -lat, -prism_variable, -buffer), names_to = "year") %>% 
  pivot_wider(names_from = prism_variable, values_from = value) 

write_csv(PRISM_yearly_df, "PRISM_yearly_df.csv")
PRISM_yearly_df <- read_csv("PRISM_yearly_df.csv")
