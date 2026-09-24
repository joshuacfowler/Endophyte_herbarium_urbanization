# Purpose: Evaluates temporally specific covariate data, and re-fits models of endophyte prevalence 
# Authors: Joshua Fowler 
# Updated: Sep 24, 2026


library(tidyverse)

library(raster)
library(terra)
library(exactextractr)
library(sf)

#reading in the TREND-nitrogen dataset, which is county level. 
# also have to read in the county info to match up to their id scheme
atm_ox <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/Atmospheric_Oxidized.txt", sep = ",")
atm_red <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/Atmospheric_Reduced.txt", sep = ",")
nit_surplus <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/NSurplus.txt", sep = ",")

county_shape <- st_read("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/County_Boundaries_2017/Byrnesetal_TREND_CountyBoundaries.shp") %>% 
  mutate(GEOID = as.integer(GEOID))


state_fips <- read.csv("./Analyses/Supp_Nitrogen_Data_Sources/state_and_county_fips_master.csv") %>% rename(GEOID = fips)

# plot(st_geometry(county_shape))

atm_join_ox <- left_join(county_shape, atm_ox) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NOx", names_to = "year") %>% mutate(year = parse_number(year))
atm_join_red <- left_join(county_shape, atm_red) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NHx", names_to = "year") %>% mutate(year = parse_number(year))
nit_surplus_join <- left_join(county_shape, nit_surplus) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NSurplus", names_to = "year") %>% mutate(year = parse_number(year))


ggplot(atm_join_ox)+
  geom_line(aes(x = year, y = NOx, group = GEOID), alpha = .01)
ggplot(atm_join_red)+
  geom_line(aes(x = year, y = NHx, group = GEOID), alpha = .01)
ggplot(nit_surplus_join)+
  geom_line(aes(x = year, y = log(log(NSurplus + 100)), group = GEOID), alpha = .01) 



# Reading in the gTrend data for the locations in our study

######################################################
##### Read in the endophyte data ###########
######################################################
endo_herb_data <- read_csv(file = "/Users/joshuacfowler/Documents/R_projects/Endophyte_herbarium_urbanization/Analyses/endo_herb_nit.csv")



# converting the lat long to same crs as the rasters are stored in
# define a crs
crs <- paste("PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n DATUM[\"North American Datum 1983\",\n    ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n   LENGTHUNIT[\"metre\",1]],\n     ID[\"EPSG\",6269]],\n     PRIMEM[\"Greenwich\",0,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8901]]],\n     CONVERSION[\"unknown\",\n     METHOD[\"Albers Equal Area\",\n     ID[\"EPSG\",9822]],\n     PARAMETER[\"Latitude of false origin\",23,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8821]],\n     PARAMETER[\"Longitude of false origin\",-96,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8822]],\n     PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8823]],\n     PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8824]],\n     PARAMETER[\"Easting at false origin\",0,\n     LENGTHUNIT[\"metre\",1],\n     ID[\"EPSG\",8826]],\n     PARAMETER[\"Northing at false origin\",0,\n     LENGTHUNIT[\"metre\",1],\n     ID[\"EPSG\",8827]]],\n     CS[\Cartesian,2],\n     AXIS[\"(E)\",east,\n     ORDER[\1],\n     LENGTHUNIT[\"metre\",1,\n     ID[\"EPSG\",9001]]],\n     AXIS[\"(N)\",north,\n     ORDER[\2],\n     LENGTHUNIT[\"metre\",1,\n     ID[\"EPSG\",9001]]]] \n     wkt:\n     [\"unknown\",\n     BASEGEOGCRS[\"unknown\",\n     DATUM[\"North American Datum 1983\",\n     ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n     LENGTHUNIT[\"metre\",1]],\n     ID[\"EPSG\",6269]],\n     PRIMEM[\"Greenwich\",0,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8901]]],\n     CONVERSION[\"unknown\",\n     METHOD[\"Albers Equal Area\",\n     ID[\"EPSG\",9822]],\n     PARAMETER[\"Latitude of false origin\",23,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8821]],\n     PARAMETER[\"Longitude of false origin\",-96,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8822]],\n     PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8823]],\n     PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n     ANGLEUNIT[\"degree\",0.0174532925199433],\n     ID[\"EPSG\",8824]],\n     PARAMETER[\"Easting at false origin\",0,\n     LENGTHUNIT[\"metre\",1],\n     ID[\"EPSG\",8826]],\n     PARAMETER[\"Northing at false origin\",0,\n     LENGTHUNIT[\"metre\",1],\n     ID[\"EPSG\",8827]]],\n     CS[\Cartesian,2],\n     AXIS[\"(E)\",east,\n     ORDER[\1],\n     LENGTHUNIT[\"metre\",1,\n     ID[\"EPSG\",9001]]],\n     AXIS[\"(N)\",north,\n     ORDER[\2],\n     LENGTHUNIT[\"metre\",1,\n     ID[\"EPSG\",9001]]]]\n")
# epsg6703km <- paste(
#   "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
#   "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
#   "+units=km +no_defs"
# )

locations<- endo_herb_data %>% 
  st_as_sf(coords = c("lon", "lat"), crs = crs, remove = FALSE)  %>% 
  st_transform(crs) %>% 
  dplyr::select(lon,lat, Sample_id, Institution_specimen_id, new_id, Spp_code, year, month, day)

# reading in tif files that have 250m resolution rasters of nitrogen surplus components
dir <- "./Analyses/Supp_Nitrogen_Data_Sources/"
gTrends_atm_ox_files <- list.files(path = paste0(dir,"Atmospheric_Oxidized/" ), pattern = ".tif$")
gTrends_nSurplus_files <- list.files(path = paste0(dir,"Surplus/" ), pattern = ".tif$")

# load the raster files as a stack 
# this takes alittle bit, especially if using the full set of years

atm_ox_stack <- terra::rast(stack(paste0(dir,"Atmospheric_Oxidized/", gTrends_atm_ox_files))) 

nSurplus_stack <- terra::rast(stack(paste0(dir,"Surplus/", gTrends_nSurplus_files))) 






# extracting the annual values at our coordinates
coords_df <- locations %>% 
  dplyr::select(lon,lat) %>% 
  distinct() 

coords <- SpatialPoints(cbind(coords_df$lon, coords_df$lat), proj4string = CRS(crs))
buffers_10 <- st_buffer(st_as_sf(coords), dist = 10000) # 10 km buffer
# buffers_30 <- st_buffer(st_as_sf(coords), dist = 30000) # 30 km buffer

# extract each monthly measurement as the mean of values within each buffer (Takes a long)
terra::gdalCache(3276) # need to expand the "gdal cache" when using the full raster stack.

atm_ox_10km <- exactextractr::exact_extract(atm_ox_stack, buffers_10, fun = "mean")
writeRDS(atm_ox_10km, "atm_ox_10km.Rds")

nSurplus_10km <- exactextractr::exact_extract(nSurplus_stack, buffers_10, fun = "mean")

