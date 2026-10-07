# Purpose: Evaluates temporally specific covariate data, and re-fits models of endophyte prevalence 
# Authors: Joshua Fowler 
# Updated: Sep 24, 2026


library(devtools)
# library("devtools")
# devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE, force = TRUE)
#INLA relies on Rgraphviz (and other packages, you can use Bioconductor to help install)
library(dplyr)
library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
# devtools::install_github('timcdlucas/INLAutils')
# library(INLAutils) # supposedly has a function to plot residuals, might not need?
library(inlabru)
library(fmesher)

library(sf)
library(rmapshaper)
library(terra)
library(tidyterra)
#issue with vctrs namespace, .0.60 is loaded but need updated version

library(tidybayes) # using this for dotplots
library(GGally)
library(patchwork)
library(metR)
library(egg) # for labelling panels
library(ggmap)
library(pROC)
library(ggplot2)
library(maps)



library(raster)
library(exactextractr)


invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")

species_codes <- c("AGHY", "AGPE", "ELVI")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")





#reading in the TREND-nitrogen dataset, which is county level. 
# also have to read in the county info to match up to their id scheme

# Note that units are in kilograms of nitrogen per hectare per year, while our main analysis uses kilograms of nitrogen per km ^2 per year
# There are 100 hectare in 1km^2
atm_ox <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/Atmospheric_Oxidized.txt", sep = ",")
atm_red <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/Atmospheric_Reduced.txt", sep = ",")
nit_surplus <- read.delim("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/NSurplus.txt", sep = ",")

county_shape <- st_read("./Analyses/Supp_Nitrogen_Data_Sources/TREND-nitrogen/TREND-nitrogen/County_Boundaries_2017/Byrnesetal_TREND_CountyBoundaries.shp") %>% 
  mutate(GEOID = as.integer(GEOID))


state_fips <- read.csv("./Analyses/Supp_Nitrogen_Data_Sources/state_and_county_fips_master.csv") %>% rename(GEOID = fips)

# plot(st_geometry(county_shape))

atm_ox_join <- left_join(county_shape, atm_ox) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NOx", names_to = "year") %>% mutate(year = parse_number(year))
atm_red_join <- left_join(county_shape, atm_red) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NHx", names_to = "year") %>% mutate(year = parse_number(year))
nit_surplus_join <- left_join(county_shape, nit_surplus) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "NSurplus", names_to = "year") %>% mutate(year = parse_number(year))

NDep_join <- atm_ox_join %>% st_drop_geometry() %>% left_join(atm_red_join %>% st_drop_geometry()) %>% left_join(nit_surplus_join %>% st_drop_geometry()) %>% 
  mutate(NDep = NOx + NHx) %>% mutate(across(c(NDep, NOx, NHx, NSurplus), ~ .x*100))

ggplot(NDep_join)+
  geom_line(aes(x = year, y = NOx, group = GEOID), alpha = .01)
ggplot(NDep_join)+
  geom_line(aes(x = year, y = NHx, group = GEOID), alpha = .01)
ggplot(NDep_join)+
  geom_line(aes(x = year, y = NDep, group = GEOID), alpha = .01)
ggplot(NDep_join)+
  geom_line(aes(x = year, y = NSurplus, group = GEOID), alpha = .01) 



NDep_normal <- NDep_join %>% 
  filter(year >= 1999) %>% 
  group_by(STATEFP, COUNTYFP, COUNTYNS, AFFGEOID, GEOID, NAME, LSAD, StateFIPS, CountyFIPS, AREA..HA., name,   state) %>% 
  summarize(NOx_Normal = mean(NOx),
            NHx_Normal = mean(NHx),
            NDep_Normal = mean(NDep),
            NSurplus_Normal = mean(NSurplus)) 

NDep_change <- NDep_join %>% 
  ungroup() %>% 
  filter(year == max(year) | year == min(year)) %>% 
  pivot_wider(names_from = year, values_from = c(NOx, NHx, NDep, NSurplus)) %>% 
  mutate(NOx_change = NOx_2017 - NOx_1930,
         NHx_change = NHx_2017 - NHx_1930,
         NDep_change = NDep_2017 - NDep_1930,
         NSurplus_change = NSurplus_2017 - NSurplus_1930) %>% 
  left_join(NDep_normal)



ggplot(NDep_change)+
  geom_point(aes(x = NOx_2017, y = NOx_change))
ggplot(NDep_change)+
  geom_point(aes(x = NOx_Normal, y = NOx_change))

ggplot(NDep_change)+
  geom_point(aes(x = NHx_2017, y = NHx_change))
ggplot(NDep_change)+
  geom_point(aes(x = NHx_Normal, y = NHx_change))

ggplot(NDep_change)+
  geom_point(aes(x = NDep_2017, y = NDep_change))
ggplot(NDep_change)+
  geom_point(aes(x = NDep_Normal, y = NDep_change))


ggplot(NDep_change)+
  geom_point(aes(x = NSurplus_2017, y = NSurplus_change))
ggplot(NDep_change)+
  geom_point(aes(x = NSurplus_Normal, y = NSurplus_change))


# # Reading in the gTrend data for the locations in our study
# 
# ######################################################
# ##### Read in the endophyte data ###########
# ######################################################
# endo_herb_data <- read_csv(file = "/Users/joshuacfowler/Documents/R_projects/Endophyte_herbarium_urbanization/Analyses/endo_herb_nit.csv")
# 
# 
# 
# # converting the lat long to same crs as the rasters are stored in
# # define a crs
# # crs <- paste("PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"North American Datum 1983\",\n            ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n                LENGTHUNIT[\"metre\",1]],\n            ID[\"EPSG\",6269]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8901]]],\n    CONVERSION[\"unknown\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")
# crs <- paste("GEOGCRS[\"unknown\",\n    DATUM[\"North American Datum 1983\",\n        ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n            LENGTHUNIT[\"metre\",1]],\n        ID[\"EPSG\",6269]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433],\n        ID[\"EPSG\",8901]],\n    CS[ellipsoidal,2],\n        AXIS[\"longitude\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]],\n        AXIS[\"latitude\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]]")
# 
# # epsg6703km <- paste(
# #   "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
# #   "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
# #   "+units=km +no_defs"
# # )
# 
# locations<- endo_herb_data %>% 
#   st_as_sf(coords = c("lon", "lat"), crs = crs, remove = FALSE)  %>% 
#   st_transform(crs) %>% 
#   dplyr::select(lon,lat, Sample_id, Institution_specimen_id, new_id, Spp_code, year, month, day)
# 
# # reading in tif files that have 250m resolution rasters of nitrogen surplus components
# dir <- "./Analyses/Supp_Nitrogen_Data_Sources/gTREND-nitrogen/"
# gTrends_atm_ox_files <- list.files(path = paste0(dir,"Atmospheric_Oxidized/" ), pattern = ".tif$")
# gTrends_atm_red_files <- list.files(path = paste0(dir,"Atmospheric_Reduced/" ), pattern = ".tif$")
# gTrends_nSurplus_files <- list.files(path = paste0(dir,"Surplus/" ), pattern = ".tif$")
# 
# # load the raster files as a stack 
# # this takes alittle bit, especially if using the full set of years
# 
# atm_ox_stack <- terra::rast(stack(paste0(dir,"Atmospheric_Oxidized/", gTrends_atm_ox_files))) 
# atm_red_stack <- terra::rast(stack(paste0(dir,"Atmospheric_Reduced/", gTrends_atm_red_files))) 
# nSurplus_stack <- terra::rast(stack(paste0(dir,"Surplus/", gTrends_nSurplus_files)))
# 
# 
# 
# 
# 
# 
# # extracting the annual values at our coordinates
# coords_df <- locations %>% 
#   dplyr::select(lon,lat) %>% 
#   distinct() 
# 
# coords <- SpatialPoints(cbind(coords_df$lon, coords_df$lat), proj4string = CRS(crs))
# buffers_10 <- st_buffer(st_as_sf(coords), dist = 10000) # 10 km buffer
# # buffers_30 <- st_buffer(st_as_sf(coords), dist = 30000) # 30 km buffer
# 
# # extract each monthly measurement as the mean of values within each buffer (Takes a long)
# terra::gdalCache(3276) # need to expand the "gdal cache" when using the full raster stack.
# 
# atm_ox_10km <- exactextractr::exact_extract(atm_ox_stack, buffers_10, fun = "mean")
# saveRDS(atm_ox_10km, paste0(dir,"gTREND-nitrogen_extracted/atm_ox_10km.Rds")
# 
# 
# atm_red_10km <- exactextractr::exact_extract(atm_red_stack, buffers_10, fun = "mean")
# saveRDS(atm_red_10km, "gTREND-nitrogen_extracted/atm_red_10km.Rds")
# 
# nSurplus_10km <- exactextractr::exact_extract(nSurplus_stack, buffers_10, fun = "mean")
# saveRDS(nSurplus_10km, "gTREND-nitrogen_extracted/nSurplus_10km.Rds")
# 





# Reading in the HISLAND-US data for the locations in our study

######################################################
##### Read in the land cover data for collection locations ###########
######################################################

# reading in tif files that have 250m resolution rasters of nitrogen surplus components
dir <- "./Analyses/SUPP_LULC_Data_Sources/"
lulc_files <- list.files(path = paste0(dir,"conus_lulc_boolean/" ), pattern = ".tif$")

# Looking only at rasters after 1895, but actual supplement will use 1930 to match nitrogen data
lulc_files <- lulc_files[(1895-1629):(2020-1629)]
# load the raster files as a stack 
# this takes alittle bit, especially if using the full set of years

lulc_stack <- terra::rast(stack(paste0(dir,"conus_lulc_boolean/", lulc_files)))



endo_herb_data <- read_csv(file = "/Users/joshuacfowler/Documents/R_projects/Endophyte_herbarium_urbanization/Analyses/endo_herb_nit.csv")



# converting the lat long to same crs as the rasters are stored in
# define a crs
# crs <- paste("PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"unknown\",\n        DATUM[\"North American Datum 1983\",\n            ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n                LENGTHUNIT[\"metre\",1]],\n            ID[\"EPSG\",6269]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8901]]],\n    CONVERSION[\"unknown\",\n        METHOD[\"Albers Equal Area\",\n            ID[\"EPSG\",9822]],\n        PARAMETER[\"Latitude of false origin\",23,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8821]],\n        PARAMETER[\"Longitude of false origin\",-96,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8822]],\n        PARAMETER[\"Latitude of 1st standard parallel\",29.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8823]],\n        PARAMETER[\"Latitude of 2nd standard parallel\",45.5,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8824]],\n        PARAMETER[\"Easting at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8826]],\n        PARAMETER[\"Northing at false origin\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8827]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1,\n                ID[\"EPSG\",9001]]]]")
# crs <- paste("GEOGCRS[\"unknown\",\n    DATUM[\"North American Datum 1983\",\n        ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n            LENGTHUNIT[\"metre\",1]],\n        ID[\"EPSG\",6269]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433],\n        ID[\"EPSG\",8901]],\n    CS[ellipsoidal,2],\n        AXIS[\"longitude\",east,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]],\n        AXIS[\"latitude\",north,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]]")

crs <- crs(lulc_stack[[1]])
locations<- endo_herb_data %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)  %>% 
  st_transform(crs) %>%
  dplyr::select(lon,lat, Sample_id, Institution_specimen_id, new_id, Spp_code, year, month, day)

# extracting the annual values at our coordinates
coords_df <- locations %>%
  dplyr::select(lon,lat) %>%
  distinct()

buffers_10 <- st_buffer(coords_df, dist = 10000) # 10 km buffer
# buffers_30 <- st_buffer(st_as_sf(coords), dist = 30000) # 30 km buffer

# extract each monthly measurement as the mean of values within each buffer (Takes a long)
terra::gdalCache(3276) # need to expand the "gdal cache" when using the full raster stack.
# target_val <- 6
lulc_10km <- exactextractr::exact_extract(lulc_stack, buffers_10, fun =  "frac")


# LULC type:
# 0 nodata value
# 1 urban
# 2 crop
# 3 pasture
# 4 forest
# 5 shrub
# 6 grassland
# 7 wetland
# 8 water
# 9 barren

lulc_10km_pivot <- lulc_10km %>% 
  bind_cols(coords_df) %>% 
  pivot_longer(cols = -c(lon, lat, geometry), names_to = c("type", "year"), names_sep = "\\.", values_to = "fraction") %>% 
  mutate(across(c(year), ~ parse_number(.x))) %>% 
  pivot_wider(names_from = type, values_from = fraction) %>% 
  mutate(lulc_PercentUrban = frac_1, 
         lulc_PercentAg = frac_2 + frac_3)
                  
ggplot(lulc_10km_pivot)+
  geom_point(aes(x = lulc_PercentAg, y =lulc_PercentUrban, color = year))

ggplot(lulc_10km_pivot)+
  geom_line(aes(x = year, y =lulc_PercentUrban, group = geometry), alpha = .1)

ggplot(lulc_10km_pivot)+
  geom_line(aes(x = year, y =lulc_PercentAg, group = geometry), alpha = .1)


write_csv(lulc_10km_pivot, paste0(dir,"lulc_extracted/lulc_10km_pivot.csv"))





################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################
# This is where I'm loading in the version of the data set without land cover and nitrogen data, but we ought to be able to replace this easily


Mallorypath <- "C:/Users/malpa/OneDrive/Documents/EndoHerbQGIS/"
Joshpath <- "Analyses/"
path <- Joshpath


endo_herb_georef <- read_csv(file = paste0(path, "full_Zonalhist_NLCD_2001_10km.csv")) %>%
  filter(Country != "Canada") %>%
  filter(Country != "CA") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>%
  mutate(species_index = as.factor(case_when(Spp_code == "AGHY" ~ "1",
                                             Spp_code == "AGPE" ~ "2",
                                             Spp_code == "ELVI" ~ "3"))) %>%
  mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
                             Spp_code == "AGPE" ~ "A. perennans",
                             Spp_code == "ELVI" ~ "E. virginicus")) %>%
  mutate(decade = floor(year/10)*10)%>%
  mutate(DevelopedOpenSpace = HISTO_21,
         DevelopedLowIntensity = HISTO_22,
         MediumDeveloped = HISTO_23,
         HighDeveloped = HISTO_24,
         PastureHay = HISTO_81,
         CultivatedCrops = HISTO_82,
         TotalAg = PastureHay + CultivatedCrops,
         TotalPixels = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24 + HISTO_0 + HISTO_11 +HISTO_12 + HISTO_31 +HISTO_41 + HISTO_42 + HISTO_43 + HISTO_52 + HISTO_71+ HISTO_90 + TotalAg + HISTO_95,
         TotalDeveloped = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24,
         OtherLC = (TotalPixels - (TotalAg + TotalDeveloped))/TotalPixels *100,
         PercentUrban = TotalDeveloped/TotalPixels * 100,
         PercentAg = TotalAg/TotalPixels * 100)

#fixing column names in yearly_endo and filtering
yearly_endo <- read_csv(file = paste0(path,"endo_herb_yearly_nlcd.csv"))%>%
  filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>%
  mutate(species_index = as.factor(case_when(Spp_code == "AGHY" ~ "1",
                                             Spp_code == "AGPE" ~ "2",
                                             Spp_code == "ELVI" ~ "3"))) %>%
  mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
                             Spp_code == "AGPE" ~ "A. perennans",
                             Spp_code == "ELVI" ~ "E. virginicus")) %>%
  mutate(decade = floor(year/10)*10)%>%
  mutate(DevelopedOpenSpace = HISTO_21,
         DevelopedLowIntensity = HISTO_22,
         MediumDeveloped = HISTO_23,
         HighDeveloped = HISTO_24,
         PastureHay = HISTO_81,
         CultivatedCrops = HISTO_82,
         TotalAg = PastureHay + CultivatedCrops,
         TotalPixels = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24 + HISTO_11 +HISTO_12 + HISTO_31 +HISTO_41 + HISTO_42 + HISTO_43 + HISTO_52 + HISTO_71+ HISTO_90 + TotalAg + HISTO_95,
         TotalDeveloped = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24,
         OtherLC = (TotalPixels - (TotalAg + TotalDeveloped))/TotalPixels *100,
         spec_PercentUrban = TotalDeveloped/TotalPixels * 100,
         spec_PercentAg = TotalAg/TotalPixels * 100)%>%
  dplyr::select(Sample_id,spec_PercentUrban, spec_PercentAg)



endo_herb_georef1 <- left_join(endo_herb_georef, yearly_endo, by = "Sample_id")


# Doing some filtering to remove NA's and some data points that probably aren't accurate species id's
endo_herb_merge1 <- endo_herb_georef1 %>%
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>%
  filter(!is.na(lon) & !is.na(year)) %>%
  filter(lon>-110 ) %>%
  filter(Country != "Canada" ) %>%
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>%
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+")) 
#loading in nitrogen data too
# nit <- read.csv(file = "endo_herb_nit.csv") %>%
#   select(Sample_id, NO3_mean, NH4_mean, TIN_mean)

nit_avgs <- read.csv(file = paste0(path, "nitrogen_mean_df.csv")) %>% 
  dplyr::select(lon, lat, buffer, mean_TIN, mean_NO3, mean_NH4) %>% 
  pivot_wider(id_cols = c(lon, lat), names_from = c("buffer"), values_from = c("mean_TIN", "mean_NO3", "mean_NH4"))
endo_herb_merge2 <- left_join(endo_herb_merge1, nit_avgs, by = c("lon", "lat"))


nit_yearly <- read.csv(file = paste0(path,"nitrogen_yearly_df.csv")) %>% 
  dplyr::select(lon, lat, year, buffer, TIN, NO3, NH4) %>% 
  pivot_wider(id_cols = c(lon, lat, year), names_from = c("buffer"), values_from = c("TIN", "NO3", "NH4"))

endo_herb <- left_join(endo_herb_merge2, nit_yearly, by = c("lon", "lat", "year")) %>% 
  mutate(sample_temp = Sample_id) %>%
  separate(sample_temp, into = c("Herb_code", "spp_code", "specimen_code", "tissue_code")) %>%
  mutate(species_index = as.factor(case_when(spp_code == "AGHY" ~ "1",
                                             spp_code == "AGPE" ~ "2",
                                             spp_code == "ELVI" ~ "3"))) %>%
  mutate(species = case_when(spp_code == "AGHY" ~ "A. hyemalis",
                             spp_code == "AGPE" ~ "A. perennans",
                             spp_code == "ELVI" ~ "E. virginicus")) %>%
  filter(scorer_id != "Scorer26") %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(spp_code)) %>%
  filter(!is.na(lon) & !is.na(year)) %>%
  filter(!is.na(PercentAg), !is.na(mean_NO3_10km)) 




# Creating scorer and collector levels
scorer_levels <- levels(as.factor(endo_herb$scorer_id))
scorer_no <- paste0("Scorer",1:nlevels(as.factor(endo_herb$scorer_id)))

endo_herb$scorer_factor <- scorer_no[match(as.factor(endo_herb$scorer_id), scorer_levels)]


collector_levels <- levels(as.factor(endo_herb$collector_string))
collector_no <- paste0("Collector",1:nlevels(as.factor(endo_herb$collector_string)))

endo_herb$collector_factor <- collector_no[match(as.factor(endo_herb$collector_string), collector_levels)]

# converting the lat long to epsg 6703km in km
# define a crs
epsg6703km <- paste(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5",
  "+lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83",
  "+units=km +no_defs"
)

endo_herb_sf<- endo_herb %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) %>% 
  mutate(scorer_index = parse_number(scorer_factor),
         collector_index = parse_number(collector_factor)) 



climate <- read_csv(file = paste0(path,"PRISM_yearly_df.csv")) %>% 
  pivot_wider(id_cols = c(lon, lat, year), names_from = c("buffer"), values_from = c("tmean", "ppt"))

endo_herb <- left_join(endo_herb_sf, climate, by = c("lon", "lat", "year")) %>% 
  filter(year>=1895) %>% 
  filter(seed_scored != 0) %>% 
  filter(Sample_id != "BRIT_AGHY_506") %>% # dropping this because the georeferenceing is wrong
  filter(!(Sample_id == "AM_ELVI_122" & is.na(Municipality))) %>% 
  dplyr::distinct(Sample_id, score_number, .keep_all = TRUE)



####################################################################################
###### merging endophyte scoring dataset with the TREND-Nitrogen data ##############
####################################################################################
county_in_crs <- st_transform(county_shape, crs = crs(endo_herb))


atm_ox_join <- left_join(county_in_crs, atm_ox) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "TREND_NOx", names_to = "year") %>% mutate(year = parse_number(year), TREND_NOx = 100*TREND_NOx)
atm_red_join <- left_join(county_in_crs, atm_red) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "TREND_NHx", names_to = "year") %>% mutate(year = parse_number(year), TREND_NHx = 100*TREND_NHx)
nit_surplus_join <- left_join(county_in_crs, nit_surplus) %>% left_join(state_fips) %>% 
  pivot_longer(cols = starts_with("y"), values_to = "TREND_NSurplus", names_to = "year") %>% mutate(year = parse_number(year), TREND_NSurplus = 100*TREND_NSurplus)



endo_herb_TREND <- endo_herb %>% st_join(atm_ox_join) %>%  # join = join_by("year" == "year", "County_fixed" == "name"))
  filter(year.x == year.y) %>% rename(year = year.x) %>% dplyr::select(-year.y) %>% 
  st_join(atm_red_join) %>% 
  filter(year.x == year.y) %>% rename(year = year.x) %>% dplyr::select(-year.y) %>% 
  st_join(nit_surplus_join) %>% 
  filter(year.x == year.y) %>% rename(year = year.x) %>% dplyr::select(-year.y) %>% 
  mutate(TREND_NDep = TREND_NOx + TREND_NHx) %>% 
  left_join(NDep_change)
  
  
endo_herb_check <- endo_herb_TREND %>% dplyr::select(year, mean_TIN_10km, TREND_NOx, TREND_NHx,TREND_NDep, TREND_NSurplus, County, County_fixed, name, State, state)


ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = TREND_NOx, color = as.numeric(as.factor(name)))) + guides(color = "none")
ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = TREND_NHx, color = as.numeric(as.factor(name)))) + guides(color = "none")
ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = TREND_NDep, color = as.numeric(as.factor(name)))) + guides(color = "none")
ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = TREND_NSurplus, color = as.numeric(as.factor(name)))) + guides(color = "none")


ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = NDep_change)) + guides(color = "none")

ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = NDep_Normal)) + guides(color = "none")
ggplot(endo_herb_TREND)+
  geom_point(aes(x = mean_TIN_10km, y = NOx_Normal)) + guides(color = "none")




####################################################################################
###### merging endophyte scoring dataset with the extracted lulc dataset ##############
####################################################################################

lulc_10km_pivot <- read_csv(file = "./Analyses/SUPP_LULC_Data_Sources/lulc_extracted/lulc_10km_pivot.csv") %>% 
  st_drop_geometry() %>% dplyr::select(-geometry)

endo_herb_TREND_lulc <- endo_herb_TREND %>% 
  left_join(lulc_10km_pivot) 
  # dplyr::select(Sample_id, County, State, lat, lon, year, mean_TIN_10km, NDep_Normal, frac_1, PercentUrban, lulc_PercentUrban )

ggplot(endo_herb_TREND_lulc)+
  geom_point(aes(x = PercentUrban, lulc_PercentUrban))
ggplot(endo_herb_TREND_lulc)+
  geom_point(aes(x = PercentAg, lulc_PercentAg))
##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
data_summary <- endo_herb_TREND_lulc %>% 
  dplyr::summarize(year = mean(year, na.rm = T),
                   lulc_PercentUrban = mean(lulc_PercentUrban, na.rm = T),
                   lulc_PercentAg = mean(lulc_PercentAg, na.rm = T),
                   TREND_NOx= mean(TREND_NOx, na.rm = T),
                   TREND_NHx= mean(TREND_NHx, na.rm = T),
                   TREND_NDep= mean(TREND_NDep, na.rm = T),
                   TREND_NSurplus= mean(TREND_NSurplus, na.rm = T),
                   tmean_10km = mean(tmean_10km, na.rm = T),
                   ppt_10km = mean(ppt_10km, na.rm = T))
data <- endo_herb_TREND_lulc %>% 
  mutate(year = year - data_summary$year,
         lulc_PercentUrban = lulc_PercentUrban - data_summary$lulc_PercentUrban,
         lulc_PercentAg = lulc_PercentAg - data_summary$lulc_PercentAg,
         TREND_NOx = TREND_NOx - data_summary$TREND_NOx,
         TREND_NHx = TREND_NHx - data_summary$TREND_NHx,
         TREND_NDep = TREND_NDep - data_summary$TREND_NDep,
         TREND_NSurplus = TREND_NSurplus - data_summary$TREND_NSurplus,
         tmean_10km = tmean_10km - data_summary$tmean_10km,
         ppt_10km = ppt_10km - data_summary$ppt_10km) %>%  
  mutate(Spp_index = as.numeric(as.factor(Spp_code))) %>% 
  filter(!is.na(ppt_10km))
# spp_count <- data %>% 
#   filter(score_number == 1) %>% 
#   group_by(species) %>% summarize(n())
#  county_count <- data %>% 
#    filter(score_number == 1, is.na(Municipality) & !is.na(County)) %>% 
#    summarize(n())
# city_count <- data %>% 
#    filter(score_number == 1, !is.na(Municipality) & !is.na(County)) %>% 
#    summarize(n())

# endo_status_summary <- data %>%
#   group_by(species) %>% 
#   summarize(mean = mean(Endo_status_liberal),
#             sd = sd(Endo_status_liberal)^2)


# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
coords <- cbind(data$easting, data$northing)

non_convex_bdry <- fm_extensions(
  data$geometry,
  convex = c(250, 500),
  concave = c(250, 500),
  crs = fm_crs(data)
)

coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE))) %>% st_transform(epsg6703km)
# plot(coastline)




bdry <- st_intersection(coastline$geom, non_convex_bdry[[1]])

# plot(bdry)

bdry_polygon <- st_cast(st_zm(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")

non_convex_bdry[[1]] <- bdry_polygon



max.edge = diff(range(coords[,1]))/(100)


mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge*2, max.edge*8), # km inside and outside
  cutoff = max.edge,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge
# plot it
# plot(mesh)


mesh_plot <- ggplot() +
  gg(data = mesh) +
  geom_point(data = data, aes(x = easting, y = northing, col = species), size = .8) +
  coord_sf()+
  theme_bw() +
  labs(x = "", y = "", color = "Species")+
  theme(legend.text = element_text(face = "italic"))
# mesh_plot
# ggsave(mesh_plot, filename = "Plots/mesh_plot.png", width = 6, height = 5)



# make spde (stochastic partial differential equation)

# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- 1

# The priors from online tutorials are :   # P(practic.range < 0.05) = 0.01 # P(sigma > 1) = 0.01
# For ESA presentation, I used the following which at least "converged" but seem sensitive to choices
# for AGHY =  P(practic.range < 0.1) = 0.01 # P(sigma > 1) = 0.01
# for AGPE =  P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
# for ELVI = P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
spde <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.5),
  
  prior.sigma = c(prior_sigma, 0.5)
)

# inlabru makes making spatial effects simpler compared to "INLA" because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# setting the random effects prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))

# This is the model formula with a spatial effect (spatially varying intercept). To this, we can add predictor variables

# formula
# version for each species separately
# s_components <- ~ Intercept(1) +
#   year +
#   space_int(coords, model = spde)

# View(model.matrix(~0 + Spp_code*PercentAg, data))


# version with all species in one model. Note that we remove the intercept, and then we have to specify that the species is a factor 

# comparing different levels of interactions




# s_components.NOx <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(TREND_NOx + PercentAg + PercentUrban + ppt_10km + tmean_10km), model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde)

# s_components.NSurplus <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(TREND_NSurplus + PercentAg + PercentUrban + ppt_10km + tmean_10km), model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde)

s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(TREND_NDep + lulc_PercentAg + lulc_PercentUrban + ppt_10km + tmean_10km), model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.year <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(TREND_NDep*year + lulc_PercentAg*year + lulc_PercentUrban*year + ppt_10km  + tmean_10km), model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



s_formula <- Endo_status_liberal ~ .


# Now run the model



fit <- bru(s_components,
           like(
             formula = s_formula,
             family = "binomial",
             Ntrials = 1,
             data = data
           ),
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE
           )
)



fit.year <- bru(s_components.year,
                like(
                  formula = s_formula,
                  family = "binomial",
                  Ntrials = 1,
                  data = data
                ),
                options = list(
                  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
                  control.inla = list(int.strategy = "eb"),
                  verbose = TRUE
                )
)


# DIC
fit$dic$dic
fit.year$dic$dic

# WAIC
fit$waic$waic
fit.year$waic$waic

# cpo
-mean(log(fit$cpo$cpo))
-mean(log(fit.year$cpo$cpo))

# checking convergence
fit$mode$mode.status
fit.year$mode$mode.status



################################################################################################################################
##########  Plotting the prediction without year effects for NDep ###############
################################################################################################################################

min_ag<- min(data$lulc_PercentAg)
max_ag <- max(data$lulc_PercentAg)

min_urb<- min(data$lulc_PercentUrban)
max_urb<- max(data$lulc_PercentUrban)

min_nit<- min(data$TREND_NDep)
max_nit<- max(data$TREND_NDep)

preddata.1 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     lulc_PercentAg = rep(seq(min_ag, max_ag, length.out = 50), times = 3),
                     lulc_PercentUrban = 0,
                     TREND_NDep = 0,
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))
preddata.2 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     lulc_PercentAg = 0,
                     lulc_PercentUrban = rep(seq(min_urb, max_urb, length.out = 50), times = 3),
                     TREND_NDep = 0,
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

preddata.3 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     lulc_PercentAg = 0,
                     lulc_PercentUrban = 0,
                     TREND_NDep = rep(seq(min_nit, max_nit, length.out = 50), times = 3),
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

ag.pred <- predict(
  fit,
  newdata = preddata.1,
  formula = ~ invlogit(fixed),#  + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(lulc_PercentAg = lulc_PercentAg + data_summary$lulc_PercentAg) # undoing mean centering

urb.pred <- predict(
  fit,
  newdata = preddata.2,
  formula = ~ invlogit(fixed),#  + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(lulc_PercentUrban = lulc_PercentUrban + data_summary$lulc_PercentUrban) # undoing mean centering



nit.pred <- predict(
  fit,
  newdata = preddata.3,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(TREND_NDep = TREND_NDep + data_summary$TREND_NDep) # undoing mean centering



values <-  c("#b2abd2", "#5e3c99")


ag_binned <- endo_herb_TREND_lulc %>% 
  mutate(ag_bin = cut(lulc_PercentAg, breaks = 30)) %>% 
  group_by(species, ag_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_ag = mean(lulc_PercentAg),
            sample = n())

urb_binned <- endo_herb_TREND_lulc %>% 
  mutate(urb_bin = cut(lulc_PercentUrban, breaks = 30)) %>% 
  group_by(species, urb_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_urb = mean(lulc_PercentUrban),
            sample = n())

nit_binned <- endo_herb_TREND_lulc %>% 
  mutate(nit_bin = cut(TREND_NDep, breaks = 30)) %>% 
  group_by(species, nit_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_nit = mean(TREND_NDep),
            sample = n())

ag_trend <- ggplot(ag.pred) +
  geom_line(aes(x = lulc_PercentAg, mean)) +
  geom_ribbon(aes(lulc_PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#B38600") +
  geom_ribbon(aes(lulc_PercentAg, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = ag_binned, aes(x = mean_ag, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,400))+
  facet_wrap(~species,  ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Ag. (%)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))#+
# lims(y = c(0,1), x = c(0, 100))

urb_trend <- ggplot(urb.pred) +
  geom_line(aes(lulc_PercentUrban, mean)) +
  geom_ribbon(aes(lulc_PercentUrban, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#021475") +
  geom_ribbon(aes(lulc_PercentUrban, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = urb_binned, aes(x = mean_urb, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,400))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))


nit_trend <- ggplot(nit.pred) +
  geom_line(aes(TREND_NDep, mean)) +
  geom_ribbon(aes(TREND_NDep, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#BF00A0") +
  geom_ribbon(aes(TREND_NDep, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = nit_binned, aes(x = mean_nit, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,400))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Est. of Historic Nit. Dep", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1.1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))



ag_trend <- tag_facet(ag_trend)
urb_trend <- tag_facet(urb_trend, tag_pool =  letters[-(1:3)])
nit_trend <- tag_facet(nit_trend, tag_pool =  letters[-(1:6)])
fig2 <-   ag_trend + urb_trend + nit_trend + plot_layout(ncol = 3, guides = "collect") 
ggsave(fig2, file = "Plots/Figure_2_temporal_data_sources_supplement.png", width = 10, height = 8)




################################################################################################################################
##########  Plotting the posteriors from the model without year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
param_names <- fit$summary.random$fixed$ID

n_draws <- 1000

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors <- generate(
  fit,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))


posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")

colnames(posteriors_df) <- sub("Spp_code", "", colnames(posteriors_df))
colnames(posteriors_df) <- gsub(":", ".", colnames(posteriors_df))

# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "lulc_PercentAg" = "Agr.",
                                                         "lulc_PercentUrban" = "Urb.",
                                                         "TREND_NDep" = "Nit.",
                                                         "ppt_10km" = "PPT.",
                                                         "tmean_10km" = "Temp.")),
                          levels = c("Intercept","Nit.","Agr.","Urb.","PPT.","Temp.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))







posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", ncol = 6)+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/posterior_hist_temporal_data_sources_Supp.png", width = 8, height = 3)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            median = median(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)
# write.csv(effects_summary, file = "Posterior_prob_results.csv")




################################################################################################################################
##########  Plotting the prediction without year effects for NOx ###############
################################################################################################################################

min_ag<- min(data$PercentAg)
max_ag <- max(data$PercentAg)

min_urb<- min(data$PercentUrban)
max_urb<- max(data$PercentUrban)

min_nit<- min(data$TREND_NOx)
max_nit<- max(data$TREND_NOx)

preddata.1 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = rep(seq(min_ag, max_ag, length.out = 50), times = 3),
                     PercentUrban = 0,
                     TREND_NOx = 0,
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))
preddata.2 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = 0,
                     PercentUrban = rep(seq(min_urb, max_urb, length.out = 50), times = 3),
                     TREND_NOx = 0,
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

preddata.3 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = 0,
                     PercentUrban = 0,
                     TREND_NOx = rep(seq(min_nit, max_nit, length.out = 50), times = 3),
                     ppt_10km = 0,
                     tmean_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

ag.pred <- predict(
  fit.NOx,
  newdata = preddata.1,
  formula = ~ invlogit(fixed),#  + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(PercentAg = PercentAg + data_summary$PercentAg) # undoing mean centering

urb.pred <- predict(
  fit.NOx,
  newdata = preddata.2,
  formula = ~ invlogit(fixed),#  + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(PercentUrban = PercentUrban + data_summary$PercentUrban) # undoing mean centering



nit.pred <- predict(
  fit.NOx,
  newdata = preddata.3,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(TREND_NOx = TREND_NOx + data_summary$TREND_NOx) # undoing mean centering



values <-  c("#b2abd2", "#5e3c99")


ag_binned <- endo_herb_TREND %>% 
  mutate(ag_bin = cut(PercentAg, breaks = 30)) %>% 
  group_by(species, ag_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_ag = mean(PercentAg),
            sample = n())

urb_binned <- endo_herb_TREND %>% 
  mutate(urb_bin = cut(PercentUrban, breaks = 30)) %>% 
  group_by(species, urb_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_urb = mean(PercentUrban),
            sample = n())

nit_binned <- endo_herb_TREND %>% 
  mutate(nit_bin = cut(TREND_NOx, breaks = 30)) %>% 
  group_by(species, nit_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_nit = mean(TREND_NOx),
            sample = n())

ag_trend <- ggplot(ag.pred) +
  geom_line(aes(x = PercentAg, mean)) +
  geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#B38600") +
  geom_ribbon(aes(PercentAg, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = ag_binned, aes(x = mean_ag, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,310))+
  facet_wrap(~species,  ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Ag. (%)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))#+
# lims(y = c(0,1), x = c(0, 100))

urb_trend <- ggplot(urb.pred) +
  geom_line(aes(PercentUrban, mean)) +
  geom_ribbon(aes(PercentUrban, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#021475") +
  geom_ribbon(aes(PercentUrban, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = urb_binned, aes(x = mean_urb, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,310))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))


nit_trend <- ggplot(nit.pred) +
  geom_line(aes(TREND_NOx, mean)) +
  geom_ribbon(aes(TREND_NOx, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#BF00A0") +
  geom_ribbon(aes(TREND_NOx, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = nit_binned, aes(x = mean_nit, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,310))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Est. of Historic Nit. Dep. (kg N/km^2/year)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1.1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))



ag_trend <- tag_facet(ag_trend)
urb_trend <- tag_facet(urb_trend, tag_pool =  letters[-(1:3)])
nit_trend <- tag_facet(nit_trend, tag_pool =  letters[-(1:6)])
fig2 <-   ag_trend + urb_trend + nit_trend + plot_layout(ncol = 3, guides = "collect") 
ggsave(fig2, file = "Plots/Figure_2_TREND_Nitrogen.png", width = 10, height = 8)



# 
# ################################################################################################################################
# ##########  Plotting the posteriors from the model without year effect for NOx###############
# ################################################################################################################################
# 
# # param_names <- fit.4$summary.random$fixed$ID
# param_names <- fit.NOx$summary.random$fixed$ID
# 
# n_draws <- 1000
# 
# # we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
# posteriors <- generate(
#   fit.NOx,
#   formula = ~ fixed_latent,
#   n.samples = n_draws) 
# rownames(posteriors) <- param_names
# colnames(posteriors) <- c( paste0("iter",1:n_draws))
# 
# 
# posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")
# 
# colnames(posteriors_df) <- sub("Spp_code", "", colnames(posteriors_df))
# colnames(posteriors_df) <- gsub(":", ".", colnames(posteriors_df))
# 
# # Calculate the effects of the predictor, given that the reference level is for AGHY
# effects_df <- posteriors_df %>% 
#   rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
#   pivot_longer( cols = -c(iteration), names_to = "param") %>% 
#   mutate(model = "No Year") %>% 
#   mutate(param_label = sub("^[^.]+.", "", param),
#          spp_label = sub("\\..*","", param)) %>% 
#   mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
#                                                          "Int" = "Intercept",
#                                                          "PercentAg" = "Agr.",
#                                                          "PercentUrban" = "Urb.",
#                                                          "TREND_NOx" = "Nit.",
#                                                          "ppt_10km" = "PPT.",
#                                                          "tmean_10km" = "Temp.")),
#                           levels = c("Intercept","Nit.","Agr.","Urb.","PPT.","Temp.")),
#          spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
#                         levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))
# 
# 
# 
# 
# 
# 
# 
# posterior_hist <- ggplot(effects_df)+
#   stat_halfeye(aes(x = value, y = spp_f, fill  = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
#   
#   # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
#   # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
#   # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
#   # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
#   
#   geom_vline(xintercept = 0)+
#   facet_wrap(~param_f, scales = "free_x", ncol = 6)+
#   labs(x = "Posterior Est.", y = "Species")+
#   guides(fill = "none")+
#   scale_color_manual(values = species_colors)+
#   scale_fill_manual(values = species_colors)+
#   scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
#   theme_bw() + theme(axis.text.y = element_text(face = "italic"),
#                      axis.text.x = element_text(size = rel(.8)))
# 
# # posterior_hist
# ggsave(posterior_hist, filename = "Plots/posterior_hist.png", width = 8, height = 3)
# 
# 
# 
# 
# effects_summary <- effects_df %>% 
#   group_by(param, param_label, spp_label) %>% 
#   summarize(mean = mean(value), 
#             median = median(value), 
#             lwr = quantile(value, .025),
#             upr = quantile(value, .975),
#             prob_pos = sum(value>0)/1000,
#             prob_neg = 1-prob_pos)
# # write.csv(effects_summary, file = "Posterior_prob_results.csv")
# 
# 

################################################################################################################################
##########  Plotting the posteriors from the model with year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
param_names <- fit.year$summary.random$fixed$ID
param_names <- gsub(":year:lulc_PercentAg", ":lulc_PercentAg:year", param_names)
param_names <- gsub(":year:lulc_PercentUrban", ":lulc_PercentUrban:year", param_names)

n_draws <- 1000

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors <- generate(
  fit.year,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))


posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")


colnames(posteriors_df) <- sub("Spp_code", "", colnames(posteriors_df))
colnames(posteriors_df) <- gsub(":", ".", colnames(posteriors_df))

# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "Year") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "year" = "Year",
                                                         "lulc_PercentAg" = "Agr.",
                                                         "lulc_PercentUrban" = "Urb.",
                                                         "TREND_NDep" = "Nit.",
                                                         "ppt_10km" = "Ppt.",
                                                         "tmean_10km" = "Temp.")),
                          levels = c("Intercept"  ,"Nit.","Agr.","Urb.","Ppt.", "Temp.", "Year", "Agr. X Year", "Urb. X Year", "Nit. X Year")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))




posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", ncol = 6)+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/posterior_hist_TREND_lulc_yearmodel.png", width = 9, height = 6)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            median = median(value),
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)



###############################################################################################################################
##########  simpler plot of temporal trends ###############
################################################################################################################################
# note that this model now has dynamic values of anthropogenic covariates and a year effect
# The year effect is essentially capturing how much temporal change occured in prevalence beyond that driven by the underlying shift in each covariate. 
# This is not the same as in the main analysis, where the equivalent figure asks how regions with high vs low anthropogenic covariate had relatively higher/lower rates of change

#### Generating plot showing predicted prevalence in the past (1930) and in the present (2017) given the covariates at that time to assess the overall change over time

lulc_start_end <- lulc_10km_pivot %>% 
  filter(year %in% c(1930, 2017)) %>% 
  pivot_wider(id_cols = c(lon, lat), names_from = year, values_from = c(lulc_PercentUrban, lulc_PercentAg))
  
NDep_start_end <- NDep_change %>% dplyr::select("STATEFP", "COUNTYFP", "COUNTYNS", "AFFGEOID",   
                                                "GEOID", "NAME", "LSAD", "StateFIPS",
                                                "CountyFIPS", "AREA..HA.", "name", "state",   
                                                "NDep_1930", "NDep_2017")
# We can use NDep_change for the same info about nit dep, and this is already merged into the original data dataframe
data_change <- data %>% left_join(lulc_start_end) 

preddata_1930 <- data_change %>% 
  dplyr::select(Spp_code, year,lon, lat, TREND_NDep, lulc_PercentAg, lulc_PercentAg, ppt_10km, tmean_10km, NDep_1930, NDep_2017, lulc_PercentAg_1930, lulc_PercentAg_2017, lulc_PercentUrban_1930, lulc_PercentUrban_2017) %>% 
  st_drop_geometry() %>% 
  mutate(TREND_NDep = NDep_1930 - data_summary$TREND_NDep,
         lulc_PercentAg = lulc_PercentAg_1930 - data_summary$lulc_PercentAg,
         lulc_PercentUrban = lulc_PercentUrban_1930 - data_summary$lulc_PercentUrban,
         year = 1930 - data_summary$year,
         ppt_10km = 0, 
         tmean_10km = 0,
         collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         year_label = case_when(year >=0 ~ "max_year",
                                year < 0 ~ "min_year")) %>% dplyr::select(-c(NDep_1930, NDep_2017, lulc_PercentAg_1930, lulc_PercentAg_2017, lulc_PercentUrban_1930, lulc_PercentUrban_2017))

preddata_2017 <- data_change %>% 
  dplyr::select(Spp_code, year, lon, lat, TREND_NDep, lulc_PercentAg, lulc_PercentAg, ppt_10km, tmean_10km, NDep_1930, NDep_2017, lulc_PercentAg_1930, lulc_PercentAg_2017, lulc_PercentUrban_1930, lulc_PercentUrban_2017) %>% 
  st_drop_geometry() %>% 
  mutate(TREND_NDep = NDep_2017 - data_summary$TREND_NDep,
         lulc_PercentAg = lulc_PercentAg_2017 - data_summary$lulc_PercentAg,
         lulc_PercentUrban = lulc_PercentUrban_2017 - data_summary$lulc_PercentUrban,
         year = 2017 - data_summary$year,
         ppt_10km = 0, 
         tmean_10km = 0,
         collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         year_label = case_when(year >=0 ~ "max_year",
                                year < 0 ~ "min_year")) %>% dplyr::select(-c(NDep_1930, NDep_2017, lulc_PercentAg_1930, lulc_PercentAg_2017, lulc_PercentUrban_1930, lulc_PercentUrban_2017))

preddata <- bind_rows(preddata_1930, preddata_2017) %>% distinct()

# gennerating predictions and back-transforming the standardized year variable



year.pred <- generate(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  # probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 500) 
colnames(year.pred) <- paste0("iter", 1:500) 


avg_change <- tibble(preddata, as_tibble(year.pred)) %>% 
  pivot_longer(cols = iter1:iter500, names_to = "iteration", values_to = "posterior") %>% 
  dplyr::select(-year, -TREND_NDep, -lulc_PercentAg, -lulc_PercentUrban, -ppt_10km, -tmean_10km) %>% 
  pivot_wider(id_cols = c(Spp_code, lon, lat, species, iteration), names_from = c(year_label), values_from = posterior, names_prefix = "post.") %>% 
  mutate(diff = (post.max_year - post.min_year)*100) %>% 
  group_by(Spp_code, species, lon, lat) %>% 
  dplyr::summarise(diff_mean = mean(diff),
                   diff_median = median(diff),
                   lwr = quantile(diff, .025),
                   upr = quantile(diff, .975),
                   prob_pos = sum(diff>0)/500,
                   diff_prob = max((sum(diff<0)/500),(sum(diff>0)/500)),
                   diff_prob_threshold = case_when(diff_prob>=.90 ~ .8,#">=.90", 
                                                   TRUE ~ .2)) %>% left_join(lulc_start_end) %>% left_join(data %>% dplyr::select(lat, lon, NDep_1930, NDep_2017) %>% distinct() %>% st_drop_geometry()) 

##### univariate version of change plots #####

ag_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray60")+
  # geom_smooth(aes(x = lulc_PercentAg_2017-lulc_PercentAg_1930, y = diff_mean), method = "lm", color = "black")+
  geom_point(aes(x = lulc_PercentAg_2017-lulc_PercentAg_1930, y = diff_mean, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 3)+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Ag. Land Cover (%)", y = "Change in % Prevalence / Year")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_blank())
# ag_trend

urb_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray60")+
  # geom_smooth(aes(x = lulc_PercentUrban_2017-lulc_PercentUrban_1930, y = diff_mean), method = "lm", color = "black")+
  geom_point(aes(x = lulc_PercentUrban_2017-lulc_PercentUrban_1930, y = diff_mean, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 3)+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Urb. Land Cover (%)", y = "Change in % Prevalence / Year")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_blank())
# urb_trend

nit_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray60")+
  # geom_smooth(aes(x = NDep_2017-NDep_1930, y = diff_mean), method = "lm", color = "black")+
  geom_point(aes(x = NDep_2017-NDep_1930, y = diff_mean, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 3)+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Nit. Dep. (kg N/km^2)", y = "Change in % Prevalence / Year")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_text(face = "italic", angle = 0))
# nit_trend


ag_trend_tag <- tag_facet(ag_trend)
urb_trend_tag <- tag_facet(urb_trend, tag_pool =  letters[-(1:3)])
nit_trend_tag <- tag_facet(nit_trend, tag_pool =  letters[-(1:6)])

Fig_change_per_change_supp <- ag_trend_tag + urb_trend_tag + nit_trend_tag + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Rate of Change in Endophyte Prevalence")& theme(plot.title = element_text(size = rel(1.5)))
ggsave(Fig_change_per_change_supp, filename = "Plots/Fig_change_per_change_supp.png", width = 14, height = 10)

# 

###### Version with pairs plots of change #####
ag_nit_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = NDep_2017-NDep_1930, y = (lulc_PercentAg_2017-lulc_PercentAg_1930)*100, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 2,) +
  # geom_point(data = endo_herb, aes(x = PercentAg,y = mean_TIN_10km), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Nit. Dep. (kg N/km^2)", y = "1930 to 2017 Change in Ag. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_blank())
# ag_nit_trend

ag_nit_prob <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = NDep_2017-NDep_1930, y = (lulc_PercentAg_2017-lulc_PercentAg_1930)*100, fill = diff_prob), shape = 21, size = 2,  alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "1930 to 2017 Change in Nit. Dep. (kg N/km^2)", y = "1930 to 2017 Change in Ag. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_nit_prob

# ggplot(avg_change)+
#   geom_point(aes(x = lulc_PercentAg_2017, y = lulc_PercentAg_2017-lulc_PercentAg_1930))

ag_urb_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = (lulc_PercentAg_2017-lulc_PercentAg_1930)*100, y = (lulc_PercentUrban_2017-lulc_PercentUrban_1930)*100, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 2) +
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Ag. Land Cover (%)", y = "1930 to 2017 Change in Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_urb_trend

ag_urb_prob <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = (lulc_PercentAg_2017-lulc_PercentAg_1930)*100, y = (lulc_PercentUrban_2017-lulc_PercentUrban_1930)*100, fill = diff_prob), shape = 21, size = 2, alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "1930 to 2017 Change Ag. in Land Cover (%)", y = "1930 to 2017 Change in Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_urb_prob



urb_nit_trend <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = NDep_2017-NDep_1930, y = (lulc_PercentUrban_2017-lulc_PercentUrban_1930)*100, fill = diff_mean , alpha = diff_prob_threshold), shape = 21, size = 2) +
  # geom_point(data = endo_herb, aes(x = PercentAg,y = mean_TIN_10km), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "1930 to 2017 Change in Nitrogen Deposition (kg N/km^2)", y = "1930 to 2017 Change in Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_text(face = "italic", angle = 0))
# urb_nit_trend

urb_nit_prob <- ggplot(avg_change)+
  geom_hline(yintercept = 0, color = "gray75")+
  geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = NDep_2017-NDep_1930, y = (lulc_PercentUrban_2017-lulc_PercentUrban_1930)*100, fill = diff_prob), shape = 21, size = 2, alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "1930 to 2017 Change in Nitrogen Deposition (kg N/km^2)", y = "1930 to 2017 Change in Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_text(face = "italic", angle = 0))
# urb_nit_prob



ag_urb_trend_tag <- tag_facet(ag_urb_trend)
ag_nit_trend_tag <- tag_facet(ag_nit_trend, tag_pool =  letters[-(1:3)])
urb_nit_trend_tag <- tag_facet(urb_nit_trend, tag_pool =  letters[-(1:6)])

ag_urb_prob_tag <- tag_facet(ag_urb_prob, tag_pool =  letters[-(1:9)])
ag_nit_prob_tag <- tag_facet(ag_nit_prob, tag_pool =  letters[-(1:12)])
urb_nit_prob_tag <- tag_facet(urb_nit_prob, tag_pool =  letters[-(1:15)])



Fig_trends_suppA <- ag_urb_trend_tag + ag_nit_trend_tag + urb_nit_trend_tag + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Rate of Change in Endophyte Prevalence")& theme(plot.title = element_text(size = rel(1.5)))
Fig_trends_suppB<- ag_urb_prob_tag + ag_nit_prob_tag + urb_nit_prob_tag + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Posterior Probability of Effect")& theme(plot.title = element_text(size = rel(1.5)))

Fig_trends_supp <- wrap_elements(Fig_trends_suppA)/ wrap_elements(Fig_trends_suppB)

ggsave(Fig_trends_supp, filename = "Plots/Fig_trends_supp.png", width = 14, height = 14)




##### version plotted across contemporary deposition/land use ####



ag_nit_trend <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = NDep_2017, y = (lulc_PercentAg_2017)*100, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 3) +
  # geom_point(data = endo_herb, aes(x = PercentAg,y = mean_TIN_10km), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "Nit. Dep. (kg N/km^2)", y = "Ag. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_blank())
# ag_nit_trend

ag_nit_prob <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = NDep_2017, y = (lulc_PercentAg_2017)*100, fill = diff_prob), shape = 21, size = 3,  alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "Nit. Dep. (kg N/km^2)", y = "Ag. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_nit_prob

# ggplot(avg_change)+
#   geom_point(aes(x = lulc_PercentAg_2017, y = lulc_PercentAg_2017-lulc_PercentAg_1930))

ag_urb_trend <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = (lulc_PercentAg_2017)*100, y = (lulc_PercentUrban_2017)*100, fill = diff_mean, alpha = diff_prob_threshold), shape = 21, size = 3) +
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "Ag. Land Cover (%)", y = "Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_urb_trend

ag_urb_prob <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = (lulc_PercentAg_2017)*100, y = (lulc_PercentUrban_2017)*100, fill = diff_prob), shape = 21, size = 3, alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "1930 to 2027 Change Ag. in Land Cover (%)", y = "Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# ag_urb_prob



urb_nit_trend <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes(x = NDep_2017, y = (lulc_PercentUrban_2017)*100, fill = diff_mean , alpha = diff_prob_threshold), shape = 21, size = 3) +
  # geom_point(data = endo_herb, aes(x = PercentAg,y = mean_TIN_10km), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-60, 60))+
  # scale_shape_manual(values= c(16,21))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "% Prevalence / Period", x = "Nitrogen Deposition (kg N/km^2)", y = "Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right =element_text(face = "italic", angle = 0))
# urb_nit_trend

urb_nit_prob <- ggplot(avg_change)+
  # geom_hline(yintercept = 0, color = "gray75")+
  # geom_vline(xintercept = 0, color = "gray75")+
  geom_point(aes( x = NDep_2017, y = (lulc_PercentUrban_2017)*100, fill = diff_prob), shape = 21, size = 3, alpha = .8) +
  # stat_contour(aes(x = mean_TIN_10km, y = PercentAg, z = diff_prob, alpha = ..level..^4),position = "identity", breaks = seq(.5:1, by = .05), linewidth = .5, color = "black")+
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species, ncol = 1,  strip.position = "right", scales = "free")+
  scale_alpha_continuous(limits = c(0,1), guide = "none")+
  scale_fill_distiller(palette = "YlGn", direction = 1, limits = c(.5, 1))+
  # scale_fill_viridis_c(option = "turbo")+
  labs(fill = "Probability of Effect", x = "Nitrogen Deposition (kg N/km^2)", y = "Urban. Land Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text.y.right = element_blank())
# urb_nit_prob



ag_urb_trend_tag <- tag_facet(ag_urb_trend)
ag_nit_trend_tag <- tag_facet(ag_nit_trend, tag_pool =  letters[-(1:3)])
urb_nit_trend_tag <- tag_facet(urb_nit_trend, tag_pool =  letters[-(1:6)])

ag_urb_prob_tag <- tag_facet(ag_urb_prob, tag_pool =  letters[-(1:9)])
ag_nit_prob_tag <- tag_facet(ag_nit_prob, tag_pool =  letters[-(1:12)])
urb_nit_prob_tag <- tag_facet(urb_nit_prob, tag_pool =  letters[-(1:15)])



Fig_trends_2017_suppA <- ag_urb_trend_tag + ag_nit_trend_tag + urb_nit_trend_tag + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Rate of Change in Endophyte Prevalence")& theme(plot.title = element_text(size = rel(1.5)))
Fig_trends_2017_suppB<- ag_urb_prob_tag + ag_nit_prob_tag + urb_nit_prob_tag + plot_layout(nrow = 1, guides = "collect") + plot_annotation(title = "Posterior Probability of Effect")& theme(plot.title = element_text(size = rel(1.5)))

Fig_trends_2017_supp <- wrap_elements(Fig_trends_2017_suppA)/ wrap_elements(Fig_trends_2017_suppB)

ggsave(Fig_trends_2017_supp, filename = "Plots/Fig_trends_2017_supp.png", width = 13, height = 13)






























simple_trend_plot <- ggplot(avg_change)+
  # geom_hline(yintercept = 0)+
  # geom_jitter(data = avg_posteriors, aes(y = diff, x = factor(treatment_x, levels = ), fill = diff), width = .25, height = 0, color = "black",shape = 21,alpha = .7)+
  geom_linerange(aes(ymin = lwr, ymax = upr, x = (treatment_x), ), color = "black", lwd = 1)+
  geom_point(aes(y = diff_mean, x = (treatment_x), fill = diff_mean), size = 3, color = "black",shape = 21) + 
  # scale_color_distiller(palette = "RdYlBu", direction = -1)+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-99,99))+
  facet_grid(species ~ treatment_group, scales = "free")+
  guides(fill = "none")+
  labs( x= "", y= "Change in % Prevalence / Century")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))
simple_trend_plot
tagged_simple <- tag_facet2(simple_trend_plot)















avg_change <- tibble(preddata, as_tibble(year.pred)) %>% 
  mutate(treatment = case_when(PercentAg <0 ~ "Low Agr.", PercentAg>0 ~ "High Agr.",
                               PercentUrban <0 ~ "Low Urb.", PercentUrban>0 ~ "High Urb.",
                               mean_TIN_10km<0 ~ "Low Nit.", mean_TIN_10km>0 ~ "High Nit.", TRUE ~ "Avg. Covariates"),
         treatment_group = factor(case_when(grepl("Agr.", treatment) ~ "Agr. Land Cover (%)",
                                            grepl("Urb.", treatment) ~ "Urb. Land Cover (%)",
                                            grepl("Nit.", treatment) ~ "Nit. Deposition (kg N/km^2/year)"), levels = c("Agr. Land Cover (%)", "Urb. Land Cover (%)", "Nit. Deposition (kg N/km^2/year)")),
         treatment_label = case_when(grepl("Low", treatment) ~ "Low",
                                     grepl("High", treatment) ~ "High"),
         treatment_label = factor(treatment_label, level = c("Low", "High")),
         treatment_x = case_when(treatment_group == "Agr. Land Cover (%)" ~ paste0(treatment_label, "\n(",round((PercentAg + data_summary$PercentAg) , digits = 1), "%)"),
                                 treatment_group == "Urb. Land Cover (%)" ~ paste0(treatment_label, "\n(",round((PercentUrban + data_summary$PercentUrban) , digits = 1), "%)"),
                                 treatment_group == "Nit. Deposition (kg N/km^2/year)" ~ paste0(treatment_label, "\n(",round((mean_TIN_10km + data_summary$mean_TIN_10km) , digits = 1), ")")),
         treatment_x = factor(treatment_x, level = c("Low\n(0%)" ,  "Low\n(111.7)" ,  "High\n(95.6%)" ,  "High\n(99.1%)",  "High\n(639.1)"))) %>% 
  #        treatment_x = factor(treatment_x, level = c(paste0("Low", " (",round((unique(PercentAg) + data_summary$PercentAg) , digits = 1), "%)"), paste0("High", " (",round((unique(PercentAg) + data_summary$PercentAg) , digits = 1), "%)"),
  #                                                    paste0("Low", " (",round((unique(PercentUrban) + data_summary$PercentUrban) , digits = 1), "%)"), paste0("High", " (",round((unique(PercentUrban) + data_summary$PercentUrban) , digits = 1), "%)"),
  #                                                    paste0("Low", " (",round((unique(mean_TIN_10km) + data_summary$mean_TIN_10km) , digits = 1), " kg N/km^2/year)"), paste0("High", " (",unique((mean(mean_TIN_10km) + data_summary$mean_TIN_10km) , digits = 1), " kg N/km^2/year)")))) %>% 
  # # mutate(year = year + data_summary$year,
  #        PercentUrban = PercentUrban + data_summary$PercentUrban,
  #        PercentAg= PercentAg + data_summary$PercentAg,
  #        mean_TIN_10km= mean_TIN_10km + data_summary$mean_TIN_10km) %>% 
  pivot_longer(cols = iter1:iter500, names_to = "iteration", values_to = "posterior") %>% 
  pivot_wider(id_cols = c(Spp_code, PercentAg, PercentUrban, mean_TIN_10km, ppt_10km, tmean_10km, treatment, treatment_group,treatment_label, treatment_x, collector_index, scorer_index, species, iteration), names_from = c(year_label), values_from = posterior, names_prefix = "post.") %>%  
  mutate(diff = (post.max_year - post.min_year)*100) %>% 
  group_by(Spp_code, treatment, treatment_group,treatment_label, treatment_x, species) %>% 
  dplyr::summarise(diff_mean = mean(diff),
                   diff_median = median(diff),
                   lwr = quantile(diff, .025),
                   upr = quantile(diff, .975),
                   prob_pos = sum(diff>0)/500) %>% 
  filter(treatment != "Avg. Covariates")

avg_posteriors <- tibble(preddata, as_tibble(year.pred)) %>% 
  mutate(treatment = case_when(PercentAg <0 ~ "Low Agr.", PercentAg>0 ~ "High Agr.",
                               PercentUrban <0 ~ "Low Urb.", PercentUrban>0 ~ "High Urb.",
                               mean_TIN_10km<0 ~ "Low Nit.", mean_TIN_10km>0 ~ "High Nit.", TRUE ~ "Avg. Covariates"),
         treatment_group = factor(case_when(grepl("Agr.", treatment) ~ "Agr. Land Cover (%)",
                                            grepl("Urb.", treatment) ~ "Urb. Land Cover (%)",
                                            grepl("Nit.", treatment) ~ "Nit. Deposition (kg N/km^2/year)"), levels = c("Agr. Land Cover (%)", "Urb. Land Cover (%)", "Nit. Deposition (kg N/km^2/year)")),
         treatment_label = case_when(grepl("Low", treatment) ~ "Low",
                                     grepl("High", treatment) ~ "High"),
         treatment_label = factor(treatment_label, level = c("Low", "High")),
         treatment_x = case_when(treatment_group == "Agr. Land Cover (%)" ~ paste0(treatment_label, "\n(",round((PercentAg + data_summary$PercentAg) , digits = 1), "%)"),
                                 treatment_group == "Urb. Land Cover (%)" ~ paste0(treatment_label, "\n(",round((PercentUrban + data_summary$PercentUrban) , digits = 1), "%)"),
                                 treatment_group == "Nit. Deposition (kg N/km^2/year)" ~ paste0(treatment_label, "\n(",round((mean_TIN_10km + data_summary$mean_TIN_10km) , digits = 1), ")")),
         treatment_x = factor(treatment_x, level = c("Low\n(0%)" ,  "Low\n(111.7)" ,  "High\n(95.6%)" ,  "High\n(99.1%)",  "High\n(639.1)"))) %>% 
  # mutate(year = year + data_summary$year,
  #        PercentUrban = PercentUrban + data_summary$PercentUrban,
  #        PercentAg= PercentAg + data_summary$PercentAg,
  #        mean_TIN_10km= mean_TIN_10km + data_summary$mean_TIN_10km) %>% 
  pivot_longer(cols = iter1:iter500, names_to = "iteration", values_to = "posterior") %>% 
  pivot_wider(id_cols = c(Spp_code, PercentAg, PercentUrban, mean_TIN_10km, ppt_10km, tmean_10km, treatment, treatment_group, treatment_label, treatment_x, collector_index, scorer_index, species, iteration), names_from = c(year_label), values_from = posterior, names_prefix = "post.") %>% 
  mutate(diff = (post.max_year - post.min_year)*100) %>% 
  group_by(Spp_code, treatment, treatment_group, treatment_label, treatment_x, species) %>% 
  filter(treatment != "Avg. Covariates") %>% 
  sample_n(size = 100)

tag_facet2 <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                       hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
}
# colors <- c("#B38600", "#021475", "#BF00A0")
simple_trend_plot <- ggplot(avg_change)+
  geom_hline(yintercept = 0)+
  geom_jitter(data = avg_posteriors, aes(y = diff, x = factor(treatment_x, levels = ), fill = diff), width = .25, height = 0, color = "black",shape = 21,alpha = .7)+
  geom_linerange(aes(ymin = lwr, ymax = upr, x = (treatment_x), ), color = "black", lwd = 1)+
  geom_point(aes(y = diff_mean, x = (treatment_x), fill = diff_mean), size = 3, color = "black",shape = 21) + 
  # scale_color_distiller(palette = "RdYlBu", direction = -1)+
  scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-99,99))+
  facet_grid(species ~ treatment_group, scales = "free")+
  guides(fill = "none")+
  labs( x= "", y= "Change in % Prevalence / Century")+
  theme_bw()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))
simple_trend_plot
tagged_simple <- tag_facet2(simple_trend_plot)
ggsave(tagged_simple, filename = "Plots/temporal_trend_plot.png", width = 8.5, height =7)  


