# Purpose: Fits top models and vizualizes results
# Authors: Joshua Fowler and Mallory Tucker
# Updated: Dec 23, 2025



library(devtools)
# library("devtools")
# devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE, force = TRUE)
#INLA relies on Rgraphviz (and other packages, you can use Bioconductor to help install)
library(dplyr)
library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
# devtools::install_github('timcdlucas/INLAutils')
library(INLAutils) # supposedly has a function to plot residuals, might not need?
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
library(egg) # for labelling panels
library(ggmap)
library(pROC)
library(ggplot2)
library(maps)

invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")

species_codes <- c("AGHY", "AGPE", "ELVI")
species_names <- c("A. hyemalis", "A. perennans", "E. virginicus")




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
  select(Sample_id,spec_PercentUrban, spec_PercentAg)



endo_herb_georef1 <- left_join(endo_herb_georef, yearly_endo, by = c("Sample_id"), relationship = "many-to-many")


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
  select(lon, lat, buffer, mean_TIN, mean_NO3, mean_NH4) %>% 
  pivot_wider(id_cols = c(lon, lat), names_from = c("buffer"), values_from = c("mean_TIN", "mean_NO3", "mean_NH4"))
endo_herb_merge2 <- left_join(endo_herb_merge1, nit_avgs, by = c("lon", "lat"))


nit_yearly <- read.csv(file = paste0(path,"nitrogen_yearly_df.csv")) %>% 
  select(lon, lat, year, buffer, TIN, NO3, NH4) %>% 
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

endo_herb <- left_join(endo_herb_sf, climate, by = c("lon", "lat", "year"))




#load in map outline data
outline_map <- map_data("world")
states_shape <- map_data("state")
counties <- map_data("county")



##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
data_summary <- endo_herb %>% 
  dplyr::summarize(year = mean(year, na.rm = T),
                   PercentUrban = mean(PercentUrban, na.rm = T),
                   PercentAg = mean(PercentAg, na.rm = T),
                   mean_TIN_10km = mean(mean_TIN_10km, na.rm = T),
                   tmean_10km = mean(tmean_10km, na.rm = T),
                   ppt_10km = mean(ppt_10km, na.rm = T))
data <- endo_herb %>% 
  mutate(year = year - data_summary$year,
         PercentUrban = PercentUrban - data_summary$PercentUrban,
         PercentAg = PercentAg - data_summary$PercentAg,
         mean_TIN_10km = mean_TIN_10km - data_summary$mean_TIN_10km,
         tmean_10km = tmean_10km - data_summary$tmean_10km,
         ppt_10km = ppt_10km - data_summary$ppt_10km) %>% 
  filter(!is.na(tmean_10km))

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

bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")

non_convex_bdry[[1]] <- bdry_polygon



max.edge = diff(range(coords[,1]))/(100)


mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge*4, max.edge*8), # km inside and outside
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
# ggsave(mesh_plot, filename = "mesh_plot.png", width = 6, height = 5)



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
data <- data %>% 
  mutate(Spp_index = as.numeric(as.factor(Spp_code)))

s_components.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code + mean_TIN_10km + PercentAg + PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.3 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km + Spp_code*PercentAg + Spp_code*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.4 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km*PercentAg*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.5 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.6 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.7 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.8 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*tmean_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.9 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*ppt_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.nit.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*ppt_10km*mean_TIN_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)
s_components.nit.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*tmean_10km*mean_TIN_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)
s_components.nit.3 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentUrban*mean_TIN_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)
s_components.nit.4 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg*mean_TIN_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)






s_components.climate.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code + mean_TIN_10km + PercentAg + PercentUrban + tmean_10km + ppt_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.climate.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km + Spp_code*PercentAg + Spp_code*PercentUrban + Spp_code*tmean_10km + Spp_code*ppt_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)




s_components.year.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code + year + mean_TIN_10km + PercentAg + PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.year.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*year*mean_TIN_10km + Spp_code*year*PercentAg + Spp_code*year*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

s_components.year.3 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*year*mean_TIN_10km*PercentAg*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.year.climate.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code + year + mean_TIN_10km + PercentAg + PercentUrban + tmean_10km + ppt_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)


s_components.year.climate.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*year*mean_TIN_10km + Spp_code*year*PercentAg + Spp_code*year*PercentUrban + Spp_code*year*tmean_10km + Spp_code*year*ppt_10km, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



# s_components.y.rfx <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km + Spp_code*PercentAg + Spp_code*PercentUrban, model = "fixed")+
#   year(year, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$year)), hyper = list(pc_prec)) +
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde)



s_formula <- Endo_status_liberal ~ .


# Now run the model



fit.1 <- bru(s_components.1,
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

fit.2 <- bru(s_components.2,
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

fit.3 <- bru(s_components.3,
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

fit.4 <- bru(s_components.4,
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

#
fit.5 <- bru(s_components.5,
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

fit.6 <- bru(s_components.6,
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

fit.7 <- bru(s_components.7,
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

fit.8 <- bru(s_components.8,
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

fit.9 <- bru(s_components.9,
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



fit.nit.1 <- bru(s_components.nit.1,
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

fit.nit.2 <- bru(s_components.nit.2,
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

fit.nit.3 <- bru(s_components.nit.3,
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

fit.nit.4 <- bru(s_components.nit.4,
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


fit.climate.1 <- bru(s_components.climate.1,
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


fit.climate.2 <- bru(s_components.climate.2,
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





fit.year.1 <- bru(s_components.year.1,
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


fit.year.2 <- bru(s_components.year.2,
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


fit.year.3 <- bru(s_components.year.3,
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




fit.year.climate.1 <- bru(s_components.year.climate.1,
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


fit.year.climate.2 <- bru(s_components.year.climate.2,
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

# components
components <- c()

components[1] <- as.character(c(fit.1$bru_info$model$effects$fixed$main$input$input))
components[2] <- as.character(c(fit.2$bru_info$model$effects$fixed$main$input$input))
components[3] <- as.character(c(fit.3$bru_info$model$effects$fixed$main$input$input))
components[4] <- as.character(c(fit.4$bru_info$model$effects$fixed$main$input$input))
components[5] <- as.character(c(fit.5$bru_info$model$effects$fixed$main$input$input))
components[6] <- as.character(c(fit.6$bru_info$model$effects$fixed$main$input$input))
components[7] <- as.character(c(fit.7$bru_info$model$effects$fixed$main$input$input))
components[8] <- as.character(c(fit.8$bru_info$model$effects$fixed$main$input$input))
components[9] <- as.character(c(fit.9$bru_info$model$effects$fixed$main$input$input))
components[10] <- as.character(c(fit.climate.1$bru_info$model$effects$fixed$main$input$input))
components[11] <- as.character(c(fit.climate.2$bru_info$model$effects$fixed$main$input$input))
components[12] <- as.character(c(fit.year.1$bru_info$model$effects$fixed$main$input$input))
components[13] <- as.character(c(fit.year.2$bru_info$model$effects$fixed$main$input$input))
components[14] <- as.character(c(fit.year.3$bru_info$model$effects$fixed$main$input$input))
components[15] <- as.character(c(fit.year.climate.1$bru_info$model$effects$fixed$main$input$input))
components[16] <- as.character(c(fit.year.climate.2$bru_info$model$effects$fixed$main$input$input))
components[17] <- as.character(c(fit.nit.1$bru_info$model$effects$fixed$main$input$input))
components[18] <- as.character(c(fit.nit.2$bru_info$model$effects$fixed$main$input$input))
components[19] <- as.character(c(fit.nit.3$bru_info$model$effects$fixed$main$input$input))
components[20] <- as.character(c(fit.nit.4$bru_info$model$effects$fixed$main$input$input))


# DIC
dic <- c()
dic[1] <- fit.1$dic$dic
dic[2] <- fit.2$dic$dic
dic[3] <- fit.3$dic$dic
dic[4] <- fit.4$dic$dic
dic[5] <- fit.5$dic$dic
dic[6] <- fit.6$dic$dic
dic[7] <- fit.7$dic$dic
dic[8] <- fit.8$dic$dic
dic[9] <- fit.9$dic$dic
dic[10] <- fit.climate.1$dic$dic
dic[11] <- fit.climate.2$dic$dic
dic[12] <- fit.year.1$dic$dic
dic[13] <- fit.year.2$dic$dic
dic[14] <- fit.year.3$dic$dic
dic[15] <- fit.year.climate.1$dic$dic
dic[16] <- fit.year.climate.2$dic$dic
dic[17] <- fit.nit.1$dic$dic
dic[18] <- fit.nit.2$dic$dic
dic[19] <- fit.nit.3$dic$dic
dic[20] <- fit.nit.4$dic$dic


# waic
waic <- c()
waic[1] <- fit.1$waic$waic
waic[2] <- fit.2$waic$waic
waic[3] <- fit.3$waic$waic
waic[4] <- fit.4$waic$waic
waic[5] <- fit.5$waic$waic
waic[6] <- fit.6$waic$waic
waic[7] <- fit.7$waic$waic
waic[8] <- fit.8$waic$waic
waic[9] <- fit.9$waic$waic
waic[10] <- fit.climate.1$waic$waic
waic[11] <- fit.climate.2$waic$waic
waic[12] <- fit.year.1$waic$waic
waic[13] <- fit.year.2$waic$waic
waic[14] <- fit.year.3$waic$waic
waic[15] <- fit.year.climate.1$waic$waic
waic[16] <- fit.year.climate.2$waic$waic
waic[17] <- fit.nit.1$waic$waic
waic[18] <- fit.nit.2$waic$waic
waic[19] <- fit.nit.3$waic$waic
waic[20] <- fit.nit.4$waic$waic


# cpo
cpo <- c()
cpo[1] <-  -mean(log(fit.1$cpo$cpo))
cpo[2] <-  -mean(log(fit.2$cpo$cpo))
cpo[3] <-  -mean(log(fit.3$cpo$cpo))
cpo[4] <-  -mean(log(fit.4$cpo$cpo))
cpo[5] <- -mean(log(fit.5$cpo$cpo))
cpo[6] <- -mean(log(fit.6$cpo$cpo))
cpo[7] <- -mean(log(fit.7$cpo$cpo))
cpo[8] <- -mean(log(fit.8$cpo$cpo))
cpo[9] <- -mean(log(fit.9$cpo$cpo))
cpo[10] <-  -mean(log(fit.climate.1$cpo$cpo))
cpo[11] <-  -mean(log(fit.climate.2$cpo$cpo))
cpo[12] <-  -mean(log(fit.year.1$cpo$cpo))
cpo[13] <-  -mean(log(fit.year.2$cpo$cpo))
cpo[14] <-  -mean(log(fit.year.3$cpo$cpo))
cpo[15] <- -mean(log(fit.year.climate.1$cpo$cpo))
cpo[16] <- -mean(log(fit.year.climate.2$cpo$cpo))
cpo[17] <- -mean(log(fit.nit.1$cpo$cpo))
cpo[18] <- -mean(log(fit.nit.2$cpo$cpo))
cpo[19] <- -mean(log(fit.nit.3$cpo$cpo))
cpo[20] <- -mean(log(fit.nit.4$cpo$cpo))

# checking convergence
convergence <- c()
convergence[1] <- fit.1$mode$mode.status
convergence[2] <- fit.2$mode$mode.status
convergence[3] <- fit.3$mode$mode.status
convergence[4] <- fit.4$mode$mode.status
convergence[5] <- fit.5$mode$mode.status
convergence[6] <- fit.6$mode$mode.status
convergence[7] <- fit.7$mode$mode.status
convergence[8] <- fit.8$mode$mode.status
convergence[9] <- fit.9$mode$mode.status
convergence[10] <- fit.climate.1$mode$mode.status
convergence[11] <- fit.climate.2$mode$mode.status
convergence[12] <- fit.year.1$mode$mode.status
convergence[13] <- fit.year.2$mode$mode.status
convergence[14] <- fit.year.3$mode$mode.status
convergence[15] <- fit.year.climate.1$mode$mode.status
convergence[16] <- fit.year.climate.2$mode$mode.status
convergence[17] <- fit.nit.1$mode$mode.status
convergence[18] <- fit.nit.2$mode$mode.status
convergence[19] <- fit.nit.3$mode$mode.status
convergence[20] <- fit.nit.4$mode$mode.status



# calculating AUC



# predicting the training data
validation.1 <- predict(fit.1,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.2 <- predict(fit.2,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.3 <- predict(fit.3,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.4 <- predict(fit.4,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.5 <- predict(fit.5,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.6 <- predict(fit.6,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.7 <- predict(fit.7,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.8 <- predict(fit.8,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.9 <- predict(fit.9,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 

validation.climate.1 <- predict(fit.1,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.climate.2 <- predict(fit.2,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.year.1 <- predict(fit.year.1,newdata = data,
                                formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.year.2 <- predict(fit.year.2,newdata = data,
                                formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.year.3 <- predict(fit.year.3,newdata = data,
                             formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 

validation.year.climate.1 <- predict(fit.year.climate.1,newdata = data,
                             formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.year.climate.2 <- predict(fit.year.climate.2,newdata = data,
                                     formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.nit.1 <- predict(fit.nit.1,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.nit.2 <- predict(fit.nit.2,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.nit.3 <- predict(fit.nit.3,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 
validation.nit.4 <- predict(fit.nit.4,newdata = data,
                        formula = ~ invlogit(fixed + scorer + collector + space_int), n.samples = 500) 



rocobj.1 <- pROC::roc(data$Endo_status_liberal, validation.1$mean)
rocobj.2 <- pROC::roc(data$Endo_status_liberal, validation.2$mean)
rocobj.3 <- pROC::roc(data$Endo_status_liberal, validation.3$mean)
rocobj.4 <- pROC::roc(data$Endo_status_liberal, validation.4$mean)
rocobj.5 <- pROC::roc(data$Endo_status_liberal, validation.5$mean)
rocobj.6 <- pROC::roc(data$Endo_status_liberal, validation.6$mean)
rocobj.7 <- pROC::roc(data$Endo_status_liberal, validation.7$mean)
rocobj.8 <- pROC::roc(data$Endo_status_liberal, validation.8$mean)
rocobj.9 <- pROC::roc(data$Endo_status_liberal, validation.9$mean)
rocobj.10 <- pROC::roc(data$Endo_status_liberal, validation.climate.1$mean)
rocobj.11 <- pROC::roc(data$Endo_status_liberal, validation.climate.1$mean)
rocobj.12 <- pROC::roc(data$Endo_status_liberal, validation.year.1$mean)
rocobj.13 <- pROC::roc(data$Endo_status_liberal, validation.year.2$mean)
rocobj.14 <- pROC::roc(data$Endo_status_liberal, validation.year.3$mean)
rocobj.15 <- pROC::roc(data$Endo_status_liberal, validation.year.climate.1$mean)
rocobj.16 <- pROC::roc(data$Endo_status_liberal, validation.year.climate.2$mean)
rocobj.17 <- pROC::roc(data$Endo_status_liberal, validation.nit.1$mean)
rocobj.18 <- pROC::roc(data$Endo_status_liberal, validation.nit.2$mean)
rocobj.19 <- pROC::roc(data$Endo_status_liberal, validation.nit.3$mean)
rocobj.20 <- pROC::roc(data$Endo_status_liberal, validation.nit.4$mean)


# AUC values
auc <- c()
auc[1] <- rocobj.1$auc
auc[2] <- rocobj.2$auc
auc[3] <- rocobj.3$auc
auc[4] <- rocobj.4$auc
auc[5] <- rocobj.5$auc
auc[6] <- rocobj.6$auc
auc[7] <- rocobj.7$auc
auc[8] <- rocobj.8$auc
auc[9] <- rocobj.9$auc
auc[10] <- rocobj.10$auc
auc[11] <- rocobj.11$auc
auc[12] <- rocobj.12$auc
auc[13] <- rocobj.13$auc
auc[14] <- rocobj.14$auc
auc[15] <- rocobj.15$auc
auc[16] <- rocobj.16$auc
auc[17] <- rocobj.17$auc
auc[18] <- rocobj.18$auc
auc[19] <- rocobj.19$auc
auc[20] <- rocobj.20$auc

# model comparison table
table <- data.frame(components, dic, waic, cpo, auc, convergence) %>% arrange(dic, waic, cpo, auc) %>% 
  mutate(across(-c(components), ~round(.,digits = 4))) 
table <- table %>% 
  mutate(components = factor(components, levels = (table$components)[rev(order(table$dic))]),
         year = grepl("year", components))
write.csv(table, "Analyses/model.comparison.csv")


fill_colors <- RColorBrewer::brewer.pal(4,"Set1")
       
dic_plot <- ggplot(table%>%  select(components, dic, year))+
  geom_tile(aes(y = components, x = 1, fill = dic), color = "white", alpha = .8)+
  geom_text(aes(y = components, x = 1, label = dic), size = 3)+
  labs(x = "DIC",y = "", fill = "DIC",  title = "A")+
  scale_fill_distiller(palette = "Reds", direction = -1)+
  guides(fill = "none")+
  theme_void()+
  theme(axis.text.y = element_text( size = rel(.6), hjust = 1 ),
        axis.text.x = element_blank(),
        axis.title.x = element_text(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
# dic_plot


waic_plot <- ggplot(table %>% select(components, waic, year))+
  geom_tile(aes(y = components, x = 1, fill = waic),color = "white", alpha = .8)+
  geom_text(aes(y = components, x = 1, label = waic), size = 3)+
  labs(x = "WAIC",y = "", fill = "WAIC",  title = "B")+
  scale_fill_distiller(palette = "Blues", direction = -1)+
  guides(fill = "none")+
  theme_void()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
# waic_plot

cpo_plot <- ggplot(table %>% select(components, cpo, year))+
  geom_tile(aes(y = components, x = 1, fill = cpo),color = "white", alpha = .8)+
  geom_text(aes(y = components, x = 1, label = cpo), size = 3)+
  labs(x = "-log(CPO)",y = "", fill = "CPO",  title = "C")+
  scale_fill_distiller(palette = "Greens", direction = -1)+
  guides(fill = "none")+
  theme_void()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
# cpo_plot


auc_plot <- ggplot(table %>% select(components, auc, year))+
  geom_tile(aes(y = components, x = 1, fill = auc),color = "white", alpha = .8)+
  geom_text(aes(y = components, x = 1, label = auc), size = 3)+
  guides(fill = "none")+
  scale_fill_distiller(palette = "Greys", direction = -1)+
  theme_void()+
  labs(x = "AUC", y = "",fill = "AUC", title = "D")+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_text(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
# auc_plot


model_comparison_plot <- dic_plot + waic_plot + cpo_plot + auc_plot +plot_layout(nrow = 1) 
ggsave(model_comparison_plot, filename = "Plots/model_comparison_plot.png", width = 12, height = 5)

