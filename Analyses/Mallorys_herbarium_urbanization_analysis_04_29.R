# Purpose: Fits spatial model of change in Epichloe endophyte prevalence in herbarium specimens and assesses effect of different anthropogenic land-use on trends.
# Authors: Mallory Tucker and Joshua Fowler
# Updated: Feb 13, 2024

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
library(patchwork)
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
Joshpath <- "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/"
path <- Joshpath


# endo_herb_georef <- read_csv(file = paste0(path, "Zonalhist_NLCD_10km_.csv")) %>%
#   filter(Country != "Canada") %>%
#   mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
#                               grepl("ELVI", Sample_id) ~ "ELVI",
#                               grepl("AGPE", Sample_id) ~ "AGPE")) %>%
#   mutate(species_index = as.factor(case_when(Spp_code == "AGHY" ~ "1",
#                                              Spp_code == "AGPE" ~ "2",
#                                              Spp_code == "ELVI" ~ "3"))) %>%
#   mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
#                              Spp_code == "AGPE" ~ "A. perennans",
#                              Spp_code == "ELVI" ~ "E. virginicus")) %>%
#   mutate(decade = floor(year/10)*10)%>%
#   mutate(DevelopedOpenSpace = HISTO_21,
#          DevelopedLowIntensity = HISTO_22,
#          MediumDeveloped = HISTO_23,
#          HighDeveloped = HISTO_24,
#          PastureHay = HISTO_81,
#          CultivatedCrops = HISTO_82,
#          TotalAg = PastureHay + CultivatedCrops,
#          TotalPixels = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24 + HISTO_0 + HISTO_11 +HISTO_12 + HISTO_31 +HISTO_41 + HISTO_42 + HISTO_43 + HISTO_52 + HISTO_71+ HISTO_90 + TotalAg + HISTO_95,
#          TotalDeveloped = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24,
#          OtherLC = (TotalPixels - (TotalAg + TotalDeveloped))/TotalPixels *100,
#          PercentUrban = TotalDeveloped/TotalPixels * 100,
#          PercentAg = TotalAg/TotalPixels * 100)
# 
# # Doing some filtering to remove NA's and some data points that probably aren't accurate species id's
# endo_herb <- endo_herb_georef %>%
#   filter(!is.na(Endo_statu)) %>%
#   filter(!is.na(Spp_code)) %>%
#   filter(!is.na(lon) & !is.na(year)) %>%
#   filter(lon>-110 ) %>%
#   filter(Country != "Canada" ) %>%
#   mutate(year_bin = case_when(year<1970 ~ "pre-1970",
#                               year>=1970 ~ "post-1970")) %>%
#   mutate(endo_status_text = case_when(Endo_statu == 0 ~ "E-",
#                                       Endo_statu == 1 ~ "E+"))
# # #loading in nitrogen data too
# nit <- read.csv(file = "endo_herb_nit.csv")
# nitendoherb <- merge(nit, endo_herb)
# endo_herb <- nitendoherb
# write.csv(endo_herb, file = "EndoHerb_withNitrogen.csv")

endo_herb <- read_csv(file = "Analyses/EndoHerb_withNitrogen.csv") %>%
  mutate(sample_temp = Sample_id) %>%
  separate(sample_temp, into = c("Herb_code", "spp_code", "specimen_code", "tissue_code")) %>%
  mutate(species_index = as.factor(case_when(spp_code == "AGHY" ~ "1",
                                             spp_code == "AGPE" ~ "2",
                                             spp_code == "ELVI" ~ "3"))) %>%
  mutate(species = case_when(spp_code == "AGHY" ~ "A. hyemalis",
                             spp_code == "AGPE" ~ "A. perennans",
                             spp_code == "ELVI" ~ "E. virginicus")) %>%
  mutate(std_year = (year-mean(year, na.rm = T)),
         std_nit = (NO3_mean - mean(NO3_mean, na.rm = T))/sd(NO3_mean, na.rm = T),
         std_urb = (PercentUrban - mean(PercentUrban, na.rm = T))/sd(PercentUrban, na.rm = T),
         std_ag = (PercentAg - mean(PercentAg, na.rm = T))/sd(PercentAg, na.rm = T)) %>%  # I am mean centering but not scaling by standard deviation to preserve units for interpretation of the parameter values
  filter(scorer_id != "Scorer26") %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(spp_code)) %>%
  filter(!is.na(lon) & !is.na(year)) %>%
  filter(!is.na(PercentAg), !is.na(NO3_mean))

# Looking at just AGHY for now, but we can likely fit one model for all species
# endo_herb <- endo_herb %>% 
#   filter(Spp_code == "AGHY")

# updating the collector labels

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

endo_herb<- endo_herb %>% 
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
  st_transform(epsg6703km) %>% 
  mutate(
    easting = st_coordinates(.)[, 1],
    northing = st_coordinates(.)[, 2]
  ) %>% 
  mutate(scorer_index = parse_number(scorer_factor),
         collector_index = parse_number(collector_factor)) 



# register_google(key = "")
# map <- ggmap::get_map(zoom = 3, source = "google", maptype = c("satellite"))

outline_map <- map_data("world")
states_shape <- map_data("state")
counties_shape <- map_data("county")
# ggplot()+
#   geom_map(data = counties_shape, map = counties_shape, aes(long, lat, map_id = region))

collections_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = species), alpha = .7, lwd = .5)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = species_colors)+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

collections_map
# ggsave(collections_map, filename = "collections_map.png", width = 7, height = 4)

endo_status_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = endo_status_text), alpha = .7, lwd = .8)+
  facet_wrap(~factor(year_bin, levels = c("pre-1970", "post-1970"))+species)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = c(endophyte_colors[2],endophyte_colors[6]))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
endo_status_map

# ggsave(endo_status_map, filename = "endo_status_map.png", width = 10, height = 5)



# plotting the relationship between percentag and percenturban

endo_herb_disturbance <- endo_herb %>% 
  mutate(nondisturbance_ratio = PercentAg/PercentUrban,
         disturbance_ratio = log(PercentAg/PercentUrban))
disturbance <- ggplot(endo_herb_disturbance)+
  geom_point(aes(x = PercentAg, y = PercentUrban, color = disturbance_ratio))+
  scale_color_gradient2(low = "#542788", mid = "grey", high = "#b35806")+
  geom_point(aes(x = quantile(PercentAg, .05), y = quantile(PercentUrban, .95)), color = "#542788", size = 4, shape = 4)+
  geom_point(aes(y = quantile(PercentUrban, .05), x = quantile(PercentAg, .95)), color = "#b35806", size = 4, shape = 4)+
  geom_point(aes(y = quantile(PercentUrban, .05), x = quantile(PercentAg, .05)), color = "black", size = 4, shape = 4)+
  labs(x = "Percent Ag. (%)", y = "Percent Urban (%)", color = "Disturbance \n Gradient")+
  theme_bw()

disturbance
ggsave(disturbance, filename = "disturbance_dimensions.png")
cor(endo_herb$PercentAg, endo_herb$PercentUrban)

cor(endo_herb$PercentAg, endo_herb$PercentUrban, method = "spearman")
ggplot(endo_herb)+
  geom_point(aes(x = PercentAg, y = NO3_mean, color = Spp_code))


ggplot(endo_herb)+
  geom_point(aes(x = PercentUrban, y = NO3_mean))


mode <- function(codes){
  which.max(tabulate(codes))
}

summary_endo_herb <- endo_herb %>% 
  mutate(Sample_id_temp = Sample_id) %>% 
  separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
  # filter(seed_scored>0) %>% 
  filter(month<=12&month>0) %>% 
  group_by(species) %>% 
  dplyr::summarize(n(),
            avg_seed = mean(seed_scored, na.rm = T),
            avg_month = mode(as.numeric(month)))

#########################################################################################
####################### Map Figure Code ################################################
########################################################################################

#load in map outline data
outline_map <- map_data("world")
states_shape <- map_data("state")
counties <- map_data("county")

#Join map outline data with endo_herb data
endo_herb$region <- tolower(paste(endo_herb$State, endo_herb$County, sep = ","))

tiny_df_urb <- data.frame(endo_herb$region, endo_herb$PercentUrban)
counties$region <- tolower(paste(counties$region, counties$subregion, sep = ","))
colnames(tiny_df_urb) <- c("region", "PercentUrban")

counties_data <- counties %>%
  left_join(tiny_df_urb, by = "region")

counties_data <- na.omit(counties_data)

# Create the choropleth map
endo_urb_map <- ggplot(data = counties_data) +
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "darkgray", linewidth = .3, fill = "#E8E7E7")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "gray", linewidth = .2, fill = "#E8E7E7")+
  geom_polygon(aes(x = long, y = lat, group = group, fill = PercentUrban), color = "darkgray", linewidth = .1) +
  coord_cartesian(xlim = c(-109, -68), ylim = c(25, 48))+
  scale_fill_gradient(low = "lightgray", high = "#021475", na.value = NA)+
  labs(x = "Longitude", y = "Latitude", fill = "Urban Cover %")+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

###############################################################################
############################ agricultural map ################################
###############################################################################

#Join map outline data with endo_herb data
tiny_df_ag <- data.frame(endo_herb$region, endo_herb$PercentAg)
colnames(tiny_df_ag) <- c("region", "PercentAg")

counties_data_ag <- counties %>%
  left_join(tiny_df_ag, by = "region")

counties_data_ag <- na.omit(counties_data_ag)

# Create the choropleth map
endo_ag_map <- ggplot(data = counties_data_ag) +
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "darkgray", linewidth = .3, fill = "#E8E7E7")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "gray", linewidth = .2, fill = "#E8E7E7")+
  geom_polygon(aes(x = long, y = lat, group = group, fill = PercentAg), color = "darkgray", linewidth = .1) +
  # coord_fixed(ratio = 1) +
  coord_cartesian(xlim = c(-109, -68), ylim = c(25, 48))+
  scale_fill_gradient(low = "lightgray", high = "#B38600", na.value = NA)+
  theme_light() +
  theme(legend.text = element_text(face = "italic")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x = "Longitude", y = "Latitude", fill = "Agricultural %")


#################################################################################
################################## nitrogen map#################################
################################################################################

#Join map outline data with endo_herb data
tiny_df_nit <- data.frame(endo_herb$region, endo_herb$NO3_mean)
colnames(tiny_df_nit) <- c("region", "NO3_mean")

counties_data_nit <- counties %>%
  left_join(tiny_df_nit, by = "region")

counties_data_nit <- na.omit(counties_data_nit)

# Create the choropleth map
endo_nit_map <- ggplot(data = counties_data_nit) +
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "darkgray", linewidth = .3, fill = "#E8E7E7")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "gray", linewidth = .2, fill = "#E8E7E7")+
  geom_polygon(aes(x = long, y = lat, group = group, fill = NO3_mean), color = "darkgray", linewidth = .1) +
  # coord_fixed(ratio = 1) +
  coord_cartesian(xlim = c(-109, -68), ylim = c(25, 48))+
  scale_fill_gradient(low = "lightgray", high = "#BF00A0", na.value = NA)+
  theme_light() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.text = element_text(face = "italic")) +
  labs(x = "Longitude", y = "Latitude", fill = "Nitrogen Dep")

#Compile maps into one panel
mapfig <-  endo_ag_map +endo_urb_map + endo_nit_map + plot_layout(ncol = 3) + plot_annotation(tag_levels = "A")
# mapfig <-  endo_ag_map /endo_urb_map / endo_nit_map + plot_layout(ncol = 1) + plot_annotation(tag_levels = "A")

#Save map file
ggsave(mapfig, file = "Map_Figure.png", width = 5, height = 10)
ggsave(mapfig, file = "Map_Figure.png", width = 18, height = 5)


##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)

data <- endo_herb

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
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*4), # km inside and outside
  cutoff = max.edge,
  crs = fm_crs(data)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge
# plot it
#plot(mesh)


mesh_plot <- ggplot() +
  gg(data = mesh) +
  geom_point(data = data, aes(x = easting, y = northing, col = species), size = 1) +
  coord_sf()+
  theme_bw() +
  labs(x = "", y = "")
# mesh_plot
# ggsave(mesh_plot, filename = "Plots/mesh_plot.png", width = 6, height = 6)



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

# s_components.1 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg*std_year*PercentUrban, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 
#   
# s_components.2 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg*std_year + Spp_code*PercentUrban*std_year, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 
# 
s_components.3 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*std_year*NO3_mean, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)

# s_components.4 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg + Spp_code*PercentUrban + Spp_code*NO3_mean, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 


s_components.4 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg + Spp_code*PercentUrban + Spp_code*NO3_mean, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde) 


# s_components.4 <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*std_ag + Spp_code*std_urb + Spp_code*std_nit, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 

# s_components.5 <-  ~ 0 + species(main = ~ 0 + Spp_code, model = "fixed") +  
#   species.ag(main = ~ 0 + Spp_code:PercentAg, model = "fixed") + species.urb(main = ~ 0 + Spp_code:PercentUrban, model = "fixed") + species.nit(main = ~ 0 +Spp_code:NO3_mean, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 
# 
# s_components.5 <-  ~ 0 + fixed(main = ~ 0 + Spp_code*PercentAg*std_year + Spp_code*PercentUrban*std_year + Spp_code*NO3_mean*std_year, model = "fixed")+
#   scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
#   collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
#   space_int(coords, model = spde) 
# putting the components with the formula
s_formula <- Endo_status_liberal ~ .


# Now run the model

# fit.1 <- bru(s_components.1,
#            like(
#              formula = s_formula,
#              family = "binomial",
#              Ntrials = 1,
#              data = data
#            ),
#            options = list(
#              control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
#              control.inla = list(int.strategy = "eb"),
#              verbose = TRUE
#            )
# )
# 
# fit.2 <- bru(s_components.2,
#              like(
#                formula = s_formula,
#                family = "binomial",
#                Ntrials = 1,
#                data = data
#              ),
#              options = list(
#                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
#                control.inla = list(int.strategy = "eb"),
#                verbose = TRUE
#              )
# )


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

# fit.5 <- bru(s_components.5,
#              like(
#                formula = s_formula,
#                family = "binomial",
#                Ntrials = 1,
#                data = data
#              ),
#              options = list(
#                control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
#                control.inla = list(int.strategy = "eb"),
#                verbose = TRUE
#              )
# )

fit.1$dic$dic
fit.2$dic$dic
fit.3$dic$dic
fit.4$dic$dic
fit.5$dic$dic


fit.1$mode$mode.status # a 0 or low value indicates "convergence"
fit.2$mode$mode.status # a 0 or low value indicates "convergence"

fit.1$summary.fixed
fit.1$summary.random

saveRDS(fit.2, file = "fit_wo_N.rds")
fit.2 <- read_rds(file = "fit_wo_N.rds")






################################################################################################################################
##########  Plotting the prediction without year effects ###############
################################################################################################################################

min_ag<- min(data$PercentAg)
mean_ag <- mean(data$PercentAg)
max_ag <- max(data$PercentAg)

min_urb<- min(data$PercentUrban)
mean_urb <- mean(data$PercentUrban)
max_urb<- max(data$PercentUrban)

min_nit<- min(data$NO3_mean)
mean_nit <- mean(data$NO3_mean)
max_nit<- max(data$NO3_mean)

preddata.1 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = rep(seq(min_ag, max_ag, length.out = 50), times = 3),
                     PercentUrban = mean_urb,
                     NO3_mean = mean_nit,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))
preddata.2 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = mean_ag,
                     PercentUrban = rep(seq(min_urb, max_urb, length.out = 50), times = 3),
                     NO3_mean = mean_nit,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

preddata.3 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = mean_ag,
                     PercentUrban = mean_urb,
                     NO3_mean = rep(seq(min_nit, max_nit, length.out = 50), times = 3),
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

ag.pred <- predict(
  fit.4,
  newdata = preddata.1,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) 

urb.pred <- predict(
  fit.4,
  newdata = preddata.2,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) 

nit.pred <- predict(
  fit.4,
  newdata = preddata.3,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) 

values <-  c("#b2abd2", "#5e3c99")

ag_trend <- ggplot(ag.pred) +
  geom_line(aes(PercentAg, mean)) +
  geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#B38600") +
  geom_ribbon(aes(PercentAg, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  # geom_point(data = ag_data_binned, aes(x = mean_ag, y = mean_endo, size = sample, fill = year_bin), color = "black", shape = 21)+
  facet_wrap(~species,  ncol = 1, scales = "free_x", strip.position="right")+  
  # scale_color_manual(values = c("#b2abd2", "#5e3c99"))+
  # scale_fill_manual(values = c("#b2abd2", "#5e3c99"))+
  labs(y = "Endophyte Prevalence", x = "Percent Ag. (%)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1), x = c(0, 100))

urb_trend <- ggplot(urb.pred) +
  geom_line(aes(PercentUrban, mean)) +
  geom_ribbon(aes(PercentUrban, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#021475") +
  geom_ribbon(aes(PercentUrban, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  # geom_point(data = ag_data_binned, aes(x = mean_ag, y = mean_endo, size = sample, fill = year_bin), color = "black", shape = 21)+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  # scale_color_manual(values = c("#b2abd2", "#5e3c99"))+
  # scale_fill_manual(values = c("#b2abd2", "#5e3c99"))+
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))


nit_trend <- ggplot(nit.pred) +
  geom_line(aes(NO3_mean, mean)) +
  geom_ribbon(aes(NO3_mean, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#BF00A0") +
  geom_ribbon(aes(NO3_mean, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  # geom_point(data = ag_data_binned, aes(x = mean_ag, y = mean_endo, size = sample, fill = year_bin), color = "black", shape = 21)+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  # scale_color_manual(values = c("#b2abd2", "#5e3c99"))+
  # scale_fill_manual(values = c("#b2abd2", "#5e3c99"))+
  labs(y = "Endophyte Prevalence", x = "Nitrogen Deposition (kg N/km^2)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(), strip.text = element_text(face = "italic", size = rel(1.1)), strip.text.y.right = element_text(angle = 0),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))


fig1 <- ag_trend + urb_trend + nit_trend + plot_layout(ncol = 3) #+ plot_annotation(tag_levels = "A")
ggsave(fig1, file = "Figure_1.png", width = 8.5, height = 8)

################################################################################################################################
##########  Plotting the posteriors from the model without year effect ###############
################################################################################################################################

param_names <- fit.4$summary.random$fixed$ID

n_draws <- 500

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors <- generate(
  fit.4,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))


posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")



# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  mutate(NIT.AGHY = NO3_mean,
         NIT.AGPE = NO3_mean+`Spp_codeAGPE:NO3_mean`,
         NIT.ELVI = NO3_mean+`Spp_codeELVI:NO3_mean`,
         AG.AGHY = PercentAg,
         AG.AGPE = PercentAg+`Spp_codeAGPE:PercentAg`,
         AG.ELVI = PercentAg+`Spp_codeELVI:PercentAg`,
         URB.AGHY = PercentUrban,
         URB.AGPE = PercentUrban+`Spp_codeAGPE:PercentUrban`,
         URB.ELVI = PercentUrban+`Spp_codeELVI:PercentUrban`,
         INT.AGHY = Spp_codeAGHY,
         INT.AGPE = Spp_codeAGPE,
         INT.ELVI = Spp_codeELVI) %>% 
  select(-all_of(param_names)) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year") %>% 
  mutate(spp_label = sub(".*\\.", "", param),
         param_label = sub("\\..*","", param))



posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_label, fill  = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_label, scales = "free")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_without_year_intxn.png", width = 10, height = 10)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/500)

################################################################################################################################
##########  Assessing model fit     ###############
################################################################################################################################

# predicting the training data

validation.pred <- predict(
  fit.4,
  newdata = endo_herb,
  formula = ~ invlogit(fixed + scorer + collector + space_int),
  n.samples = 100) 


rocobj <- pROC::roc(endo_herb$Endo_status_liberal, validation.pred$mean)

ROC_training_plot <- ggroc(rocobj) 
ggsave(ROC_training_plot, filename = "ROC_training_plot.png", width = 4, height = 4)

# AUC values
rocobj$auc
# 0.7892



################################################################################################################################
##########  Getting and plotting prediction from NitXYear ###############
################################################################################################################################
mean_year <- mean(data$year)
summary_data <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(std_year),
            max_year = max(std_year),
            min_nit = quantile(NO3_mean, .05),
            max_nit = quantile(NO3_mean, .95),
            mean_nit = mean(NO3_mean))

preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             std_year = seq(summary_data[summary_data$Spp_code == "AGHY",]$min_year, summary_data[summary_data$Spp_code == "AGHY",]$max_year, length.out = 20), 
                             NO3_mean = c(summary_data[summary_data$Spp_code == "AGHY",]$min_nit, summary_data[summary_data$Spp_code == "AGHY",]$max_nit))
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             std_year = seq(summary_data[summary_data$Spp_code == "AGPE",]$min_year, summary_data[summary_data$Spp_code == "AGPE",]$max_year, length.out = 20), 
                             NO3_mean = c(summary_data[summary_data$Spp_code == "AGPE",]$min_nit, summary_data[summary_data$Spp_code == "AGPE",]$max_nit))
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             std_year = seq(summary_data[summary_data$Spp_code == "ELVI",]$min_year, summary_data[summary_data$Spp_code == "ELVI",]$max_year, length.out = 20), 
                             NO3_mean = c(summary_data[summary_data$Spp_code == "ELVI",]$min_nit, summary_data[summary_data$Spp_code == "ELVI",]$max_nit))

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         nit_label = case_when(NO3_mean > mean(NO3_mean) ~ "High Nitrogen",
                               NO3_mean < mean(NO3_mean) ~ "Low Nitrogen"))


# gennerating predictions and back-transforming the standardized year variable



year.pred <- predict(
  fit.3,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_nit = cut(NO3_mean, breaks = 2),
         binned_year = cut(year, breaks = 40)) %>%
  group_by(Spp_code, species,binned_nit, binned_year) %>%
  summarise(mean_nit = mean(NO3_mean),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(nit_bin = case_when(mean_nit>=17.6~ "High Nitrogen",
                            mean_nit<17.6 ~ "Low Nitrogen"))




# histogram_data <- endo_herb%>%
#   mutate(disturbance_label = case_when(NO3_mean > mean(NO3_mean) ~ "Nitrogen",
#                                     NO3_mean < mean(NO3_mean) ~ "Undisturbed"))
# histogram_data2 <- endo_herb%>%
#   mutate(nit_label = case_when(NO3_mean > mean(NO3_mean) ~ "High Nitrogen",
#                                        NO3_mean < mean(NO3_mean) ~ "Low Nitrogen"))

nit_yr_trend <- ggplot(year.pred)+
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = nit_label, fill = nit_label), alpha = 0.3) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = nit_label, fill = nit_label), alpha = 0.3) +
  geom_line(aes(year, mean, group = nit_label), color = "black") +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, fill = nit_bin, size = sample), shape = 21, alpha = .6)+
  facet_wrap(~species, nrow = 3, scales = "free")+
  guides(color = "none")+
  scale_color_manual(values = c("High Nitrogen" = "#BF00A0",
                               "Low Nitrogen"="darkgray"))+
  scale_fill_manual(values = c("High Nitrogen" = "#BF00A0",
                                "Low Nitrogen"="darkgray"))+
  labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size", fill = "Nitrogen Deposition")+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        strip.text = element_text(face = "italic", size = rel(1.1)), strip.text.y.right = element_text(angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1), x = c(1832, 2020))

nit_yr_trend
ggsave(nit_yr_trend, filename = "Nitrogen_and_Year.png", width = 6, height = 8)



# year_trend_facet <- ggplot(year.pred) +
#   geom_dots(data = histogram_data, aes(x = year, y = Endo_status_liberal, color = disturbance_label, side = ifelse(Endo_status_liberal == 1, "bottom", "top")),binwidth = 1.5)+
#   # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = ag_bin))+
#   geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
#   geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
#   geom_line(aes(year, mean, group = disturbance_label, color = disturbance_label)) +
#   facet_wrap(~species+disturbance_label, nrow = 3)+
#   # scale_color_manual(values = c("#b35806", "black", "#542788"))+
#   # scale_fill_manual(values = c("#b35806", "black", "#542788"))+
#   guides(color = "none", fill = "none")+
#   labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size")+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         legend.text = element_text(face = "italic"),
#         plot.margin = unit(c(0,.1,.1,.1), "line"))+
#   lims(y = c(0,1))
# 
# year_trend_facet
# 
# ggsave(year_trend_facet, filename = "year_trend_facet_N.png", width = 10, height = 10)

################################################################################################################################
##########  Plotting the posteriors from NitXYear model ###############
################################################################################################################################


param_names <- fit.3$summary.random$fixed$ID

n_draws <- 500

posteriors <- generate(
  fit.3,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))

posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")

# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  mutate(NIT.AGHY = NO3_mean,
         NIT.AGPE = NO3_mean+`Spp_codeAGPE:NO3_mean`,
         NIT.ELVI = NO3_mean+`Spp_codeELVI:NO3_mean`,
         YEAR.AGHY = std_year,
         YEAR.AGPE = std_year+`Spp_codeAGPE:std_year`,
         YEAR.ELVI = std_year+`Spp_codeELVI:std_year`,
         YEARxNIT.AGHY = `NO3_mean:std_year`,
         YEARxNIT.AGPE = `NO3_mean:std_year`+`Spp_codeAGPE:NO3_mean:std_year`,
         YEARxNIT.ELVI = `NO3_mean:std_year`+`Spp_codeELVI:NO3_mean:std_year`,
         INT.AGHY = Spp_codeAGHY,
         INT.AGPE = Spp_codeAGPE,
         INT.ELVI = Spp_codeELVI) %>% 
  select(-all_of(param_names)) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "NitXYear") %>% 
  mutate(spp_label = sub(".*\\.", "", param),
         param_label = sub("\\..*","", param))


posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_label, scales = "free")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_with_year.png")



effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/500,
            prob_neg = sum(value<0)/500)


# trying to get an estimate of what is the difference in slope between the two plotted prediction lines for high and low nitrogen


posteriors_trend <- generate(
  fit.3,
  newdata = preddata,
  formula = ~ fixed,
  n.samples = n_draws) 
colnames(posteriors_trend) <- c( paste0("iter",1:n_draws))

posteriors_trend_df <- preddata %>% bind_cols(as_tibble(posteriors_trend)) %>% 
  pivot_longer(cols = contains("iter"), names_to = "iteration") %>% 
  mutate(year = std_year + mean_year)

  


ggplot(data = posteriors_trend_df)+
  geom_line(aes(x = year, y = value, color = nit_label, group = interaction(nit_label, iteration)))+
  facet_wrap(~Spp_code)

posteriors_trend_summary <- posteriors_trend_df %>% 
  group_by(Spp_code, nit_label, iteration) %>% 
  filter(year == min(year) | year == max(year)) %>% 
  mutate(year_label = case_when(year == min(year)~"start",
                                year == max(year)~"end")) %>% select(-year, -std_year) %>% 
  pivot_wider(names_from = year_label, values_from = value) %>% 
  mutate(slope = end-start) 

ggplot(data = posteriors_trend_summary)+
  geom_histogram(aes(x = slope, fill = nit_label))+
  facet_wrap(~Spp_code)

prob_superior <- posteriors_trend_summary %>% 
  select(-start, -end, -NO3_mean) %>% 
  pivot_wider(names_from = nit_label, names_glue = "slope_{nit_label}", values_from = slope) %>% 
  mutate(slope_diff = `slope_High Nitrogen` - `slope_Low Nitrogen`) %>% 
  group_by(Spp_code) %>% 
  summarize(mean = mean(slope_diff), 
            lwr = quantile(slope_diff, .025),
            upr = quantile(slope_diff, .975),
            prob_pos = sum(slope_diff>0)/500,
            prob_neg = sum(slope_diff<0)/500)















################################################################################################################################
########### old stuff #########
################################################################################################################################

################################################################################################################################
##########  Getting and plotting prediction across percent agriculte from the model  ###############
################################################################################################################################

min_year <- min(data$std_year)
max_year <- max(data$std_year)

mean_year <- mean(data$year)


min_ag<- min(data$PercentAg)
max_ag <- max(data$PercentAg)


urb <- quantile(data$PercentUrban, 0.05)

preddata <- expand.grid(Spp_code = c("AGHY","AGPE","ELVI"), std_year = c( 1900, 1980)-mean_year, PercentAg = seq(min_ag, max_ag, length.out = 100), PercentUrban = urb,
                        
                        collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

# gennerating predictions and back-transforming the standardized year variable



ag.pred <- predict(
  fit.2,
  newdata = preddata,
  formula = ~ invlogit(fixed), #+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)


# binning the data for plotting
ag_data_binned <- endo_herb %>% 
  mutate(binned_ag = cut(PercentAg, breaks = 25),
         binned_year = cut(year, breaks = 2)) %>%
  group_by(Spp_code, species,binned_ag, binned_year) %>%   
  summarise(mean_ag = mean(PercentAg),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_liberal)/sqrt(sample)) %>% 
  mutate(year_bin = case_when(mean_year>=1940~ paste("1980"),
                             mean_year<1940 ~ paste("1900")))


ag_trend <- ggplot(ag.pred) +
  geom_line(aes(PercentAg, mean, group = std_year, color = as.factor(year))) +
  geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975, group = std_year, fill = as.factor(year)), alpha = 0.2) +
  geom_ribbon(aes(PercentAg, ymin = q0.25, ymax = q0.75, group = std_year, fill = as.factor(year)), alpha = 0.2) +
  geom_point(data = ag_data_binned, aes(x = mean_ag, y = mean_endo, size = sample, fill = year_bin), color = "black", shape = 21)+
  facet_wrap(~species)+  
  scale_color_manual(values = c("#b2abd2", "#5e3c99"))+
  scale_fill_manual(values = c("#b2abd2", "#5e3c99"))+
  labs(y = "Endophyte Prevalence", x = "Percent Ag (%)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))

ag_trend

ggsave(ag_trend, filename = "ag_trend.png")



################################################################################################################################
##########  Getting and plotting prediction across percent urban from the model  ###############
################################################################################################################################

min_year <- min(data$std_year)
max_year <- max(data$std_year)

mean_year <- mean(data$year)


min_urb<- min(data$PercentUrban)
max_urb <- max(data$PercentUrban)


ag<- quantile(data$PercentAg, 0.05)

preddata <- expand.grid(Spp_code = c("AGHY","AGPE","ELVI"), std_year = c( 1900, 1980)-mean_year, PercentUrban = seq(min_urb, max_urb, length.out = 100), PercentAg = ag,
                        
                        collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

# gennerating predictions and back-transforming the standardized year variable



urb.pred <- predict(
  fit.2,
  newdata = preddata,
  formula = ~ invlogit(fixed), #+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)


# binning the data for plotting
urb_data_binned <- endo_herb %>% 
  mutate(binned_urb = cut(PercentUrban, breaks = 15),
         binned_year = cut(year, breaks = 2)) %>%
  group_by(Spp_code, species,binned_urb, binned_year) %>%   
  summarise(mean_urb = mean(PercentUrban),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_liberal)/sqrt(sample)) %>% 
  mutate(year_bin = case_when(mean_year>=1940~ paste("1980"),
                              mean_year<1940 ~ paste("1900")))


urb_trend <- ggplot(urb.pred) +
  geom_line(aes(PercentUrban, mean, group = std_year, color = as.factor(year))) +
  geom_ribbon(aes(PercentUrban, ymin = q0.025, ymax = q0.975, group = std_year, fill = as.factor(year)), alpha = 0.2) +
  geom_ribbon(aes(PercentUrban, ymin = q0.25, ymax = q0.75, group = std_year, fill = as.factor(year)), alpha = 0.2) +
  geom_point(data = urb_data_binned, aes(x = mean_urb, y = mean_endo, size = sample, fill = year_bin), color = "black", shape = 21)+
  facet_wrap(~species)+
  scale_shape_manual(values = c( 1, 16))+
  scale_color_manual(values = c("#b2abd2", "#5e3c99"))+
  scale_fill_manual(values = c("#b2abd2", "#5e3c99"))+
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)",shape = "Year", color = "Year", fill = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))

urb_trend

ggsave(urb_trend, filename = "urb_trend.png")







################################################################################################################################
##########  Getting and plotting prediction across year from the model  ###############
################################################################################################################################

min_year <- min(data$std_year)
max_year <- max(data$std_year)

mean_year <- mean(data$year)


min_ag<- quantile(data$PercentAg, .05)
max_ag <- quantile(data$PercentAg, .95)

min_urb <- quantile(data$PercentUrban, .05)
max_urb <- quantile(data$PercentUrban, .95)

preddata <- expand.grid(Spp_code = c("AGHY","AGPE","ELVI"), std_year = seq( min_year, max_year, length.out = 20), PercentAg = c(min_ag, max_ag), PercentUrban = c(min_urb, max_urb),
                        
                        collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3])) %>% 
  filter(PercentAg == min_ag & PercentUrban == max_urb | PercentAg == max_ag & PercentUrban == min_urb | PercentAg == min_ag & PercentUrban == min_urb) %>% 
  mutate(disturbance_label = case_when(PercentAg == min_ag & PercentUrban == min_urb ~ "Undisturbed",
                                       PercentAg == max_ag ~ "Agricultural",
                                       PercentUrban == max_urb ~ "Urban"))

# gennerating predictions and back-transforming the standardized year variable



year.pred <- predict(
  fit.2,
  newdata = preddata,
  formula = ~ invlogit(fixed),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = std_year + mean_year)


# binning the data for plotting
endo_herb_binned <- endo_herb %>% 
  mutate(binned_ag = cut(PercentAg, breaks = 2),
         binned_year = cut(year, breaks = 7)) %>%
  group_by(Spp_code, species,binned_ag, binned_year) %>%   
  summarise(mean_ag = mean(PercentAg),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_liberal)/sqrt(sample)) %>% 
  mutate(ag_bin = case_when(mean_ag>=40~ paste("20"),
                              mean_ag<40 ~ paste("80")))


histogram_data <- endo_herb_disturbance %>%
  filter(disturbance_ratio>-Inf) %>% 
  mutate(disturbance_label = case_when(disturbance_ratio > -9 & disturbance_ratio < -1.5 ~ "Urban",
                                       disturbance_ratio >= -1.5 & disturbance_ratio <= 1.5 ~ "Undisturbed",
                                       disturbance_ratio > 1.5 ~ "Agricultural"))
         
  


# year_trend <- ggplot(year.pred) +
#   # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = ag_bin))+
#   geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
#   geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
#   geom_line(aes(year, mean, group = disturbance_label, color = disturbance_label)) +
#   facet_wrap(~species)+
#   # scale_color_manual(values = c("#33a02c", "black", "#1f78b4"))+
#   # scale_fill_manual(values = c("#33a02c", "black", "#1f78b4"))+
#   labs(y = "Endophyte Prevalence", x = "Year", color = "Percent Ag.",fill = "Percent Ag.", size = "Sample Size")+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         legend.text = element_text(face = "italic"),
#         plot.margin = unit(c(0,.1,.1,.1), "line"))+
#   lims(y = c(0,1))
# 
# year_trend
# ggsave(year_trend, filename = "year_trend.png")

year_trend_facet <- ggplot(year.pred) +
  geom_dots(data = histogram_data, aes(x = year, y = Endo_status_liberal, color = disturbance_label, side = ifelse(Endo_status_liberal == 1, "bottom", "top")),binwidth = 2)+
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample, color = ag_bin))+
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = disturbance_label, fill = disturbance_label), alpha = 0.3) +
  geom_line(aes(year, mean, group = disturbance_label, color = disturbance_label)) +
  facet_wrap(~species+disturbance_label)+
  scale_color_manual(values = c("#b35806", "black", "#542788"))+
  scale_fill_manual(values = c("#b35806", "black", "#542788"))+
  guides(color = "none", fill = "none")+
  labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))

year_trend_facet

ggsave(year_trend_facet, filename = "year_trend_facet.png")













################################################################################################################################
##########  Plotting the posteriors from the model ###############
################################################################################################################################

param_names <- fit.2$summary.random$fixed$ID

n_draws <- 500

posteriors <- generate(
  fit.2,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))


posteriors_df.not_N <- as_tibble(posteriors, rownames = "param") %>% 
  mutate(param_label = case_when(param == param_names[1] ~ "Int AGHY",
                               param == param_names[2] ~ "Int AGPE",
                               param == param_names[3] ~ "Int ELVI",
                               param == param_names[4] ~ "Ag AGHY",
                               param == param_names[5] ~ "Year AGHY",
                               param == param_names[6] ~ "Urban AGHY",
                               param == param_names[7] ~ "Ag AGPE",
                               param == param_names[8] ~ "Ag ELVI",
                               param == param_names[9] ~ "Year AGPE",
                               param == param_names[10] ~ "Year ELVI",
                               param == param_names[11] ~ "Year:Ag AGHY",
                               param == param_names[12] ~ "Urban AGPE",
                               param == param_names[13] ~ "Urban ELVI",
                               param == param_names[14] ~ "Year:Urban AGHY",
                               param == param_names[15] ~ "Year:Ag AGPE",
                               param == param_names[16] ~ "Year:Ag ELVI",
                               param == param_names[17] ~ "Year:Urban AGPE",
                               param == param_names[18] ~ "Year:Urban ELVI"),
         spp_label = sub(".* ", "", param_label),
         param_type = sub(" .*", "", param_label)) %>% 
  pivot_longer( cols = -c(param, param_label, spp_label, param_type), names_to = "iteration") %>% 
  mutate(model = "No Nitrogen")

posteriors_summary <- posteriors_df.not_N %>% 
  group_by(param, param_label, spp_label, param_type) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .1),
            upr = quantile(value, .9))

posterior_hist <- ggplot(posteriors_df.not_N)+
  # geom_density(aes(x = value, fill = spp_label), color = "black", alpha = .6)+
  geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_type, scales = "free")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_without_nitrogen.png")


# another way to make the plot, which is nice, but it is a lot harder for me to control the layout of this plot compared to working with raw draws from the posteriors
# flist <- vector("list", NROW(fit$summary.random$fixed))
# 
# for (i in seq_along(flist)) flist[[i]] <- plot(fit, varname = "fixed", index = i) +geom_vline(aes(xintercept = 0))
# multiplot(plotlist = flist, cols = 3)
# 
# 
# 
# 
# # This is probably the graph to show the probability that temporal trends differ between more and less agricultural areas
# 
# 
# agXyear_posterior <- ggplot(filter(posteriors_df, grepl("PercentAg:std_year", posteriors_df$param)))+
#   geom_histogram(aes(x = value))+
#   geom_vline(xintercept = 0)+
#   facet_wrap(~param, scales = "free")
# agXyear_posterior
# 
# ggsave(agXyear_posterior, filename = "agXyear_posterior.png")












################################################################################################################################
##########  Now re-fitting the model to include nitrogen deposition as a predictor ###############
################################################################################################################################

# First making a plot of the relationship between ag and urb and nitrogen


disturbanceXnitrogen <- ggplot(endo_herb_disturbance)+
  geom_point(aes(x = PercentAg, y = NO3_mean, color = NO3_mean))+
  scale_color_viridis_c()+
  geom_point(aes(x = quantile(PercentAg, .05), y = quantile(NO3_mean, .95)), color = "#542788", size = 4, shape = 4)+
  geom_point(aes(y = quantile(NO3_mean, .05), x = quantile(PercentAg, .95)), color = "#b35806", size = 4, shape = 4)+
  geom_point(aes(y = quantile(NO3_mean, .05), x = quantile(PercentAg, .05)), color = "black", size = 4, shape = 4)+
  # labs(x = "Percent Ag. (%)", y = "Percent Urban (%)", color = "Nitrogen \n Deposition")+
  theme_bw()


disturbanceXnitrogen
ggsave(disturbanceXnitrogen, filename = "disturbanceXnitrogen_dimensions.png")
cor(endo_herb$PercentAg, endo_herb$PercentUrban)
cor(endo_herb$PercentAg, endo_herb$NO3_mean)
cor(endo_herb$PercentUrban, endo_herb$NO3_mean)


# version with all species in one model. Note that we remove the intercept, and then we have to specify that the species is in the fixed effects model

s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg*std_year + Spp_code*PercentUrban*std_year +  Spp_code*NO3_mean*std_year, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde) 

s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentAg + Spp_code*PercentUrban +  Spp_code*NO3_mean, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde) 
# s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*PercentUrban*std_year, model = "fixed")+
#   space_int(coords, model = spde) 
# formula, with "." meaning "add all the model components":
s_formula <- Endo_status_liberal ~ .


# Now run the model

fit.N <- bru(s_components,
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

fit.N$dic$dic
fit.N$mode$mode.status # a 0 or low value indicates "convergence"

fit.N$summary.fixed
fit.N$summary.random

saveRDS(fit.N, file = "fit_with_N.rds")
fit.2 <- read_rds( file = "fit_wo_N.rds")

fit.N <- read_rds(file = "fit_with_N.rds")




################################################################################################################################
##########  Getting and plotting prediction across nitrogen depostion from the model  ###############
################################################################################################################################

min_nit<- min(data$NO3_mean)
max_nit <- max(data$NO3_mean)
mean_nit <- mean(data$NO3_mean)



min_year <- min(data$std_year)
max_year <- max(data$std_year)

mean_year <- mean(data$year)


min_ag<- quantile(data$PercentAg, .05)
max_ag <- quantile(data$PercentAg, .95)

min_urb <- quantile(data$PercentUrban, .05)
max_urb <- quantile(data$PercentUrban, .95)



preddata <- expand.grid(Spp_code = c("AGHY","AGPE","ELVI"), #std_year = c(1900, 1980)-mean_year, 
                        PercentAg = seq(min_ag, max_ag,  length.out = 20), PercentUrban = seq(min_urb, max_urb,  length.out = 20), NO3_mean = seq(min_nit, max_nit, length.out = 20),
                        
                        collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3])) %>% 
  filter(PercentAg == min_ag & PercentUrban == max_urb | PercentAg == max_ag & PercentUrban == min_urb | PercentAg == min_ag & PercentUrban == min_urb) %>% 
  mutate(disturbance_label = case_when(PercentAg == min_ag & PercentUrban == min_urb ~ "Undisturbed",
                                       PercentAg == max_ag ~ "Agricultural",
                                       PercentUrban == max_urb ~ "Urban"))


# gennerating predictions and back-transforming the standardized year variable


# fit.N <- fit.4
nit.pred <- predict(
  fit.N,
  newdata = preddata,
  formula = ~ invlogit(fixed),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) #%>% 
  #mutate(year = std_year + mean_year)


# binning the data for plotting
endo_herb_binned_nitrogen <- endo_herb_disturbance %>%
  filter(disturbance_ratio>-Inf) %>% 
  mutate(binned_nit = cut(NO3_mean, breaks = 8),
         binned_disturb = cut(disturbance_ratio, breaks = 3),
         #binned_year = cut(year, breaks = 2)
         ) %>% 
  group_by(Spp_code, species,binned_nit, binned_disturb) %>% #, binned_year) %>%
  summarise(mean_nit = mean(NO3_mean),
            mean_disturb = mean(disturbance_ratio),
            #mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n(),
            se_endo = sd(Endo_status_liberal)/sqrt(sample)) %>% 
  mutate(disturbance_label = case_when(mean_disturb > -9 & mean_disturb < -1.5 ~ "Urban",
                                       mean_disturb >= -1.5& mean_disturb <= 1.5 ~ "Undisturbed",
                                       mean_disturb > 1.5 ~ "Agricultural"),
         
         #year_bin = case_when(mean_year>=1940~ paste("1980"),
          #                    mean_year<1940 ~ paste("1900"))
)

nit_trend <- ggplot(nit.pred) +
  geom_point(data = endo_herb_binned_nitrogen, aes(x = mean_nit, y = mean_endo, size = sample))+
  geom_line(aes(NO3_mean, mean, )) +
  geom_ribbon(aes(NO3_mean, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(NO3_mean, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  facet_wrap(~species + disturbance_label)+
  scale_color_manual(values = c("blue", "red"))+
  scale_fill_manual(values = c("blue", "red"))+
  
  labs(y = "Endophyte Prevalence", x = "Percent Nitrogen (%)", color = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.text = element_text(face = "italic"),
        plot.margin = unit(c(0,.1,.1,.1), "line"))+
  lims(y = c(0,1))

# nit_trend <- ggplot(nit.pred) +
#   geom_point(data = endo_herb_binned_nitrogen, aes(x = mean_nit, y = mean_endo, size = sample, color = year_bin))+
#   geom_line(aes(NO3_mean, mean, group = as.factor(year), color = as.factor(year))) +
#   geom_ribbon(aes(NO3_mean, ymin = q0.025, ymax = q0.975, group = as.factor(year), fill = as.factor(year)), alpha = 0.2) +
#   geom_ribbon(aes(NO3_mean, ymin = q0.25, ymax = q0.75,group = as.factor(year), fill = as.factor(year)), alpha = 0.2) +
#   facet_wrap(~species + disturbance_label)+
#   scale_color_manual(values = c("blue", "red"))+
#   scale_fill_manual(values = c("blue", "red"))+
#   
#   labs(y = "Endophyte Prevalence", x = "Percent Nitrogen (%)", color = "Year", size = "Sample Size")+
#   theme_classic()+
#   theme(strip.background = element_blank(),
#         legend.text = element_text(face = "italic"),
#         plot.margin = unit(c(0,.1,.1,.1), "line"))+
#   lims(y = c(0,1))

nit_trend

ggsave(nit_trend, filename = "model_with_nitrogen.png")








################################################################################################################################
##########  Plotting the posteriors from the model with nitrogen ###############
################################################################################################################################

param_names <- fit.N$summary.random$fixed$ID

n_draws <- 500

posteriors.N <- generate(
  fit.N,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors.N) <- param_names
colnames(posteriors.N) <- c( paste0("iter",1:n_draws))


posteriors_df.N <- as_tibble(posteriors.N, rownames = "param") %>% 
  mutate(param_label = case_when(param == param_names[1] ~ "Int AGHY",
                                 param == param_names[2] ~ "Int AGPE",
                                 param == param_names[3] ~ "Int ELVI",
                                 param == param_names[4] ~ "Ag AGHY",
                                 param == param_names[5] ~ "Year AGHY",
                                 param == param_names[6] ~ "Nit AGHY",
                                 param == param_names[7] ~ "Urban AGHY",
                                 param == param_names[8] ~  "Ag AGPE",
                                 param == param_names[9] ~  "Ag ELVI",
                                 param == param_names[10] ~ "Year AGPE",
                                 param == param_names[11] ~ "Year ELVI",
                                 param == param_names[12] ~ "Year:Ag AGHY",
                                 param == param_names[13] ~ "Nit AGPE",
                                 param == param_names[14] ~ "Nit ELVI",
                                 param == param_names[15] ~ "Ag:Nit AGHY",
                                 param == param_names[16] ~ "Year:Nit AGHY",
                                 param == param_names[17] ~  "Urban AGPE",
                                 param == param_names[18] ~ "Urban ELVI",
                                 param == param_names[19] ~ "Year:Urban AGHY",
                                 param == param_names[20] ~ "Urban:Nit AGHY",
                                 param == param_names[21] ~ "Year:Ag AGPE",
                                 param == param_names[22] ~ "Year:Ag ELVI",
                                 param == param_names[23] ~ "Ag:Nit AGPE",
                                 param == param_names[24] ~ "Ag:Nit ELVI",
                                 param == param_names[25] ~ "Year:Nit AGPE",
                                 param == param_names[26] ~ "Year:Nit ELVI",
                                 param == param_names[27] ~ "Year:Ag:Nit AGHY",
                                 param == param_names[28] ~ "Year:Urban AGPE",
                                 param == param_names[29] ~ "Year:Urban ELVI",
                                 param == param_names[30] ~ "Urban:Nit AGPE",
                                 param == param_names[31] ~ "Urban:Nit ELVI",
                                 param == param_names[32] ~ "Year:Urban:Nit AGHY",
                                 param == param_names[33] ~ "Year:Ag:Nit AGPE",
                                 param == param_names[34] ~ "Year:Ag:Nit ELVI",
                                 param == param_names[35] ~ "Year:Urban:Nit AGPE",
                                 param == param_names[36] ~ "Year:Urban:Nit ELVI"),
         spp_label = sub(".* ", "", param_label),
         param_type = sub(" .*", "", param_label)) %>% 
  pivot_longer( cols = -c(param, param_label, spp_label, param_type), names_to = "iteration") %>% 
  mutate(model = "Nitrogen")


posteriors_df <- posteriors_df.not_N %>% full_join(posteriors_df.N)

posteriors_summary <- posteriors_df %>% 
  group_by(param, param_label, spp_label, param_type, model) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .1),
            upr = quantile(value, .9))

posterior_hist <- ggplot(posteriors_df)+
  # geom_density(aes(x = value, fill = spp_label), color = "black", alpha = .6)+
  geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = model), position = position_dodge2(w = 0.75))+
  geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = model), position = position_dodge2(w = 0.75))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_type, scales = "free")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw()

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_with_nitrogen.png")


# another way to make the plot, which is nice, but it is a lot harder for me to control the layout of this plot compared to working with raw draws from the posteriors
flist <- vector("list", NROW(fit$summary.random$fixed))

for (i in seq_along(flist)) flist[[i]] <- plot(fit, varname = "fixed", index = i) +geom_vline(aes(xintercept = 0))
multiplot(plotlist = flist, cols = 3)




# This is probably the graph to show the probability that temporal trends differ between more and less agricultural areas


agXyear_posterior <- ggplot(filter(posteriors_df, grepl("PercentAg:std_year:NO3_mean", posteriors_df$param)))+
  geom_histogram(aes(x = value))+
  geom_vline(xintercept = 0)+
  facet_wrap(~param, scales = "free")
agXyear_posterior

ggsave(agXyear_posterior, filename = "agXyear_posterior.png")

NO3_mean_posterior <- ggplot(filter(posteriors_df, grepl("NO3_mean", posteriors_df$param)))+
  geom_histogram(aes(x = value))+
  geom_vline(xintercept = 0)+
  facet_wrap(~param, scales = "free")
agXyear_posterior

ggsave(agXyear_posterior, filename = "agXyear_posterior.png")
















# Mapping the coeffients

# get easting and northing limits

xlim <- range(mesh$loc[, 1])
ylim <- range(mesh$loc[, 2])
grd_dims <- round(c(x = diff(range(xlim)), y = diff(range(ylim))))

# # make mesh projector to get model summaries from the mesh to the mapping grid
# mesh_proj <- fm_evaluator(
#   mesh,
#   xlim = xlim, ylim = ylim, dims = grd_dims
# )
# 
# 
# 
# space_int <- data.frame(
#   median = invlogit(fit$summary.random$space_int$"0.5quant"),
#   range95 = invlogit(fit$summary.random$space_int$"0.975quant") -
#     invlogit(fit$summary.random$space_int$"0.025quant")
# )
# 
# 
# 
# # loop to get estimates on a mapping grid
# pred_grids <- lapply(
#   list(space_int = space_int),
#  fun <- function(x) as.matrix(fm_evaluate(mesh_proj, x))
# )
# 
# 
# # make a terra raster stack with the posterior median and range95
# out_stk <- rast()
# for (j in 1:1) {
#   mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
#                   Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
#   )
#   mean_j <- rast(mean_j)
#   range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
#                      Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
#   )
#   range95_j <- rast(range95_j)
#   out_j <- c(mean_j, range95_j)
#   terra::add(out_stk) <- out_j
# }
# names(out_stk) <- c(
#   "space_median", "space_range95"
# )
# #Masking the raster to our boundary
# out_stk <- terra::mask(out_stk,bdry_st)
# 
# make_plot_field <- function(data_stk, scale_label) {
#   ggplot() +
#     geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
#     geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
#     geom_sf(fill = NA) +
#     coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
#     geom_spatraster(data = data_stk) +
#     labs(x = "", y = "") +
#     scale_fill_viridis_c(name = scale_label, option = "turbo", na.value = "transparent")+
#     theme(text = element_text(size = 2))+
#     theme_bw()
# }
# 
# 
# 
# # This is the plot of the spatial intercept
# ps <- make_plot_field(
#   data_stk = out_stk[["space_median"]],
#   scale_label = "Spatial\nPosterior Mean"
# )
# ps
# 
# ps_r <- make_plot_field(
#   data_stk = out_stk[["space_range95"]],
#   scale_label = "Spatial\nPosterior 95 CI\n"
# )
# ps_r
# 
# # Taking alot of material from this blog post: https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html
# 
# # plotting the spde range and sd posteriors
# spde.range <- spde.posterior(fit, "space_int", what = "range")
# spde.logvar <- spde.posterior(fit, "space_int", what = "log.variance")
# spde.var <- spde.posterior(fit, "space_int", what = "variance")
# 
# range.plot <- plot(spde.range)
# var.plot <- plot(spde.var)
# 
# 
# range.plot
# var.plot
# 
# # and plot the matern covariance (our spatial decay effect)
# 
# cov.plot <- plot(spde.posterior(fit, "space_int", what = "matern.covariance"))
# cov.plot
# 
# 
# 
# # Making a plot of the marginal posterior of the year slope. We can see the posterior has a small positive slope for AGHY
# flist <- vector("list", NROW(fit$summary.fixed))
# for (i in seq_along(flist)) {
#   flist[[i]] <- plot(fit, rownames(fit$summary.fixed)[i])
# }
# multiplot(plotlist = flist, cols = 2)
# 
# 
# 
# ###### Getting and plotting prediction from model #####
# # Here I am showing u
# min_year <- min(endo_herb$year)
# max_year <- max(endo_herb$year)
# 
# year.pred <- predict(
#   fit,
#   data.frame(year = seq(min_year, max_year)),
#   formula = ~ invlogit(Intercept + year))
# 
# 
# ggplot(year.pred) +
#   geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal))+
#   geom_line(aes(year, mean)) +
#   geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
#   geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
#   lims(y = c(0,1))
# 
############################################################################################### 
########################################### Model with Ag*yr, Urb*yr, and Nit*yr############################
#######################################################################################################
a_components <- ~ Intercept(1) +
  NO3_mean + yearEffect(year, model = "linear") + yearXnitEffect(year*NO3_mean, model = "linear") +
  space_int(coords, model = spde) 

a_formula <- Endo_statu ~ .

a_fit <- bru(a_components,
             like(
               formula = a_formula,
               family = "binomial",
               Ntrials = 1,
               data = endo_herb
             ),
             options = list(
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
               control.inla = list(int.strategy = "eb"),
               verbose = TRUE
             )
)



###### Getting and plotting prediction from model #####
# Here I am showing u
min_urb <- min(endo_herb$PercentUrban, na.rm = TRUE)
max_urb <- max(endo_herb$PercentUrban, na.rm = TRUE)
preddata <- data.frame(PercentUrban = seq(min_urb, max_urb))
#preddata <- data.frame(PercentUrban = seq(min_urb, max_urb))


nit.pred <- predict(
  fit,
  newdata = preddata,
  formula = ~ invlogit(Intercept + year + yearXnitEffect))


# binning the data for plotting
endo_herb_binned <- endo_herb %>% 
  mutate(binned_year = cut(year, breaks = 50)) %>%
  group_by(Spp_code, species,binned_year) %>%   
  summarise(mean_year = mean(year),
            mean_endo = mean(Endo_statu),
            mean_lon = mean(lon),
            mean_lat = mean(lat),
            sample = n()) %>% 
  mutate(lat_bin = case_when(mean_lat>=35 ~ paste("43"),
                             mean_lat<35 ~ paste("35")),
         lon_bin = case_when(mean_lon<=-94 ~ paste("-90"),
                             mean_lon>-94 ~ paste("-80") ))


ggplot(urb.pred) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  # geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(PercentUrban, mean)) +
  geom_ribbon(aes(PercentUrban, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(PercentUrban, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  lims(y = c(0,1))


endo_herb <- endo_herb %>% mutate_at(c('PercentUrban', 'PercentAg', 'NO3_mean', "year"), ~(scale(.) %>% as.vector))

min_ag <- min(endo_herb$PercentAg, na.rm = TRUE)
max_ag <- max(endo_herb$PercentAg, na.rm = TRUE)
min_nit <- min(endo_herb$NO3_mean, na.rm = TRUE)
max_nit <- max(endo_herb$NO3_mean, na.rm = TRUE)

all_preddata <- expand.grid(PercentUrban = seq(min_urb, max_urb), PercentAg = mean(endo_herb$PercentAg, na.rm = TRUE), NO3_mean = c(min_nit, max_nit), year = NA)

#all_preddata <- data.frame(PercentUrban = seq(min_urb, max_urb), PercentAg = seq(min_ag, max_ag, length.out = 100), NO3_mean = seq(min_nit, max_nit, length.out =100), = seq(min_year, max_year, length.out = 100))

all.pred <- predict(
  fit,
  newdata = all_preddata,
  formula = ~ invlogit(Intercept +PercentUrban + NO3_mean + year))


ggplot(all.pred) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  # geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(PercentUrban, mean, color = NO3_mean, group = NO3_mean)) +
 # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(PercentUrban, ymin = mean - 1 * sd, ymax = mean + 1 * sd, group = NO3_mean, fill = NO3_mean), alpha = 0.2) 
  lims(y = c(0,1))

  

#######################################################################################################################
################################################## straightening this up ##############################################
#######################################################################################################################
#creating full model
#and then model with NO3_mean and ag together


#4.22.24, fixed effect thing josh sent me  
  #removed: + PercentAg*NO3_mean,
s_components <- ~ my_fixed_effects(
  main = PercentAg*year*NO3_mean,
  model = "fixed")+
  space_int(coords, model = spde) 

  

# formula, with "." meaning "add all the model components":
s_formula <- Endo_statu ~ .


# Now run the model

fit <- bru(s_components,
           like(
             formula = s_formula,
             family = "binomial",
             Ntrials = 1,
             data = endo_herb
           ),
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE
           )
)

#now plot this full model
#Here, make the sequences and data frame for each variable, can use these to make prediction data frame, can also hopefully use for later predictions too?
min_urb <- min(endo_herb$PercentUrban, na.rm = TRUE)
max_urb <- max(endo_herb$PercentUrban, na.rm = TRUE)

min_ag <- min(endo_herb$PercentAg, na.rm = TRUE)
max_ag <- max(endo_herb$PercentAg, na.rm = TRUE)

min_nit <- min(endo_herb$NO3_mean, na.rm = TRUE)
max_nit <- max(endo_herb$NO3_mean, na.rm = TRUE)

#prediction dataframe. removed: NO3_mean = seq(min_nit, max_nit)
allpreddata <- expand.grid(PercentAg = seq(min_ag, max_ag), NO3_mean= c(min_nit, max_nit), year = c(1900,2000))
#which goes into this:
allpreddata <- predict(
  fit,
  newdata = allpreddata,
  formula = ~ invlogit(Intercept + my_fixed_effects)) 
allpreddata <- allpreddata %>%
  mutate(NO3_mean = as.factor(NO3_mean), year = as.factor(year))
  

#now plot:
ggplot(allpreddata) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
 # geom_point(data =endo_herb, aes(x = PercentAg, y = Endo_statu), shape = "|")+
  geom_line(aes(PercentAg, mean, color = as.factor(year))) +
  # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975, fill = as.factor(year), group = year), alpha = 0.2)+
  facet_wrap(~ NO3_mean)
#lims(y = c(0,1))
  
################################## Ag*yr model ###########################
# b_components <- 
#   fixed_effects(main = ~ spp_code + year +PercentAg:year,) +
#   space_int(coords, model = spde) 
#tried a few different ways to make fixed effects thing but did not yet succeed

b_components <- ~ Intercept(1) +
  PercentAg:year +
  space_int(coords, model = spde) 



# formula, with "." meaning "add all the model components":
b_formula <- Endo_statu ~ .


# Now run the model

b_fit <- bru(b_components,
           like(
             formula = b_formula,
             family = "binomial",
             Ntrials = 1,
             data = endo_herb
           ),
           options = list(
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.inla = list(int.strategy = "eb"),
             verbose = TRUE
           )
)


#prediction dataframe
agpreddata <- expand.grid(Agseq = seq(min_ag, max_ag), year = c(1900,2000))
#which goes into this:
agpreddata <- predict(
  fit,
  newdata = agpreddata,
  formula = ~ invlogit(Intercept + year + Agseq)) 

#now plot:
#will troubleshoot later
ggplot(agpreddata) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(NO3_mean, mean, color = year, group = year)) +
  # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(NO3_mean, ymin = mean - 1 * sd, ymax = mean + 1 * sd, group = NO3_mean, fill = NO3_mean), alpha = 0.2)+
  lims(y = c(0,1), x = c(1900,2000))

################################## Ag*yr + ag*nit  model ###########################

c_components <- ~ Intercept(1) +
  PercentAg:year + PercentAg*NO3_mean +
  space_int(coords, model = spde) 



# formula, with "." meaning "add all the model components":
c_formula <- Endo_statu ~ .


# Now run the model

c_fit <- bru(c_components,
             like(
               formula = c_formula,
               family = "binomial",
               Ntrials = 1,
               data = endo_herb
             ),
             options = list(
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
               control.inla = list(int.strategy = "eb"),
               verbose = TRUE
             )
)


#prediction dataframe
agnitpreddata <- expand.grid(Agseq = seq(min_ag, max_ag), year = c(1900,2000),NO3_mean = seq(min_nit, max_nit))
#which goes into this:
agnitpreddata <- predict(
  fit,
  newdata = agnitpreddata,
  formula = ~ invlogit(Intercept + year + Agseq)) 

#now plot:
#will troubleshoot later
ggplot(agnitpreddata) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(NO3_mean, mean, color = year, group = year)) +
  # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(NO3_mean, ymin = mean - 1 * sd, ymax = mean + 1 * sd, group = NO3_mean, fill = NO3_mean), alpha = 0.2)+
  lims(y = c(0,1), x = c(1900,2000))

################################## nit*yr model ###########################
d_components <- ~ Intercept(1) +
  NO3_mean:year +
  space_int(coords, model = spde) 



# formula, with "." meaning "add all the model components":
d_formula <- Endo_statu ~ .


# Now run the model

d_fit <- bru(d_components,
             like(
               formula = d_formula,
               family = "binomial",
               Ntrials = 1,
               data = endo_herb
             ),
             options = list(
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
               control.inla = list(int.strategy = "eb"),
               verbose = TRUE
             )
)


#prediction dataframe
nitpreddata <- expand.grid(year = c(1900,2000), NO3_mean = seq(min_nit, max_nit))
#which goes into this:
nitpreddata <- predict(
  fit,
  newdata = agpreddata,
  formula = ~ invlogit(Intercept + year + NO3_mean)) 

#now plot:
#will troubleshoot later
ggplot(nitpreddata) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(NO3_mean, mean, color = year, group = year)) +
  # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(NO3_mean, ymin = mean - 1 * sd, ymax = mean + 1 * sd, group = NO3_mean, fill = NO3_mean), alpha = 0.2)+
  lims(y = c(0,1), x = c(1900,2000))

