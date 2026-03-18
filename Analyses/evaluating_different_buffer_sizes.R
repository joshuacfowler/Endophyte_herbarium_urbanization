# Purpose: Fits endophyte models using NO3, NH4 and TIN as alternate predictors
# Authors: Mallory Tucker and Joshua Fowler
# Updated: May 21, 2025
# Purpose: Fits endophyte models using liberal and conservative endophyte scores
# Authors: Mallory Tucker and Joshua Fowler
# Updated: June 24, 2025


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
library(alphahull)

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

# adding in the 30km buffer land cover data
landcover_30km <- read_csv(file = paste0(path,"full_Zonalhist_NLCD_2001_30km.csv")) %>% 
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
         TotalPixels = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24 + HISTO_NODATA + HISTO_11 +HISTO_12 + HISTO_31 +HISTO_41 + HISTO_42 + HISTO_43 + HISTO_52 + HISTO_71+ HISTO_90 + TotalAg + HISTO_95,
         TotalDeveloped = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24,
         OtherLC = (TotalPixels - (TotalAg + TotalDeveloped))/TotalPixels *100,
         PercentUrban_30km = TotalDeveloped/TotalPixels * 100,
         PercentAg_30km = TotalAg/TotalPixels * 100) %>% 
  select(Sample_id, Institution_specimen_id, Spp_code, new_id, score_number, year, month, day, lat, lon, PercentUrban_30km, PercentAg_30km)


endo_herb_georef <- endo_herb_georef %>% left_join(landcover_30km)

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

endo_herb <- left_join(endo_herb_sf, climate, by = c("lon", "lat", "year")) %>% 
  filter(year>=1895) %>% 
  filter(seed_scored != 0) %>% 
  filter(Sample_id != "BRIT_AGHY_506") %>% # dropping this because the georeferenceing is wrong
  filter(!(Sample_id == "AM_ELVI_122" & is.na(Municipality))) %>% 
  dplyr::distinct(Sample_id, score_number, .keep_all = TRUE) %>% 
  mutate(Endo_status_conservative = case_when(Endo_status_conservative == 4 ~ 1, TRUE ~ Endo_status_conservative))



#load in map outline data
outline_map <- map_data("world")
states_shape <- map_data("state")
counties <- map_data("county")


##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
data_summary <- endo_herb %>% 
  dplyr::summarize(year = mean(year, na.rm = T),
                   PercentUrban = mean(PercentUrban, na.rm = T),
                   PercentUrban_30km = mean(PercentUrban_30km, na.rm = T),
                   PercentAg = mean(PercentAg, na.rm = T),
                   PercentAg_30km = mean(PercentAg_30km, na.rm = T),
                   mean_TIN_10km = mean(mean_TIN_10km, na.rm = T),
                   mean_TIN_30km = mean(mean_TIN_30km, na.rm = T),
                   tmean_10km = mean(tmean_10km, na.rm = T),
                   tmean_30km = mean(tmean_30km, na.rm = T),
                   ppt_10km = mean(ppt_10km, na.rm = T),
                   ppt_30km = mean(ppt_30km, na.rm = T))

data <- endo_herb %>% 
  mutate(year = year - data_summary$year,
         PercentUrban = PercentUrban - data_summary$PercentUrban,
         PercentAg = PercentAg - data_summary$PercentAg,
         mean_TIN_10km = mean_TIN_10km - data_summary$mean_TIN_10km,
         tmean_10km = tmean_10km - data_summary$tmean_10km,
         ppt_10km = ppt_10km - data_summary$ppt_10km,
         PercentUrban_30km = PercentUrban_30km - data_summary$PercentUrban_30km,
         PercentAg_30km = PercentAg_30km - data_summary$PercentAg_30km,
         mean_TIN_30km = mean_TIN_30km - data_summary$mean_TIN_30km,
         tmean_30km = tmean_30km - data_summary$tmean_30km,
         ppt_30km = ppt_30km - data_summary$ppt_30km) %>%  
  mutate(Spp_index = as.numeric(as.factor(Spp_code))) %>% 
  filter(!is.na(ppt_10km))

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




s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(mean_TIN_10km + PercentAg + PercentUrban + ppt_10km)^4, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



s_components.year <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(mean_TIN_10km + PercentAg + PercentUrban + year + ppt_10km)^2, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



s_components_30km <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(mean_TIN_30km + PercentAg_30km + PercentUrban_30km + ppt_30km)^4, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



s_components.year_30km <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(mean_TIN_30km + PercentAg_30km + PercentUrban_30km + year + ppt_30km)^2, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)




# putting the components with the formula
s_formula <- Endo_status_liberal ~ .


# Now run the model

fit.10km <- bru(s_components,
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




fit.30km <- bru(s_components_30km,
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

fit.10km.year <- bru(s_components.year,
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

fit.30km.year <- bru(s_components.year_30km,
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
fit.10km$dic$dic
fit.30km$dic$dic
fit.10km.year$dic$dic
fit.30km.year$dic$dic


# checking convergence
fit.10km$mode$mode.status
fit.30km$mode$mode.status
fit.10km.year$mode$mode.status
fit.30km.year$mode$mode.status




################################################################################################################################
##########  Plotting the posteriors from the model without year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
param_names.10km <- fit.10km$summary.random$fixed$ID
param_names.30km <- fit.30km$summary.random$fixed$ID

n_draws <- 1000

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors.10km <- generate(
  fit.10km,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors.10km) <- param_names.10km
colnames(posteriors.10km) <- c( paste0("iter",1:n_draws))

posteriors.30km <- generate(
  fit.30km,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors.30km) <- param_names.30km
colnames(posteriors.30km) <- c( paste0("iter",1:n_draws))


posteriors.10km_df <- as_tibble(t(posteriors.10km), rownames = "iteration")
colnames(posteriors.10km_df) <- sub("Spp_code", "", colnames(posteriors.10km_df))
colnames(posteriors.10km_df) <- gsub(":", ".", colnames(posteriors.10km_df))

posteriors.30km_df <- as_tibble(t(posteriors.30km), rownames = "iteration")
colnames(posteriors.30km_df) <- sub("Spp_code", "", colnames(posteriors.30km_df))
colnames(posteriors.30km_df) <- gsub(":", ".", colnames(posteriors.30km_df))


# Calculate the effects of the predictor, given that the reference level is for AGHY
effects.10km <- posteriors.10km_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year",
         data = "10km") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "PercentAg" = "Agr.",
                                                         "PercentUrban" = "Urb.",
                                                         "mean_TIN_10km" = "Nit.",
                                                         "ppt_10km" = "PPT.")),
                          levels = c("Intercept","Nit.","Agr.","Urb.","PPT.","Nit. X Agr.","Nit. X Urb.","Nit. X PPT.","Agr. X Urb.","Agr. X PPT.","Urb. X PPT.","Nit. X Agr. X Urb.","Nit. X Agr. X PPT.","Nit. X Urb. X PPT.","Agr. X Urb. X PPT.","Nit. X Agr. X Urb. X PPT.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))

effects.30km <- posteriors.30km_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year",
         data = "30km") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "PercentAg_30km" = "Agr.",
                                                         "PercentUrban_30km" = "Urb.",
                                                         "mean_TIN_30km" = "Nit.",
                                                         "ppt_30km" = "PPT.")),
                          levels = c("Intercept","Nit.","Agr.","Urb.","PPT.","Nit. X Agr.","Nit. X Urb.","Nit. X PPT.","Agr. X Urb.","Agr. X PPT.","Urb. X PPT.","Nit. X Agr. X Urb.","Nit. X Agr. X PPT.","Nit. X Urb. X PPT.","Agr. X Urb. X PPT.","Nit. X Agr. X Urb. X PPT.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))



effects_df <- bind_rows(effects.10km, effects.30km)

buffer_colors <- c("#808080", "#2ca25f")


posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = data, point_color = data, group =  data), breaks = 50, normalize = "panels", alpha = .4)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", ncol = 4)+
  labs(x = "Posterior Est.", y = "Species", fill = "Buffer Radius", point_color = "Buffer Radius")+
  scale_color_manual(values = buffer_colors, aesthetics = "point_color")+
  scale_fill_manual(values = buffer_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/different_buffer_posteriors.png", width = 9, height = 9)





################################################################################################################################
##########  Plotting the posteriors from the model with year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
param_names.year.10km <- fit.10km.year$summary.random$fixed$ID
param_names.year.30km <- fit.30km.year$summary.random$fixed$ID

n_draws <- 1000

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors.year.10km <- generate(
  fit.10km.year,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors.year.10km) <- param_names.year.10km
colnames(posteriors.year.10km) <- c( paste0("iter",1:n_draws))

posteriors.year.30km <- generate(
  fit.30km.year,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors.year.30km) <- param_names.year.30km
colnames(posteriors.year.30km) <- c( paste0("iter",1:n_draws))


posteriors.year.10km_df <- as_tibble(t(posteriors.year.10km), rownames = "iteration")
colnames(posteriors.year.10km_df) <- sub("Spp_code", "", colnames(posteriors.year.10km_df))
colnames(posteriors.year.10km_df) <- gsub(":", ".", colnames(posteriors.year.10km_df))

posteriors.year.30km_df <- as_tibble(t(posteriors.year.30km), rownames = "iteration")
colnames(posteriors.year.30km_df) <- sub("Spp_code", "", colnames(posteriors.year.30km_df))
colnames(posteriors.year.30km_df) <- gsub(":", ".", colnames(posteriors.year.30km_df))


# Calculate the effects of the predictor, given that the reference level is for AGHY
effects.year.10km <- posteriors.year.10km_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "Year",
         data = "10km") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "year" = "Year",
                                                         "PercentAg" = "Agr.",
                                                         "PercentUrban" = "Urb.",
                                                         "mean_TIN_10km" = "Nit.",
                                                         "ppt_10km" = "Ppt.")),
                          levels = c("Intercept"  ,"Nit.","Agr.","Urb.","Year","Ppt.","Nit. X Agr.","Nit. X Urb.", "Nit. X Year","Nit. X Ppt.","Agr. X Urb.","Agr. X Year","Agr. X Ppt.","Urb. X Year","Urb. X Ppt.","Year X Ppt.","Nit. X Agr. X Urb.","Nit. X Agr. X Year","Nit. X Agr. X Ppt.","Nit. X Urb. X Year","Nit. X Urb. X Ppt.","Nit. X Year X Ppt.","Agr. X Urb. X Year","Agr. X Urb. X Ppt.","Agr. X Year X Ppt.","Urb. X Year X Ppt.","Nit. X Agr. X Urb. X Year","Nit. X Agr. X Urb. X Ppt.","Nit. X Agr. X Year X Ppt." ,"Nit. X Urb. X Year X Ppt.","Agr. X Urb. X Year X Ppt." )),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))

effects.year.30km <- posteriors.year.30km_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "Year",
         data = "30km") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "year" = "Year",
                                                         "PercentAg_30km" = "Agr.",
                                                         "PercentUrban_30km" = "Urb.",
                                                         "mean_TIN_30km" = "Nit.",
                                                         "ppt_30km" = "Ppt.")),
                          levels = c("Intercept"  ,"Nit.","Agr.","Urb.","Year","Ppt.","Nit. X Agr.","Nit. X Urb.", "Nit. X Year","Nit. X Ppt.","Agr. X Urb.","Agr. X Year","Agr. X Ppt.","Urb. X Year","Urb. X Ppt.","Year X Ppt.","Nit. X Agr. X Urb.","Nit. X Agr. X Year","Nit. X Agr. X Ppt.","Nit. X Urb. X Year","Nit. X Urb. X Ppt.","Nit. X Year X Ppt.","Agr. X Urb. X Year","Agr. X Urb. X Ppt.","Agr. X Year X Ppt.","Urb. X Year X Ppt.","Nit. X Agr. X Urb. X Year","Nit. X Agr. X Urb. X Ppt.","Nit. X Agr. X Year X Ppt." ,"Nit. X Urb. X Year X Ppt.","Agr. X Urb. X Year X Ppt." )),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))




effects_df <- bind_rows(effects.year.10km, effects.year.30km)





posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = data, point_color = data, group =  data), breaks = 50, normalize = "panels", alpha = .4)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", ncol = 4)+
  labs(x = "Posterior Est.", y = "Species", fill = "Buffer Radius", point_color = "Buffer Radius")+
  scale_color_manual(values = buffer_colors, aesthetics = "point_color")+
  scale_fill_manual(values = buffer_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/different_buffer_yearmodel_posteriors.png", width = 9, height = 11)

effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)
# write.csv(effects_summary, file = "Posterior_prob_results.csv")



