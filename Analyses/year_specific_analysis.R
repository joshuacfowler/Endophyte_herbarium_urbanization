# Purpose: Analysis using year Specific co-variate data
# Authors: Mallory Tucker and Joshua Fowler
# Updated: Dec 18, 2025


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
  dplyr::distinct(Sample_id, score_number, .keep_all = TRUE)



contemp_endo_herb <- endo_herb %>% 
  filter(!is.na(spec_PercentAg)) 

#################################################################################################
############ plotting relationship between the contemporary data and the "normal" data ##########
#################################################################################################

contemp_R2 <- contemp_endo_herb %>% 
  as_tibble() %>% select(-geometry) %>% 
  summarize(ag_R2 = cor(spec_PercentAg, PercentAg),
            urb_R2 = cor(spec_PercentUrban, PercentUrban))

TIN_R2 <- contemp_endo_herb %>% 
  as_tibble() %>% select(-geometry) %>% 
  filter(!is.na(TIN_10km)) %>% 
  summarize(TIN_R2 = cor(TIN_10km, mean_TIN_10km))

ag_contemp_R2 <- ggplot(contemp_endo_herb)+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x = spec_PercentAg, y = PercentAg), alpha = .6)+
  labs(x = "Agr. % (year specific)", y = "Agr. % (Normal)")+
  lims(x = c(0,100), y = c(0,100))+
  coord_equal()+
  labs(title = paste0("R^2", " = ", round(contemp_R2$ag_R2, digits = 3)))+
  theme_classic()
urb_contemp_R2 <- ggplot(contemp_endo_herb )+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x = spec_PercentUrban, y = PercentUrban), alpha = .6)+
  labs(x = "Urb. % (year specific)", y = "Urb. % (Normal)")+
  lims(x = c(0,100), y = c(0,100))+
  coord_equal()+
  labs(title = paste0("R^2", " = ", round(contemp_R2$urb_R2, digits = 3)))+
  theme_classic()
nit_contemp_R2 <- ggplot(contemp_endo_herb %>%   filter(!is.na(TIN_10km)))+
  geom_abline(intercept = 0, slope = 1)+
  geom_point(aes(x = TIN_10km, y = mean_TIN_10km), alpha = .6)+
  labs(x = "TIN (year specific)", y = "TIN (Normal)")+
  lims(x = c(150,900), y = c(150,900))+
  coord_equal()+
  labs(title = paste0("R^2", " = ", round(TIN_R2$TIN_R2, digits = 3)))+
  theme_classic()

contemp_R2_plot <- ag_contemp_R2 + urb_contemp_R2 + nit_contemp_R2 + plot_annotation(tag_levels = "A")
ggsave(contemp_R2_plot, filename = "Plots/contemp_R2.png", width = 8, height = 4)
##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
data_summary <- contemp_endo_herb %>% 
  filter(!is.na(ppt_10km)) %>% 
  filter(!is.na(TIN_10km)) %>% 
  dplyr::summarize(year = mean(year, na.rm = T),
                   spec_PercentAg = mean(spec_PercentAg, na.rm = T),
                   spec_PercentUrban = mean(spec_PercentUrban, na.rm = T),
                   TIN_10km = mean(mean_TIN_10km, na.rm = T),
                   tmean_10km = mean(tmean_10km, na.rm = T),
                   ppt_10km = mean(ppt_10km, na.rm = T))

data <- contemp_endo_herb %>% 
  filter(!is.na(ppt_10km)) %>% 
  filter(!is.na(TIN_10km)) %>% 
  mutate(year = year - data_summary$year,
         spec_PercentUrban = spec_PercentUrban - data_summary$spec_PercentUrban,
         spec_PercentAg = spec_PercentAg - data_summary$spec_PercentAg,
         TIN_10km = TIN_10km - data_summary$TIN_10km,
         tmean_10km = tmean_10km - data_summary$tmean_10km,
         ppt_10km = ppt_10km - data_summary$ppt_10km) %>%  
  mutate(Spp_index = as.numeric(as.factor(Spp_code))) 


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
  boundary = non_convex_bdry, max.edge = c(max.edge*4, max.edge*16), # km inside and outside
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

s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(TIN_10km + spec_PercentAg + spec_PercentUrban + ppt_10km)^3, model = "fixed")
  # scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  # collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  # space_int(coords, model = spde)



s_components.year <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(TIN_10km + spec_PercentAg + spec_PercentUrban + year + ppt_10km)^2, model = "fixed")+
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
##########  Plotting the prediction without year effects ###############
################################################################################################################################

min_ag<- min(data$spec_PercentAg)
mean_ag <- mean(data$spec_PercentAg)
max_ag <- max(data$spec_PercentAg)

min_urb<- min(data$spec_PercentUrban)
mean_urb <- mean(data$spec_PercentUrban)
max_urb<- max(data$spec_PercentUrban)

min_nit<- min(data$TIN_10km)
mean_nit <- mean(data$TIN_10km)
max_nit<- max(data$TIN_10km)

preddata.1 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     spec_PercentAg = rep(seq(min_ag, max_ag, length.out = 50), times = 3),
                     spec_PercentUrban = 0,
                     TIN_10km = 0,
                     ppt_10km = 0,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))
preddata.2 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     spec_PercentAg = 0,
                     spec_PercentUrban = rep(seq(min_urb, max_urb, length.out = 50), times = 3),
                     TIN_10km = 0,
                     ppt_10km = 0,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

preddata.3 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     spec_PercentAg = 0,
                     spec_PercentUrban = 0,
                     TIN_10km = rep(seq(min_nit, max_nit, length.out = 50), times = 3),
                     ppt_10km = 0,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

ag.pred <- predict(
  fit,
  newdata = preddata.1,
  formula = ~ invlogit(fixed),# + year_eval(year_index)),#+ collector_eval(collector_index) + scorer_eval(scorer_index) + year_eval(year_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 500) %>% 
  mutate(spec_PercentAg = spec_PercentAg + data_summary$spec_PercentAg) # undoing mean centering


urb.pred <- predict(
  fit,
  newdata = preddata.2,
  formula = ~ invlogit(fixed),# + year_eval(year_index)),#+ collector_eval(collector_index) + scorer_eval(scorer_index) + year_eval(year_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 500) %>% 
  mutate(spec_PercentUrban = spec_PercentUrban + data_summary$spec_PercentUrban) # undoing mean centering



nit.pred <- predict(
  fit,
  newdata = preddata.3,
  formula = ~ invlogit(fixed),# + year_eval(year_index)),#+ collector_eval(collector_index) + scorer_eval(scorer_index) + year_eval(year_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 500) %>% 
  mutate(TIN_10km = TIN_10km + data_summary$TIN_10km) # undoing mean centering






ag_binned <- data %>% 
  mutate(spec_PercentAg = spec_PercentAg + data_summary$spec_PercentAg) %>% # undoing mean centering
  mutate(ag_bin = cut(spec_PercentAg, breaks = 5)) %>% 
  group_by(species, ag_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_ag = mean(spec_PercentAg),
            sample = n())

urb_binned <- data %>% 
  mutate(spec_PercentUrban = spec_PercentUrban + data_summary$spec_PercentUrban) %>% # undoing mean centering
  mutate(urb_bin = cut(spec_PercentUrban, breaks = 5)) %>% 
  group_by(species, urb_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_urb = mean(spec_PercentUrban),
            sample = n())

nit_binned <- data %>% 
  mutate(TIN_10km = TIN_10km + data_summary$TIN_10km) %>% # undoing mean centering
  mutate(nit_bin = cut(TIN_10km, breaks = 5)) %>% 
  group_by(species, nit_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_nit = mean(TIN_10km),
            sample = n())

ag_trend <- ggplot(ag.pred) +
  geom_line(aes(x = spec_PercentAg, q0.5)) +
  geom_ribbon(aes(spec_PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#B38600") +
  geom_ribbon(aes(spec_PercentAg, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = ag_binned, aes(x = mean_ag, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,260))+
  facet_wrap(~species,  ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Ag. (%)",  size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))#+
# lims(y = c(0,1), x = c(0, 100))

urb_trend <- ggplot(urb.pred) +
  geom_line(aes(spec_PercentUrban, q0.5)) +
  geom_ribbon(aes(spec_PercentUrban, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#021475") +
  geom_ribbon(aes(spec_PercentUrban, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = urb_binned, aes(x = mean_urb, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,260))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))


nit_trend <- ggplot(nit.pred) +
  geom_line(aes(TIN_10km, q0.5)) +
  geom_ribbon(aes(TIN_10km, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#BF00A0") +
  geom_ribbon(aes(TIN_10km, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = nit_binned, aes(x = mean_nit, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,260))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Nitrogen Deposition (kg N/km^2)", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1.1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))



ag_trend <- tag_facet(ag_trend)
urb_trend <- tag_facet(urb_trend, tag_pool =  letters[-(1:3)])
nit_trend <- tag_facet(nit_trend, tag_pool =  letters[-(1:6)])
fig1 <-   ag_trend + urb_trend + nit_trend + plot_layout(ncol = 3, guides = "collect") 
ggsave(fig1, file = "Plots/year_specific_model1.png", width = 10, height = 8)



################################################################################################################################
##########  Plotting the posteriors from the model without year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
param_names <- fit.y.rfx$summary.random$fixed$ID

n_draws <- 500

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors <- generate(
  fit.y.rfx,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors) <- param_names
colnames(posteriors) <- c( paste0("iter",1:n_draws))


posteriors_df <- as_tibble(t(posteriors), rownames = "iteration")



# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  mutate(NIT.AGHY = mean_TIN_10km,
         NIT.AGPE = mean_TIN_10km+`Spp_codeAGPE:mean_TIN_10km`,
         NIT.ELVI = mean_TIN_10km+`Spp_codeELVI:mean_TIN_10km`,
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
         param_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(case_when(param_label == "INT" ~ "Intercept",
                                    param_label == "AG" ~ "Agric. Cover",
                                    param_label == "URB" ~ "Urban Cover",
                                    param_label == "NIT" ~ "Nit. Dep."), levels = c("Intercept", "Agric. Cover", "Urban Cover", "Nit. Dep.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus")))
  )





posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x")+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"))

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_yr.rfx.png", width = 6, height = 6)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/500,
            prob_neg = 1-prob_pos)
write.csv(effects_summary, file = "Posterior_prob_results.csv")

################################################################################################################################
##########  Assessing model fit     ###############
################################################################################################################################

# predicting the training data

validation.pred <- predict(
  fit.y.rfx,
  newdata = data,
  formula = ~ invlogit(fixed + year + scorer + collector + space_int),
  n.samples = 500) 


rocobj <- pROC::roc(data$Endo_status_liberal, validation.pred$mean)

ROC_training_plot <- ggroc(rocobj) 
ggsave(ROC_training_plot, filename = "ROC_training_plot.png", width = 4, height = 4)

# AUC values
rocobj$auc
# 0.7487


# generating posterior samples of each parameter
post.pred <- generate(
  fit.y.rfx,
  newdata = data,
  formula = ~ invlogit(fixed + year + scorer + collector + space_int),
  n.samples = 250) 

posterior_samples <- bind_cols(data, post.pred) %>% 
  select(...109:...358) %>% st_drop_geometry %>% as.matrix()

# simulating datasets from the posterior samples
n_post_draws <- 250
y_sim <- matrix(NA,n_post_draws,length(data$Endo_status_liberal))

for(i in 1:n_post_draws){
  y_sim[i,] <- rbinom(n = length(endo_herb$Endo_status_liberal), size = 1, prob = posterior_samples[,i])
}

# saveRDS(y_sim, file = "y_sim.rds")
# y_sim <- readRDS("y_sim.rds")
y_sim_df <- t(y_sim)
colnames(y_sim_df) <- paste("iter", 1:n_post_draws)
y_sim_df <- as_tibble(y_sim_df, .name_repair = "minimal") %>% 
  pivot_longer(cols = everything())

overlay_plot <- ggplot(data)+
  geom_line(data = y_sim_df, aes(x = value, group = name), stat="density", color = "red", alpha = .3) +
  geom_line(aes(x = Endo_status_liberal), stat="density") + 
  labs(x = "Endophyte Status", y = "Density")+
  theme_minimal()
overlay_plot
ggsave(overlay_plot, filename = "overlay_plot.png", width = 4, height = 4)





################################################################################################################################
##########  Getting and plotting prediction from NitXYear ###############
################################################################################################################################
mean_year <- mean(data$year)
summary_data <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_nit = quantile(mean_TIN_10km, .05),
            max_nit = quantile(mean_TIN_10km, .95),
            mean_nit = mean(mean_TIN_10km))

preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             Spp_index = 1,
                             # year = sort(unique(data$year)),
                             year = summary_data[summary_data$Spp_code == "AGHY",]$min_year:summary_data[summary_data$Spp_code == "AGHY",]$max_year,
                             mean_TIN_10km = c(summary_data[summary_data$Spp_code == "AGHY",]$min_nit, summary_data[summary_data$Spp_code == "AGHY",]$max_nit),
                             PercentAg = mean_ag,
                             PercentUrban = mean_urb)
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             Spp_index = 2,
                             # year = sort(unique(data$year)),
                             year = summary_data[summary_data$Spp_code == "AGPE",]$min_year:summary_data[summary_data$Spp_code == "AGPE",]$max_year,
                             mean_TIN_10km = c(summary_data[summary_data$Spp_code == "AGPE",]$min_nit, summary_data[summary_data$Spp_code == "AGPE",]$max_nit),
                             PercentAg = mean_ag,
                             PercentUrban = mean_urb)
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             Spp_index = 3,
                             # year = sort(unique(data$year)),
                             year = summary_data[summary_data$Spp_code == "ELVI",]$min_year:summary_data[summary_data$Spp_code == "ELVI",]$max_year,
                             mean_TIN_10km = c(summary_data[summary_data$Spp_code == "ELVI",]$min_nit, summary_data[summary_data$Spp_code == "ELVI",]$max_nit),
                             PercentAg = mean_ag,
                             PercentUrban = mean_urb)

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         nit_label = case_when(mean_TIN_10km > mean_nit ~ "High Nitrogen",
                               mean_TIN_10km < mean_nit ~ "Low Nitrogen"))


# gennerating predictions and back-transforming the standardized year variable



year.nit.pred <- predict(
  fit.y.rfx,
  newdata = preddata,
  formula = ~ invlogit( year), #+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100)  


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_nit = cut(mean_TIN_10km, breaks = 2),
         binned_year = cut(year, breaks = 50)) %>%
  group_by(Spp_code, species,binned_nit, binned_year) %>%
  summarise(mean_nit = mean(mean_TIN_10km),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(nit_bin = case_when(mean_nit>=0~ "High Nitrogen",
                             mean_nit<0 ~ "Low Nitrogen"))




# histogram_data <- endo_herb%>%
#   mutate(disturbance_label = case_when(NO3_mean > mean(NO3_mean) ~ "Nitrogen",
#                                     NO3_mean < mean(NO3_mean) ~ "Undisturbed"))
# histogram_data2 <- endo_herb%>%
#   mutate(nit_label = case_when(NO3_mean > mean(NO3_mean) ~ "High Nitrogen",
#                                        NO3_mean < mean(NO3_mean) ~ "Low Nitrogen"))

nit_yr_trend <- ggplot(year.nit.pred)+
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = nit_label, fill = nit_label), alpha = 0.3) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = nit_label, fill = nit_label), alpha = 0.3) +
  geom_line(aes(year, mean, group = nit_label), color = "black") +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, fill = nit_bin, size = sample), shape = 21, alpha = .6)+
  facet_wrap(~species, nrow = 3, scales = "free", strip.position="right")+
  guides(color = "none",
         fill = guide_legend(title.position="top", title.hjust = 0.5, order = 1),
         size = guide_legend(title.position="top", title.hjust = 0.5, order = 2))+
  scale_color_manual(values = c("High Nitrogen" = "#BF00A0",
                                "Low Nitrogen"="darkgray"))+
  scale_fill_manual(values = c("High Nitrogen" = "#BF00A0",
                               "Low Nitrogen"="darkgray"))+
  labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size", fill = "Nitrogen Deposition")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_text( size = rel(1.1)), strip.text.y.right = element_text(face = "italic",angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"),
        legend.position = "bottom", legend.box = "vertical")+
  lims(y = c(0,1), x = c(1832, 2020))

nit_yr_trend
# ggsave(nit_yr_trend, filename = "Nitrogen_and_Year.png", width = 6, height = 8)




################################################################################################################################
##########  Getting and plotting prediction from AgXYear ###############
################################################################################################################################
mean_year <- mean(data$year)
summary_data <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_ag = quantile(PercentAg, .05),
            max_ag = quantile(PercentAg, .95),
            mean_ag = mean(PercentAg))

preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             year = seq(summary_data[summary_data$Spp_code == "AGHY",]$min_year, summary_data[summary_data$Spp_code == "AGHY",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = c(summary_data[summary_data$Spp_code == "AGHY",]$min_ag, summary_data[summary_data$Spp_code == "AGHY",]$max_ag),
                             PercentUrban = 0)
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             year = seq(summary_data[summary_data$Spp_code == "AGPE",]$min_year, summary_data[summary_data$Spp_code == "AGPE",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = c(summary_data[summary_data$Spp_code == "AGPE",]$min_ag, summary_data[summary_data$Spp_code == "AGPE",]$max_ag),
                             PercentUrban = 0)
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             year = seq(summary_data[summary_data$Spp_code == "ELVI",]$min_year, summary_data[summary_data$Spp_code == "ELVI",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = c(summary_data[summary_data$Spp_code == "ELVI",]$min_ag, summary_data[summary_data$Spp_code == "ELVI",]$max_ag),
                             PercentUrban = 0)

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         ag_label = case_when(PercentAg >=0 ~ "High Agr. Cover",
                              PercentAg < 0 ~ "Low Agr. Cover"))


# gennerating predictions and back-transforming the standardized year variable



year.ag.pred <- predict(
  fit.3,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + mean_year)


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_ag = cut(PercentAg, breaks = 2),
         binned_year = cut(year, breaks = 50)) %>% 
  group_by(Spp_code, species,binned_ag, binned_year) %>%
  summarise(mean_ag = mean(PercentAg),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(ag_bin = case_when(mean_ag>=0~ "High Agr. Cover",
                            mean_ag<0 ~ "Low Agr. Cover"))




# histogram_data <- endo_herb%>%
#   mutate(disturbance_label = case_when(NO3_mean > mean(NO3_mean) ~ "Nitrogen",
#                                     NO3_mean < mean(NO3_mean) ~ "Undisturbed"))
# histogram_data2 <- endo_herb%>%
#   mutate(nit_label = case_when(NO3_mean > mean(NO3_mean) ~ "High Nitrogen",
#                                        NO3_mean < mean(NO3_mean) ~ "Low Nitrogen"))

ag_yr_trend <- ggplot(year.ag.pred)+
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = ag_label, fill = ag_label), alpha = 0.3) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = ag_label, fill = ag_label), alpha = 0.3) +
  geom_line(aes(year, mean, group = ag_label), color = "black") +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, fill = ag_bin, size = sample), shape = 21, alpha = .6)+
  facet_wrap(~species, nrow = 3, scales = "free")+
  guides(color = "none",
         fill = guide_legend(title.position="top", title.hjust = 0.5, order = 1),
         size = guide_legend(title.position="top", title.hjust = 0.5, order = 2))+
  scale_color_manual(values = c("High Agr. Cover" = "#B38600",
                                "Low Agr. Cover"="darkgray"))+
  scale_fill_manual(values = c("High Agr. Cover" = "#B38600",
                               "Low Agr. Cover"="darkgray"))+
  labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size", fill = "Agr. Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"),
        legend.position = "bottom", legend.box = "vertical")+
  lims(y = c(0,1), x = c(1832, 2020))

# ag_yr_trend
# ggsave(ag_yr_trend, filename = "Agriculture_and_Year.png", width = 6, height = 8)





################################################################################################################################
##########  Getting and plotting prediction from UrbXYear ###############
################################################################################################################################
mean_year <- mean(data$year)
summary_data <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_urb = quantile(PercentUrban, .05),
            max_urb = quantile(PercentUrban, .95),
            mean_urb = mean(PercentUrban))

preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             year = seq(summary_data[summary_data$Spp_code == "AGHY",]$min_year, summary_data[summary_data$Spp_code == "AGHY",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = 0,
                             PercentUrban = c(summary_data[summary_data$Spp_code == "AGHY",]$min_urb, summary_data[summary_data$Spp_code == "AGHY",]$max_urb))
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             year = seq(summary_data[summary_data$Spp_code == "AGPE",]$min_year, summary_data[summary_data$Spp_code == "AGPE",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = 0,
                             PercentUrban = c(summary_data[summary_data$Spp_code == "AGPE",]$min_urb, summary_data[summary_data$Spp_code == "AGPE",]$max_urb))
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             year = seq(summary_data[summary_data$Spp_code == "ELVI",]$min_year, summary_data[summary_data$Spp_code == "ELVI",]$max_year, length.out = 20), 
                             std_TIN = 0,
                             PercentAg = 0,
                             PercentUrban = c(summary_data[summary_data$Spp_code == "ELVI",]$min_urb, summary_data[summary_data$Spp_code == "ELVI",]$max_urb))

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         urb_label = case_when(PercentUrban >=0 ~ "High Urb. Cover",
                               PercentUrban < 0 ~ "Low Urb. Cover"))


# gennerating predictions and back-transforming the standardized year variable



year.urb.pred <- predict(
  fit.3,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + mean_year)


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_ag = cut(PercentUrban, breaks = 2),
         binned_year = cut(year, breaks = 50)) %>% 
  group_by(Spp_code, species,binned_ag, binned_year) %>%
  summarise(mean_urb = mean(PercentUrban),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(urb_bin =  case_when(mean_urb >=0 ~ "High Urb. Cover",
                              mean_urb < 0 ~ "Low Urb. Cover"))




# histogram_data <- endo_herb%>%
#   mutate(disturbance_label = case_when(NO3_mean > mean(NO3_mean) ~ "Nitrogen",
#                                     NO3_mean < mean(NO3_mean) ~ "Undisturbed"))
# histogram_data2 <- endo_herb%>%
#   mutate(nit_label = case_when(NO3_mean > mean(NO3_mean) ~ "High Nitrogen",
#                                        NO3_mean < mean(NO3_mean) ~ "Low Nitrogen"))

urb_yr_trend <- ggplot(year.urb.pred)+
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975, group = urb_label, fill = urb_label), alpha = 0.3) +
  geom_ribbon(aes(year, ymin = q0.25, ymax = q0.75, group = urb_label, fill = urb_label), alpha = 0.3) +
  geom_line(aes(year, mean, group = urb_label), color = "black") +
  geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, fill = urb_bin, size = sample), shape = 21, alpha = .6)+
  facet_wrap(~species, nrow = 3, scales = "free")+
  guides(color = "none",
         fill = guide_legend(title.position="top", title.hjust = 0.5, order = 1),
         size = guide_legend(title.position="top", title.hjust = 0.5, order = 2))+
  scale_color_manual(values = c("High Urb. Cover" = "#021475",
                                "Low Urb. Cover"="darkgray"))+
  scale_fill_manual(values = c("High Urb. Cover" = "#021475",
                               "Low Urb. Cover"="darkgray"))+
  labs(y = "Endophyte Prevalence", x = "Year", size = "Sample Size", fill = "Urb. Cover (%)")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"),
        legend.position = "bottom", legend.box = "vertical")+
  lims(y = c(0,1), x = c(1832, 2020))

# urb_yr_trend
# ggsave(urb_yr_trend, filename = "Urban_and_Year.png", width = 6, height = 8)


##########
##### Plotting all three driver X year plots together #####

ag_yr_trend <- tag_facet(ag_yr_trend)
urb_yr_trend <- tag_facet(urb_yr_trend, tag_pool =  letters[-(1:3)])
nit_yr_trend <- tag_facet(nit_yr_trend, tag_pool =  letters[-(1:6)])
yr_trend_plot <- ag_yr_trend + urb_yr_trend + nit_yr_trend + plot_layout(ncol = 3) 


#yr_trend_plot

ggsave(yr_trend_plot, filename = "yr_trend_plot.png", width = 10, height = 8)





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
  mutate(NIT.AGHY = std_TIN,
         NIT.AGPE = std_TIN+`Spp_codeAGPE:std_TIN`,
         NIT.ELVI = std_TIN+`Spp_codeELVI:std_TIN`,
         URB.AGHY = std_urb,
         URB.AGPE = std_urb+`Spp_codeAGPE:std_urb`,
         URB.ELVI = std_urb+`Spp_codeELVI:std_urb`,
         AG.AGHY = std_ag,
         AG.AGPE = std_ag+`Spp_codeAGPE:std_urb`,
         AG.ELVI = std_ag+`Spp_codeELVI:std_urb`,
         YEAR.AGHY = std_year,
         YEAR.AGPE = std_year+`Spp_codeAGPE:std_year`,
         YEAR.ELVI = std_year+`Spp_codeELVI:std_year`,
         YEARxNIT.AGHY = `std_year:std_TIN`,
         YEARxNIT.AGPE = `std_year:std_TIN`+`Spp_codeAGPE:std_year:std_TIN`,
         YEARxNIT.ELVI = `std_year:std_TIN`+`Spp_codeELVI:std_year:std_TIN`,
         YEARxURB.AGHY = `std_year:std_urb`,
         YEARxURB.AGPE = `std_year:std_urb`+`Spp_codeAGPE:std_year:std_urb`,
         YEARxURB.ELVI = `std_year:std_urb`+`Spp_codeELVI:std_year:std_urb`,
         YEARxAG.AGHY = `std_year:std_ag`,
         YEARxAG.AGPE = `std_year:std_ag`+`Spp_codeAGPE:std_year:std_ag`,
         YEARxAG.ELVI = `std_year:std_ag`+`Spp_codeELVI:std_year:std_ag`,
         INT.AGHY = Spp_codeAGHY,
         INT.AGPE = Spp_codeAGPE,
         INT.ELVI = Spp_codeELVI) %>% 
  select(-all_of(param_names)) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "Temporal") %>% 
  mutate(spp_label = sub(".*\\.", "", param),
         param_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(case_when(param_label == "INT" ~ "Intercept",
                                    param_label == "AG" ~ "Agric. Cover",
                                    param_label == "URB" ~ "Urban Cover",
                                    param_label == "NIT" ~ "Nit. Dep.",
                                    param_label == "YEAR" ~ "Year", 
                                    param_label == "YEARxNIT" ~ "Year x Nit.", 
                                    param_label == "YEARxURB" ~ "Year x Urb.",
                                    param_label == "YEARxAG" ~ "Year x Agric."), levels = c("Intercept", "Agric. Cover", "Urban Cover", "Nit. Dep.", "Year",
                                                                                            "Year x Agric.", "Year x Urb.", "Year x Nit.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus")))
  )



posterior_hist <- ggplot(effects_df)+
  stat_halfeye(aes(x = value, y = spp_f, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", nrow = 2)+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"))

posterior_hist
ggsave(posterior_hist, filename = "posterior_hist_with_year.png", width = 8, height = 7)



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




