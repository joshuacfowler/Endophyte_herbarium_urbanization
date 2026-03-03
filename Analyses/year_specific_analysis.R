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
                   PercentAg = mean(PercentAg, na.rm = T),
                   PercentUrban = mean(PercentUrban, na.rm = T),
                   mean_TIN_10km = mean(mean_TIN_10km, na.rm = T),
                   spec_PercentAg = mean(spec_PercentAg, na.rm = T),
                   spec_PercentUrban = mean(spec_PercentUrban, na.rm = T),
                   TIN_10km = mean(TIN_10km, na.rm = T),
                   tmean_10km = mean(tmean_10km, na.rm = T),
                   ppt_10km = mean(ppt_10km, na.rm = T))

data <- contemp_endo_herb %>% 
  filter(!is.na(ppt_10km)) %>% 
  filter(!is.na(TIN_10km)) %>% 
  mutate(year = year - data_summary$year,
         PercentUrban = PercentUrban - data_summary$PercentUrban,
         PercentAg = PercentAg - data_summary$PercentAg,
         mean_TIN_10km = mean_TIN_10km - data_summary$mean_TIN_10km,
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

s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(TIN_10km + spec_PercentAg + spec_PercentUrban + ppt_10km)^4, model = "fixed")
  # scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  # collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  # space_int(coords, model = spde)


s_components_normal <-  ~ 0 +  fixed(main = ~ 0 + Spp_code/(mean_TIN_10km + PercentAg + PercentUrban + ppt_10km)^4, model = "fixed")

s_components.year <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(TIN_10km + spec_PercentAg + spec_PercentUrban + year + ppt_10km)^2, model = "fixed")
  # scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  # collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  # space_int(coords, model = spde)
s_components.year_normal <-  ~ 0 +  fixed(main = ~ 0 + (Spp_code)/(mean_TIN_10km + PercentAg + PercentUrban + year + ppt_10km)^2, model = "fixed")


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

fit_normal <- bru(s_components_normal,
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


fit.year_normal <- bru(s_components.year_normal,
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
  mutate(ag_bin = cut(spec_PercentAg, breaks = 10)) %>% 
  group_by(species, ag_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_ag = mean(spec_PercentAg),
            sample = n())

urb_binned <- data %>% 
  mutate(spec_PercentUrban = spec_PercentUrban + data_summary$spec_PercentUrban) %>% # undoing mean centering
  mutate(urb_bin = cut(spec_PercentUrban, breaks = 10)) %>% 
  group_by(species, urb_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_urb = mean(spec_PercentUrban),
            sample = n())

nit_binned <- data %>% 
  mutate(TIN_10km = TIN_10km + data_summary$TIN_10km) %>% # undoing mean centering
  mutate(nit_bin = cut(TIN_10km, breaks = 10)) %>% 
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
# generating posterior for the year-specific model
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
  mutate(model = "No Year",
         data = "Year-specific") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "spec_PercentAg" = "Agr.",
                                                         "spec_PercentUrban" = "Urb.",
                                                         "TIN_10km" = "Nit.",
                                                         "ppt_10km" = "PPT.")),
                          levels = c("Intercept","Nit.","Agr.","Urb.","PPT.","Nit. X Agr.","Nit. X Urb.","Nit. X PPT.","Agr. X Urb.","Agr. X PPT.","Urb. X PPT.","Nit. X Agr. X Urb.","Nit. X Agr. X PPT.","Nit. X Urb. X PPT.","Agr. X Urb. X PPT.","Nit. X Agr. X Urb. X PPT.")),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = (c("A. hyemalis", "A. perennans", "E. virginicus"))))


# generating the same for the normal data model

param_names_normal <- fit_normal$summary.random$fixed$ID

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors_normal <- generate(
  fit_normal,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors_normal) <- param_names_normal
colnames(posteriors_normal) <- c( paste0("iter",1:n_draws))


posteriors_normal_df <- as_tibble(t(posteriors_normal), rownames = "iteration")

colnames(posteriors_normal_df) <- sub("Spp_code", "", colnames(posteriors_normal_df))
colnames(posteriors_normal_df) <- gsub(":", ".", colnames(posteriors_normal_df))

# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_normal_df <- posteriors_normal_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year",
         data = "Normal") %>% 
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
                        levels = (c("A. hyemalis", "A. perennans", "E. virginicus"))))





effects <- bind_rows(effects_df, effects_normal_df) 

effects_summary <- effects %>%
  group_by(param, param_label, spp_label, spp_f, param_f, data) %>%
  summarize(mean = mean(value),
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)


data_colors <- c("#2c7fb8", "#808080")


posterior_hist <- ggplot(effects_summary)+
  coord_flip()+
  geom_hline(aes(yintercept = 0), lwd = 1)+
  geom_point(aes(y = mean, x = param_f, group = data, color = data), position = position_dodge(width = .4), size =2)+
  geom_linerange(aes(ymin = lwr, ymax = upr, x = param_f, group = data, color = data),  position = position_dodge(width = 0.4))+
  ggh4x::facet_grid2( param_f ~ spp_f, remove_labels = "y", scales = "free", independent = "x")+
  labs(x = "Parameter", y = "Posterior Est.")+
  # guides(fill = data)+
  scale_color_manual(values = data_colors)+
  scale_fill_manual(values = data_colors)+
  scale_y_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.x = element_text(size = rel(.8)),
                     strip.background = element_blank(),
                     strip.text.y = element_blank(),
                     strip.text.x = element_text(face = "italic"))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/year_specific_posterior_hist.png", width = 9, height = 9)




effects_AGPE <-bind_rows(effects_df, effects_normal_df) %>% 
  filter(spp_label == "AGPE")

# plotting the same for AGPE
posterior_hist <- ggplot(effects_AGPE)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = data, group = data), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free", ncol = 4)+
  labs(x = "Posterior Est.", y = "Species")+
  # guides(fill = data)+
  scale_color_manual(values = data_colors)+
  scale_fill_manual(values = data_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

ggsave(posterior_hist, filename = "Plots/year_specific_posterior_hist_AGPE.png", width = 11, height = 6)


# 






################################################################################################################################
##########  Plotting the posteriors from the model with year effect ###############
################################################################################################################################

param_names <- fit.year$summary.random$fixed$ID

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
  mutate(model = "Year",
         data = "Year-specific") %>% 
  mutate(param_label = sub("^[^.]+.", "", param),
         spp_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(str_replace_all(param_label, c("\\." = " X ",
                                                         "Int" = "Intercept",
                                                         "year" = "Year",
                                                         "spec_PercentAg" = "Agr.",
                                                         "spec_PercentUrban" = "Urb.",
                                                         "TIN_10km" = "Nit.",
                                                         "ppt_10km" = "Ppt.")),
                          levels = c("Intercept"  ,"Nit.","Agr.","Urb.","Year","Ppt.","Nit. X Agr.","Nit. X Urb.", "Nit. X Year","Nit. X Ppt.","Agr. X Urb.","Agr. X Year","Agr. X Ppt.","Urb. X Year","Urb. X Ppt.","Year X Ppt.","Nit. X Agr. X Urb.","Nit. X Agr. X Year","Nit. X Agr. X Ppt.","Nit. X Urb. X Year","Nit. X Urb. X Ppt.","Nit. X Year X Ppt.","Agr. X Urb. X Year","Agr. X Urb. X Ppt.","Agr. X Year X Ppt.","Urb. X Year X Ppt.","Nit. X Agr. X Urb. X Year","Nit. X Agr. X Urb. X Ppt.","Nit. X Agr. X Year X Ppt." ,"Nit. X Urb. X Year X Ppt.","Agr. X Urb. X Year X Ppt." )),
         spp_f = factor(case_when(spp_label == "ELVI" ~ "E. virginicus", spp_label == "AGPE" ~ "A. perennans", spp_label == "AGHY" ~ "A. hyemalis"),
                        levels = rev(c("A. hyemalis", "A. perennans", "E. virginicus"))))




# generating the same for the normal data model

param_names_normal <- fit.year_normal$summary.random$fixed$ID

# we can sample values from the join posteriors of the parameters with the addition of "_latent" to the parameter name
posteriors_normal <- generate(
  fit.year_normal,
  formula = ~ fixed_latent,
  n.samples = n_draws) 
rownames(posteriors_normal) <- param_names_normal
colnames(posteriors_normal) <- c( paste0("iter",1:n_draws))


posteriors_normal_df <- as_tibble(t(posteriors_normal), rownames = "iteration")

colnames(posteriors_normal_df) <- sub("Spp_code", "", colnames(posteriors_normal_df))
colnames(posteriors_normal_df) <- gsub(":", ".", colnames(posteriors_normal_df))

# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_normal_df <- posteriors_normal_df %>% 
  rename(AGHY.Int = AGHY, AGPE.Int = AGPE, ELVI.Int = ELVI) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "Year",
         data = "Normal") %>% 
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






effects <- bind_rows(effects_df, effects_normal_df) %>% 
  filter(spp_label != "AGPE")

data_colors <- c("#2c7fb8", "#808080")


posterior_hist <- ggplot(effects)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = data, group = data), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free_x", ncol = 4)+
  labs(x = "Posterior Est.", y = "Species")+
  # guides(fill = "none")+
  scale_color_manual(values = data_colors)+
  scale_fill_manual(values = data_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/posterior_hist_yearmodel_year_specific.png", width = 9, height = 11)


effects_AGPE <-bind_rows(effects_df, effects_normal_df) %>% 
  filter(spp_label == "AGPE")

# plotting the same for AGPE
posterior_hist <- ggplot(effects_AGPE)+
  stat_halfeye(aes(x = value, y = spp_f, fill  = data, group = data), breaks = 50, normalize = "panels", alpha = .6)+
  
  # stat_halfeye(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, normalize = "panels", alpha = .6)+
  # stat_histinterval(aes(x = value, y = spp_label, fill = spp_label), breaks = 50, alpha = .6)+
  # geom_point(data = posteriors_summary, aes(x = mean, y = spp_label, color = spp_label))+
  # geom_linerange(data = posteriors_summary, aes(xmin = lwr, xmax = upr, y = spp_label, color = spp_label))+
  
  geom_vline(xintercept = 0)+
  facet_wrap(~param_f, scales = "free", ncol = 4)+
  labs(x = "Posterior Est.", y = "Species")+
  # guides(fill = data)+
  scale_color_manual(values = data_colors)+
  scale_fill_manual(values = data_colors)+
  scale_x_continuous(labels = scales::label_scientific(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

ggsave(posterior_hist, filename = "Plots/posterior_hist_yearmodel_year_specific_AGPE.png", width = 11, height = 6)


effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)

