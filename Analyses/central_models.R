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
         ppt_10km = ppt_10km - data_summary$ppt_10km)

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
ggsave(mesh_plot, filename = "Plots/mesh_plot.png", width = 6, height = 5)



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



s_components <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*mean_TIN_10km*PercentAg*PercentUrban, model = "fixed")+
  scorer(scorer_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$scorer_index)), hyper = list(pc_prec)) +
  collector(collector_index, model = "iid", constr = TRUE, mapper = bru_mapper_index(max(data$collector_index, na.rm = T)), hyper = list(pc_prec))+
  space_int(coords, model = spde)



s_components.year <-  ~ 0 +  fixed(main = ~ 0 + Spp_code*year*mean_TIN_10km*PercentAg*PercentUrban, model = "fixed")+
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

min_ag<- min(data$PercentAg)
max_ag <- max(data$PercentAg)

min_urb<- min(data$PercentUrban)
max_urb<- max(data$PercentUrban)

min_nit<- min(data$mean_TIN_10km)
max_nit<- max(data$mean_TIN_10km)

preddata.1 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = rep(seq(min_ag, max_ag, length.out = 50), times = 3),
                     PercentUrban = 0,
                     mean_TIN_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))
preddata.2 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = 0,
                     PercentUrban = rep(seq(min_urb, max_urb, length.out = 50), times = 3),
                     mean_TIN_10km = 0,
                     year_index = 9999,
                     collector_index = 9999, scorer_index = 9999) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]))

preddata.3 <- tibble(Spp_code = c(rep("AGHY", times = 50),rep("AGPE",times = 50),rep("ELVI",times = 50)),
                     PercentAg = 0,
                     PercentUrban = 0,
                     mean_TIN_10km = rep(seq(min_nit, max_nit, length.out = 50), times = 3),
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
  mutate(PercentAg = PercentAg + data_summary$PercentAg) # undoing mean centering

urb.pred <- predict(
  fit,
  newdata = preddata.2,
  formula = ~ invlogit(fixed),#  + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(PercentUrban = PercentUrban + data_summary$PercentUrban) # undoing mean centering



nit.pred <- predict(
  fit,
  newdata = preddata.3,
  formula = ~ invlogit(fixed),# + collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(mean_TIN_10km = mean_TIN_10km + data_summary$mean_TIN_10km) # undoing mean centering



values <-  c("#b2abd2", "#5e3c99")


ag_binned <- endo_herb %>% 
  mutate(ag_bin = cut(PercentAg, breaks = 30)) %>% 
  group_by(species, ag_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_ag = mean(PercentAg),
            sample = n())

urb_binned <- endo_herb %>% 
  mutate(urb_bin = cut(PercentUrban, breaks = 30)) %>% 
  group_by(species, urb_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_urb = mean(PercentUrban),
            sample = n())

nit_binned <- endo_herb %>% 
  mutate(nit_bin = cut(mean_TIN_10km, breaks = 30)) %>% 
  group_by(species, nit_bin) %>% 
  summarize(mean_endo = mean(Endo_status_liberal),
            mean_nit = mean(mean_TIN_10km),
            sample = n())

ag_trend <- ggplot(ag.pred) +
  geom_line(aes(x = PercentAg, mean)) +
  geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#B38600") +
  geom_ribbon(aes(PercentAg, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = ag_binned, aes(x = mean_ag, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,310))+
  facet_wrap(~species,  ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Percent Ag. (%)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
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
  labs(y = "Endophyte Prevalence", x = "Percent Urban (%)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0,.1,.1,.1), "line"))


nit_trend <- ggplot(nit.pred) +
  geom_line(aes(mean_TIN_10km, mean)) +
  geom_ribbon(aes(mean_TIN_10km, ymin = q0.025, ymax = q0.975), alpha = 0.2, fill = "#BF00A0") +
  geom_ribbon(aes(mean_TIN_10km, ymin = q0.25, ymax = q0.75), alpha = 0.2) +
  geom_point(data = nit_binned, aes(x = mean_nit, y = mean_endo, size = sample), color = "black", shape = 21)+
  scale_size_continuous(limits=c(1,310))+
  facet_wrap(~species, ncol = 1, scales = "free_x", strip.position="right")+  
  labs(y = "Endophyte Prevalence", x = "Nitrogen Deposition (kg N/km^2)", color = "Year", fill = "Year", shape = "Year", size = "Sample Size")+
  theme_classic()+
  theme(strip.background = element_blank(), 
        strip.text = element_text( size = rel(1.1)), strip.text.y.right = element_text(face = "italic", angle = 0),
        plot.margin = unit(c(0,.1,.1,.1), "line"))



ag_trend <- tag_facet(ag_trend)
urb_trend <- tag_facet(urb_trend, tag_pool =  letters[-(1:3)])
nit_trend <- tag_facet(nit_trend, tag_pool =  letters[-(1:6)])
fig2 <-   ag_trend + urb_trend + nit_trend + plot_layout(ncol = 3, guides = "collect") 
ggsave(fig2, file = "Plots/Figure_2.png", width = 10, height = 8)




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



# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  mutate(INT.AGHY = Spp_codeAGHY,
         INT.AGPE = Spp_codeAGPE,
         INT.ELVI = Spp_codeELVI,
         NIT.AGHY = mean_TIN_10km,
         NIT.AGPE = mean_TIN_10km+`Spp_codeAGPE:mean_TIN_10km`,
         NIT.ELVI = mean_TIN_10km+`Spp_codeELVI:mean_TIN_10km`,
         AG.AGHY = PercentAg,
         AG.AGPE = PercentAg+`Spp_codeAGPE:PercentAg`,
         AG.ELVI = PercentAg+`Spp_codeELVI:PercentAg`,
         URB.AGHY = PercentUrban,
         URB.AGPE = PercentUrban+`Spp_codeAGPE:PercentUrban`,
         URB.ELVI = PercentUrban+`Spp_codeELVI:PercentUrban`,
         NITxAG.AGHY = `mean_TIN_10km:PercentAg`,
         NITxAG.AGPE = `mean_TIN_10km:PercentAg` + `Spp_codeAGPE:mean_TIN_10km:PercentAg`,
         NITxAG.ELVI = `mean_TIN_10km:PercentAg` + `Spp_codeELVI:mean_TIN_10km:PercentAg`,
         NITxURB.AGHY= `mean_TIN_10km:PercentUrban`,
         NITxURB.AGPE = `mean_TIN_10km:PercentUrban` + `Spp_codeAGPE:mean_TIN_10km:PercentUrban`,
         NITxURB.ELVI = `mean_TIN_10km:PercentUrban` + `Spp_codeELVI:mean_TIN_10km:PercentUrban`,
         AGxURB.AGHY= `PercentAg:PercentUrban`,
         AGxURB.AGPE = `PercentAg:PercentUrban` + `Spp_codeAGPE:PercentAg:PercentUrban`,
         AGxURB.ELVI = `PercentAg:PercentUrban` + `Spp_codeELVI:PercentAg:PercentUrban`,
         NITxAGxURB.AGHY= `mean_TIN_10km:PercentAg:PercentUrban`,
         NITxAGxURB.AGPE= `mean_TIN_10km:PercentAg:PercentUrban` + `Spp_codeAGPE:mean_TIN_10km:PercentAg:PercentUrban`,
         NITxAGxURB.ELVI= `mean_TIN_10km:PercentAg:PercentUrban` + `Spp_codeELVI:mean_TIN_10km:PercentAg:PercentUrban`) %>% 
  select(-all_of(param_names)) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year") %>% 
  mutate(spp_label = sub(".*\\.", "", param),
         param_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(case_when(param_label == "INT" ~ "Intercept",
                                    param_label == "AG" ~ "Agric. Cover",
                                    param_label == "URB" ~ "Urban Cover",
                                    param_label == "NIT" ~ "Nit. Dep.",
                                    param_label == "NITxAG" ~ "Nit. X Agr.",
                                    param_label == "NITxURB" ~ "Nit. X Urb.",
                                    param_label == "AGxURB" ~ "Agr. X Urb.",
                                    param_label == "NITxAGxURB" ~ "Nit. X Agr. X Urb."), levels = c("Intercept", "Agric. Cover", "Urban Cover", "Nit. Dep.", "Nit. X Agr.", "Nit. X Urb.", "Agr. X Urb.", "Nit. X Agr. X Urb.")),
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
  facet_wrap(~param_f, scales = "free_x", nrow = 2)+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/posterior_hist.png", width = 9, height = 5)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)
# write.csv(effects_summary, file = "Posterior_prob_results.csv")

################################################################################################################################
##########  Assessing model fit     ###############
################################################################################################################################

# predicting the training data

validation.pred <- predict(
  fit,
  newdata = data,
  formula = ~ invlogit(fixed + scorer + collector + space_int),
  n.samples = 500) 


rocobj <- pROC::roc(data$Endo_status_liberal, validation.pred$mean)

ROC_training_plot <- ggroc(rocobj) 
ggsave(ROC_training_plot, filename = "Plots/ROC_training_plot.png", width = 4, height = 4)

# AUC values
rocobj$auc
# 0.7555


# generating posterior samples of each parameter
post.pred <- generate(
  fit,
  newdata = data,
  formula = ~ invlogit(fixed + scorer + collector + space_int),
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
ggsave(overlay_plot, filename = "Plots/overlay_plot.png", width = 4, height = 4)








################################################################################################################################
##########  Getting and plotting prediction from NitXYear ###############
################################################################################################################################

quantile_summary <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_nit = quantile(mean_TIN_10km, .05),
            max_nit = quantile(mean_TIN_10km, .95))


preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             # year = sort(unique(data$year)),
                             year = quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_year:quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_year,
                             mean_TIN_10km = c(quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_nit,quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_nit),
                             PercentAg = 0,
                             PercentUrban = 0)
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             # year = sort(unique(data$year)),
                             year = quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_year:quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_year,
                             mean_TIN_10km = c(quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_nit,quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_nit),
                             PercentAg = 0,
                             PercentUrban = 0)
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             # year = sort(unique(data$year)),
                             year = quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_year:quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_year,
                             mean_TIN_10km = c(quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_nit,quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_nit),
                             PercentAg = 0,
                             PercentUrban = 0)

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         nit_label = case_when(mean_TIN_10km > 0 ~ "High Nitrogen",
                               mean_TIN_10km < 0 ~ "Low Nitrogen"))


# gennerating predictions and back-transforming the standardized year variable



year.nit.pred <- predict(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed), #+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100)  %>% 
  mutate(year = year + data_summary$year,
         mean_TIN_10km = mean_TIN_10km + data_summary$mean_TIN_10km)


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_nit = cut(mean_TIN_10km, breaks = 2),
         binned_year = cut(year, breaks = 20)) %>%
  group_by(Spp_code, species,binned_nit, binned_year) %>%
  summarise(mean_nit = mean(mean_TIN_10km),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(nit_bin = case_when(mean_nit>=data_summary$mean_TIN_10km ~ "High Nitrogen",
                             mean_nit<data_summary$mean_TIN_10km ~ "Low Nitrogen"))




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

# nit_yr_trend
# ggsave(nit_yr_trend, filename = "Nitrogen_and_Year.png", width = 6, height = 8)




################################################################################################################################
##########  Getting and plotting prediction from AgXYear ###############
################################################################################################################################
quantile_summary <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_ag = quantile(PercentAg, .05),
            max_ag = quantile(PercentAg, .95))

preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_year, quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = c(quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_ag, quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_ag),
                             PercentUrban = 0)
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_year, quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = c(quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_ag, quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_ag),
                             PercentUrban = 0)
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_year, quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = c(quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_ag, quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_ag),
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
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + data_summary$year,
         PercentAg = PercentAg + data_summary$PercentAg)


# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_ag = cut(PercentAg, breaks = 2),
         binned_year = cut(year, breaks = 20)) %>% 
  group_by(Spp_code, species,binned_ag, binned_year) %>%
  summarise(mean_ag = mean(PercentAg),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(ag_bin = case_when(mean_ag>=data_summary$PercentAg ~ "High Agr. Cover",
                            mean_ag<data_summary$PercentAg ~ "Low Agr. Cover"))




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
quantile_summary <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_urb = quantile(PercentUrban, .05),
            max_urb = quantile(PercentUrban, .95))


preddata_aghy <- expand.grid(Spp_code = c("AGHY"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_year, quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = 0,
                             PercentUrban = c(quantile_summary[quantile_summary$Spp_code == "AGHY",]$min_urb, quantile_summary[quantile_summary$Spp_code == "AGHY",]$max_urb))
preddata_agpe <- expand.grid(Spp_code = c("AGPE"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_year, quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = 0,
                             PercentUrban = c(quantile_summary[quantile_summary$Spp_code == "AGPE",]$min_urb, quantile_summary[quantile_summary$Spp_code == "AGPE",]$max_urb))
preddata_elvi <- expand.grid(Spp_code = c("ELVI"), 
                             year = seq(quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_year, quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_year, length.out = 20), 
                             mean_TIN_10km = 0,
                             PercentAg = 0,
                             PercentUrban = c(quantile_summary[quantile_summary$Spp_code == "ELVI",]$min_urb, quantile_summary[quantile_summary$Spp_code == "ELVI",]$max_urb))

preddata <- bind_rows(preddata_aghy, preddata_agpe, preddata_elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         urb_label = case_when(PercentUrban >=0 ~ "High Urb. Cover",
                               PercentUrban < 0 ~ "Low Urb. Cover"))


# gennerating predictions and back-transforming the standardized year variable



year.urb.pred <- predict(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + data_summary$year,
         PercentUrban = PercentUrban + data_summary$PercentUrban)



# binning the data for plotting
endo_herb_binned <- endo_herb %>%
  mutate(binned_ag = cut(PercentUrban, breaks = 2),
         binned_year = cut(year, breaks = 20)) %>% 
  group_by(Spp_code, species,binned_ag, binned_year) %>%
  summarise(mean_urb = mean(PercentUrban),
            mean_year = mean(year),
            mean_endo = mean(Endo_status_liberal),
            sample = n()) %>%
  mutate(urb_bin =  case_when(mean_urb >=data_summary$PercentUrban ~ "High Urb. Cover",
                              mean_urb < data_summary$PercentUrban ~ "Low Urb. Cover"))




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

ggsave(yr_trend_plot, filename = "Plots/yr_trend_plot.png", width = 10, height = 8)



################################################################################################################################
##########  Plotting the posteriors from the model without year effect ###############
################################################################################################################################

# param_names <- fit.4$summary.random$fixed$ID
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



# Calculate the effects of the predictor, given that the reference level is for AGHY
effects_df <- posteriors_df %>% 
  mutate(INT.AGHY = Spp_codeAGHY,
         INT.AGPE = Spp_codeAGPE,
         INT.ELVI = Spp_codeELVI,
         NIT.AGHY = mean_TIN_10km,
         NIT.AGPE = mean_TIN_10km+`Spp_codeAGPE:mean_TIN_10km`,
         NIT.ELVI = mean_TIN_10km+`Spp_codeELVI:mean_TIN_10km`,
         AG.AGHY = PercentAg,
         AG.AGPE = PercentAg+`Spp_codeAGPE:PercentAg`,
         AG.ELVI = PercentAg+`Spp_codeELVI:PercentAg`,
         URB.AGHY = PercentUrban,
         URB.AGPE = PercentUrban+`Spp_codeAGPE:PercentUrban`,
         URB.ELVI = PercentUrban+`Spp_codeELVI:PercentUrban`,
         NITxAG.AGHY = `mean_TIN_10km:PercentAg`,
         NITxAG.AGPE = `mean_TIN_10km:PercentAg` + `Spp_codeAGPE:mean_TIN_10km:PercentAg`,
         NITxAG.ELVI = `mean_TIN_10km:PercentAg` + `Spp_codeELVI:mean_TIN_10km:PercentAg`,
         NITxURB.AGHY= `mean_TIN_10km:PercentUrban`,
         NITxURB.AGPE = `mean_TIN_10km:PercentUrban` + `Spp_codeAGPE:mean_TIN_10km:PercentUrban`,
         NITxURB.ELVI = `mean_TIN_10km:PercentUrban` + `Spp_codeELVI:mean_TIN_10km:PercentUrban`,
         AGxURB.AGHY= `PercentAg:PercentUrban`,
         AGxURB.AGPE = `PercentAg:PercentUrban` + `Spp_codeAGPE:PercentAg:PercentUrban`,
         AGxURB.ELVI = `PercentAg:PercentUrban` + `Spp_codeELVI:PercentAg:PercentUrban`,
         NITxAGxURB.AGHY= `mean_TIN_10km:PercentAg:PercentUrban`,
         NITxAGxURB.AGPE= `mean_TIN_10km:PercentAg:PercentUrban` + `Spp_codeAGPE:mean_TIN_10km:PercentAg:PercentUrban`,
         NITxAGxURB.ELVI= `mean_TIN_10km:PercentAg:PercentUrban` + `Spp_codeELVI:mean_TIN_10km:PercentAg:PercentUrban`,
         YEAR.AGHY = `year`,
         YEAR.AGPE = `year`+`Spp_codeAGPE:year`,
         YEAR.ELVI = `year`+`Spp_codeELVI:year`, 
         YEARxAG.AGHY = `year:PercentAg`,
         YEARxAG.AGPE = `year:PercentAg`+`Spp_codeAGPE:year:PercentAg`,
         YEARxAG.ELVI = `year:PercentAg`+`Spp_codeELVI:year:PercentAg`,  
         YEARxURB.AGHY = `year:PercentUrban`,
         YEARxURB.AGPE = `year:PercentUrban`+`Spp_codeAGPE:year:PercentUrban`,
         YEARxURB.ELVI = `year:PercentUrban`+`Spp_codeELVI:year:PercentUrban`, 
         YEARxNIT.AGHY = `year:mean_TIN_10km`,
         YEARxNIT.AGPE = `year:mean_TIN_10km`+`Spp_codeAGPE:year:mean_TIN_10km`,
         YEARxNIT.ELVI = `year:mean_TIN_10km`+`Spp_codeELVI:year:mean_TIN_10km`, 
         YEARxNITxAG.AGHY = `year:mean_TIN_10km:PercentAg`,
         YEARxNITxAG.AGPE = `year:mean_TIN_10km:PercentAg`+`Spp_codeAGPE:year:mean_TIN_10km:PercentAg`,
         YEARxNITxAG.ELVI = `year:mean_TIN_10km:PercentAg`+`Spp_codeELVI:year:mean_TIN_10km:PercentAg`, 
         YEARxNITxURB.AGHY = `year:mean_TIN_10km:PercentUrban`,
         YEARxNITxURB.AGPE = `year:mean_TIN_10km:PercentUrban`+`Spp_codeAGPE:year:mean_TIN_10km:PercentUrban`,
         YEARxNITxURB.ELVI = `year:mean_TIN_10km:PercentUrban`+`Spp_codeELVI:year:mean_TIN_10km:PercentUrban`, 
         YEARxAGxURB.AGHY = `year:PercentAg:PercentUrban`,
         YEARxAGxURB.AGPE = `year:PercentAg:PercentUrban`+`Spp_codeAGPE:year:PercentAg:PercentUrban`,
         YEARxAGxURB.ELVI = `year:PercentAg:PercentUrban`+`Spp_codeELVI:year:PercentAg:PercentUrban`, 
         YEARxNITxAGxURB.AGHY = `year:mean_TIN_10km:PercentAg:PercentUrban`,
         YEARxNITxAGxURB.AGPE = `year:mean_TIN_10km:PercentAg:PercentUrban`+`Spp_codeAGPE:year:mean_TIN_10km:PercentAg:PercentUrban`,
         YEARxNITxAGxURB.ELVI = `year:mean_TIN_10km:PercentAg:PercentUrban`+`Spp_codeELVI:year:mean_TIN_10km:PercentAg:PercentUrban`) %>% 
  select(-all_of(param_names)) %>% 
  pivot_longer( cols = -c(iteration), names_to = "param") %>% 
  mutate(model = "No Year") %>% 
  mutate(spp_label = sub(".*\\.", "", param),
         param_label = sub("\\..*","", param)) %>% 
  mutate(param_f = factor(case_when(param_label == "INT" ~ "Intercept",
                                    param_label == "AG" ~ "Agric. Cover",
                                    param_label == "URB" ~ "Urban Cover",
                                    param_label == "NIT" ~ "Nit. Dep.",
                                    param_label == "NITxAG" ~ "Nit. X Agr.",
                                    param_label == "NITxURB" ~ "Nit. X Urb.",
                                    param_label == "AGxURB" ~ "Agr. X Urb.",
                                    param_label == "NITxAGxURB" ~ "Nit. X Agr. X Urb.",
                                    param_label == "YEAR" ~ "Year",
                                    param_label == "YEARxNIT" ~ "Year X Nit.",
                                    param_label == "YEARxAG" ~ "Year X Agr.",
                                    param_label == "YEARxURB" ~ "Year X Urb.",
                                    param_label == "YEARxNITxAG" ~ "Year X Nit. X Agr.",
                                    param_label == "YEARxNITxURB" ~ "Year X Nit. X Urb.",
                                    param_label == "YEARxAGxURB" ~ "Year X Agr. X Urb.",
                                    param_label == "YEARxNITxAGxURB" ~ "Year X Nit. X Agr. X Urb."), levels = c("Intercept", "Agric. Cover", "Urban Cover", "Nit. Dep.", 
                                                                                                                "Nit. X Agr.", "Nit. X Urb.", "Agr. X Urb.", "Nit. X Agr. X Urb.",
                                                                                                                "Year", "Year X Nit.", "Year X Agr.", "Year X Urb.", "Year X Nit. X Agr.", "Year X Nit. X Urb.", "Year X Agr. X Urb.", "Year X Nit. X Agr. X Urb.")),
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
  facet_wrap(~param_f, scales = "free_x", nrow = 4)+
  labs(x = "Posterior Est.", y = "Species")+
  guides(fill = "none")+
  scale_color_manual(values = species_colors)+
  scale_fill_manual(values = species_colors)+
  scale_x_continuous(labels = scales::label_number(), guide = guide_axis(check.overlap = TRUE))+
  theme_bw() + theme(axis.text.y = element_text(face = "italic"),
                     axis.text.x = element_text(size = rel(.8)))

# posterior_hist
ggsave(posterior_hist, filename = "Plots/posterior_hist_yearmodel.png", width = 10, height = 9)




effects_summary <- effects_df %>% 
  group_by(param, param_label, spp_label) %>% 
  summarize(mean = mean(value), 
            lwr = quantile(value, .025),
            upr = quantile(value, .975),
            prob_pos = sum(value>0)/1000,
            prob_neg = 1-prob_pos)
# write.csv(effects_summary, file = "Posterior_prob_results.csv")





################################################################################################################################
##########  # Trying out a different plot of trends to visualize the interactin terms better ######
################################################################################################################################

quantile_summary <- data %>% 
  group_by(Spp_code) %>% 
  summarize(min_year = min(year),
            max_year = max(year),
            min_ag = min(PercentAg),
            max_ag = max(PercentAg),
            min_urb = min(PercentUrban),
            max_urb = max(PercentUrban),
            min_nit = min(mean_TIN_10km),
            max_nit = max(mean_TIN_10km))
# creating a hull to remove points in the covariate space that we never observe
ag_urb_coords <- data %>% 
  select(-geometry) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c('PercentAg', 'PercentUrban'), remove = FALSE) 

ag_urb_hull.aghy <- st_convex_hull(st_union(ag_urb_coords %>% filter(Spp_code == "AGHY")))
ag_urb_hull.agpe <- st_convex_hull(st_union(ag_urb_coords %>% filter(Spp_code == "AGPE")))
ag_urb_hull.elvi <- st_convex_hull(st_union(ag_urb_coords %>% filter(Spp_code == "ELVI")))

# plot(ag_urb_hull.aghy)
# agr. by urban interaction
preddata.aghy <- expand.grid(Spp_code = "AGHY",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_ag, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_ag, length.out = 100),
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_urb, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_urb, length.out = 100),
                             mean_TIN_10km = 0)  %>% 
  st_as_sf(coords = c("PercentAg", "PercentUrban"), remove = FALSE)  
aghy_within <- st_within( preddata.aghy, ag_urb_hull.aghy, sparse = FALSE)
preddata.aghy <- preddata.aghy[aghy_within,]


preddata.agpe <- expand.grid(Spp_code = "AGPE",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_ag, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_ag, length.out = 100),
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_urb, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_urb, length.out = 100),
                             mean_TIN_10km = 0) %>% 
  st_as_sf(coords = c("PercentAg", "PercentUrban"), remove = FALSE)  
agpe_within <- st_within( preddata.agpe, ag_urb_hull.agpe, sparse = FALSE)
preddata.agpe <- preddata.agpe[agpe_within,]


preddata.elvi <- expand.grid(Spp_code = "ELVI",
                             year = c(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_year, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_ag, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_ag, length.out = 100),
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_urb, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_urb, length.out = 100),
                             mean_TIN_10km = 0) %>% 
  st_as_sf(coords = c("PercentAg", "PercentUrban"), remove = FALSE)  
elvi_within <- st_within( preddata.elvi, ag_urb_hull.elvi, sparse = FALSE)
preddata.elvi <- preddata.elvi[elvi_within,]


preddata <- bind_rows(preddata.aghy, preddata.agpe, preddata.elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         year_label = case_when(year >=0 ~ "max_year",
                               year < 0 ~ "min_year")) %>% 
  as_tibble() %>% select(-geometry)


# gennerating predictiion
year.pred <- predict(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + data_summary$year,
         PercentUrban = PercentUrban + data_summary$PercentUrban,
         PercentAg= PercentAg + data_summary$PercentAg,
         mean_TIN_10km= mean_TIN_10km + data_summary$mean_TIN_10km)

year_diff <- year.pred %>% 
  select(Spp_code, species, year, year_label, PercentAg, PercentUrban, mean_TIN_10km, mean) %>% 
  pivot_wider(id_cols = c(Spp_code, species, PercentAg, PercentUrban, mean_TIN_10km), names_from = year_label, values_from = mean, names_prefix = "pred.") %>% 
  mutate(diff = pred.max_year - pred.min_year)



ag_urb_plot <- ggplot(year_diff)+
  geom_raster(aes(x = PercentAg, y = PercentUrban, fill = diff )) +
  # geom_point(data = endo_herb, aes(x = PercentAg,y = PercentUrban), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species)+
  scale_fill_viridis_c()+
  theme_classic()
ag_urb_plot

ggsave(ag_urb_plot, filename = "Plots/ag_urb_interaction_plot.png", width = 6)


######### agr. by nit interaction #####

ag_nit_coords <- data %>% 
  select(-geometry) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c('PercentAg', 'mean_TIN_10km'), remove = FALSE) 

ag_nit_hull.aghy <- st_convex_hull(st_union(ag_nit_coords %>% filter(Spp_code == "AGHY")))
ag_nit_hull.agpe <- st_convex_hull(st_union(ag_nit_coords %>% filter(Spp_code == "AGPE")))
ag_nit_hull.elvi <- st_convex_hull(st_union(ag_nit_coords %>% filter(Spp_code == "ELVI")))

preddata.aghy <- expand.grid(Spp_code = "AGHY",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_ag, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_ag, length.out = 100),
                             PercentUrban = 0,
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_nit, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentAg", "mean_TIN_10km"), remove = FALSE)  
aghy_within <- st_within( preddata.aghy, ag_nit_hull.aghy, sparse = FALSE)
preddata.aghy <- preddata.aghy[aghy_within,]


preddata.agpe <- expand.grid(Spp_code = "AGPE",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_ag, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_ag, length.out = 100),
                             PercentUrban = 0,
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_nit, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentAg", "mean_TIN_10km"), remove = FALSE)  
agpe_within <- st_within( preddata.agpe, ag_nit_hull.agpe, sparse = FALSE)
preddata.agpe <- preddata.agpe[agpe_within,]


preddata.elvi <- expand.grid(Spp_code = "ELVI",
                             year = c(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_year, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_year),
                             PercentAg = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_ag, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_ag, length.out = 100),
                             PercentUrban = 0,
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_nit, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentAg", "mean_TIN_10km"), remove = FALSE)  
elvi_within <- st_within( preddata.elvi, ag_nit_hull.elvi, sparse = FALSE)
preddata.elvi <- preddata.elvi[elvi_within,]


preddata <- bind_rows(preddata.aghy, preddata.agpe, preddata.elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         year_label = case_when(year >=0 ~ "max_year",
                                year < 0 ~ "min_year")) %>% 
  as_tibble() %>% select(-geometry)



# gennerating predictiion
year.pred <- predict(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + data_summary$year,
         PercentUrban = PercentUrban + data_summary$PercentUrban,
         PercentAg= PercentAg + data_summary$PercentAg,
         mean_TIN_10km= mean_TIN_10km + data_summary$mean_TIN_10km)

year_diff <- year.pred %>% 
  select(Spp_code, species, year, year_label, PercentAg, PercentUrban, mean_TIN_10km, mean) %>% 
  pivot_wider(id_cols = c(Spp_code, species, PercentAg, PercentUrban, mean_TIN_10km), names_from = year_label, values_from = mean, names_prefix = "pred.") %>% 
  mutate(diff = pred.max_year - pred.min_year)



ag_nit_plot <- ggplot(year_diff)+
  geom_raster(aes(x = mean_TIN_10km, y = PercentAg, fill = diff )) +
  # geom_point(data = endo_herb, aes(x = mean_TIN_10km,y = PercentAg), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species)+
  scale_fill_viridis_c()+
  theme_classic()
ag_nit_plot
ggsave(ag_nit_plot, filename = "Plots/ag_nit_interaction_plot.png", width = 6)


######### urb. by nit interaction #####

urb_nit_coords <- data %>% 
  select(-geometry) %>% 
  as_tibble() %>% 
  st_as_sf(coords = c('PercentUrban', 'mean_TIN_10km'), remove = FALSE) 

urb_nit_hull.aghy <- st_convex_hull(st_union(urb_nit_coords %>% filter(Spp_code == "AGHY")))
urb_nit_hull.agpe <- st_convex_hull(st_union(urb_nit_coords %>% filter(Spp_code == "AGPE")))
urb_nit_hull.elvi <- st_convex_hull(st_union(urb_nit_coords %>% filter(Spp_code == "ELVI")))

preddata.aghy <- expand.grid(Spp_code = "AGHY",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_year),
                             PercentAg = 0,
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_urb, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_urb, length.out = 100),
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="AGHY",]$min_nit, quantile_summary[quantile_summary$Spp_code=="AGHY",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentUrban", "mean_TIN_10km"), remove = FALSE)  
aghy_within <- st_within( preddata.aghy, urb_nit_hull.aghy, sparse = FALSE)
preddata.aghy <- preddata.aghy[aghy_within,]


preddata.agpe <- expand.grid(Spp_code = "AGPE",
                             year = c(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_year, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_year),
                             PercentAg = 0,
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_urb, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_urb, length.out = 100),
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="AGPE",]$min_nit, quantile_summary[quantile_summary$Spp_code=="AGPE",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentUrban", "mean_TIN_10km"), remove = FALSE)  
agpe_within <- st_within( preddata.agpe, urb_nit_hull.agpe, sparse = FALSE)
preddata.agpe <- preddata.agpe[agpe_within,]


preddata.elvi <- expand.grid(Spp_code = "ELVI",
                             year = c(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_year, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_year),
                             PercentAg = 0,
                             PercentUrban = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_urb, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_urb, length.out = 100),
                             mean_TIN_10km = seq(quantile_summary[quantile_summary$Spp_code=="ELVI",]$min_nit, quantile_summary[quantile_summary$Spp_code=="ELVI",]$max_nit, length.out = 100)) %>% 
  st_as_sf(coords = c("PercentUrban", "mean_TIN_10km"), remove = FALSE)  
elvi_within <- st_within( preddata.elvi, urb_nit_hull.elvi, sparse = FALSE)
preddata.elvi <- preddata.elvi[elvi_within,]


preddata <- bind_rows(preddata.aghy, preddata.agpe, preddata.elvi) %>% 
  mutate(collector_index = 9999, scorer_index = 9999,
         species = case_when(Spp_code == "AGHY" ~ species_names[1],
                             Spp_code == "AGPE" ~ species_names[2],
                             Spp_code == "ELVI" ~ species_names[3]),
         year_label = case_when(year >=0 ~ "max_year",
                                year < 0 ~ "min_year")) %>% 
  as_tibble() %>% select(-geometry)



# gennerating predictiion
year.pred <- predict(
  fit.year,
  newdata = preddata,
  formula = ~ invlogit(fixed),#+ collector_eval(collector_index) + scorer_eval(scorer_index)),
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  n.samples = 100) %>% 
  mutate(year = year + data_summary$year,
         PercentUrban = PercentUrban + data_summary$PercentUrban,
         PercentAg= PercentAg + data_summary$PercentAg,
         mean_TIN_10km= mean_TIN_10km + data_summary$mean_TIN_10km)

year_diff <- year.pred %>% 
  select(Spp_code, species, year, year_label, PercentAg, PercentUrban, mean_TIN_10km, mean) %>% 
  pivot_wider(id_cols = c(Spp_code, species, PercentAg, PercentUrban, mean_TIN_10km), names_from = year_label, values_from = mean, names_prefix = "pred.") %>% 
  mutate(diff = pred.max_year - pred.min_year)



urb_nit_plot <- ggplot(year_diff)+
  geom_raster(aes(x = mean_TIN_10km, y = PercentUrban, fill = diff )) +
  # geom_point(data = endo_herb, aes(x = mean_TIN_10km,y = PercentAg), alpha = .1)+
  # coord_sf()+
  facet_wrap(~species)+
  scale_fill_viridis_c()+
  theme_classic()
urb_nit_plot
ggsave(urb_nit_plot, filename = "Plots/urb_nit_interaction_plot.png", width = 6)



interaction_plot <- ag_urb_plot + ag_nit_plot + urb_nit_plot  + plot_layout(nrow = 3, guides = "collect")
ggsave(interaction_plot, filename = "Plots/interaction_plot.png")


