# Purpose: Fits spatial model of change in Epichloe endophyte prevalence in herbarium specimens and assesses effect of different anthropogenic land-use on trends.
# Authors: Mallory Tucker, and Joshua Fowler
# Updated: 

library(devtools)
library("devtools")
# devtools::install_github(repo = "https://github.com/hrue/r-inla", ref = "stable", subdir = "rinla", build = FALSE, force = TRUE)
#INLA relies on Rgraphviz (and other packages, you can use Bioconductor to help install)
library(dplyr)
library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
#devtools::install_github('timcdlucas/INLAutils')
library(INLAutils) # supposedly has a function to plot residuals, might not need?
library(inlabru)
library(fmesher)
library(sf)
library(rmapshaper)
library(terra)
library(tidyterra)
#issue with vctrs namespace, .0.60 is loaded but need updated version

library(patchwork)
library(ggmap)
library(ROCR)

invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################
# This is where I'm loading in the version of the data set without land cover and nitrogen data, but we ought to be able to replace this easily


Mallorypath <- "C:/Users/mnt4/Documents/sp23bios310/"
Joshpath <- "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/"
path <- Mallorypath


endo_herb_georef <- read_csv(file = ("C:/Users/malpa/Box/endo_urb/Zonal_hist_nlcd.csv")) %>%
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
         TotalPixels = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24 + HISTO_0 + HISTO_11 +HISTO_12 + HISTO_31 +HISTO_41 + HISTO_42 + HISTO_43 + HISTO_52 + HISTO_71+ HISTO_90 + TotalAg + HISTO_95,
         TotalDeveloped = HISTO_21 + HISTO_22 + HISTO_23 + HISTO_24,
         OtherLC = (TotalPixels - (TotalAg + TotalDeveloped))/TotalPixels *100,
         PercentUrban = TotalDeveloped/TotalPixels * 100,
         PercentAg = TotalAg/TotalPixels * 100)

# Doing some filtering to remove NA's and some data points that probably aren't accurate species id's
endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(Spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" ) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+"))

#loading in nitrogen data too (using N data from all years created with raster stack in merging_nitrogen script)
nit <- read.csv(file = "C:/Users/malpa/Box/endo_urb/NO3/NO3/endo_herb_nit.csv")
nitendoherb <- merge(nit, endo_herb)
endo_herb <- nitendoherb

#remove scorer 26 and mean center predictor variables
endo_herb <- endo_herb %>%
  mutate(std_year = (year-mean(year, na.rm = T))) %>%  
  mutate(std_NO3 = (NO3_mean-mean(NO3_mean, na.rm = T)))%>% 
  mutate(std_Ag = (PercentAg-mean(PercentAg, na.rm = T)))%>%
  mutate(std_Urb = (PercentUrban-mean(PercentUrban, na.rm = T)))%>%
  filter(scorer_factor != "Scorer26")


write.csv(endo_herb, file = "EndoHerb_withNitrogen.csv")

# mutate(Endo_status_liberal = case_when(Spp_code == "AGPE" & Endo_status_liberal ==0 ~ 1,
#                                         Spp_code == "AGPE" & Endo_status_liberal == 1~ 0, TRUE ~ Endo_status_liberal))

# Looking at just AGHY for now, but we can likely fit one model for all species
# endo_herb <- endo_herb %>% 
#   filter(Spp_code == "AGHY")


# updating the scorer labels
scorer_levels <- levels(as.factor(endo_herb$scorer_id))
scorer_no <- paste0("Scorer",1:nlevels(as.factor(endo_herb$scorer_id)))

endo_herb$scorer_id <- scorer_no[match(as.factor(endo_herb$scorer_id), scorer_levels)]

# updating the collector labels


collector_levels <- levels(as.factor(endo_herb$collector_))
collector_no <- paste0("Collector",1:nlevels(as.factor(endo_herb$collector_)))

#Error in `$<-.data.frame`(`*tmp*`, collector_id, value = character(0)) : replacement has 0 rows, data has 1859
# endo_herb$collector_id <- collector_no[match(as.factor(endo_herb$collector_), collector_levels)]


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
  ) 

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


mode <- function(codes){
  which.max(tabulate(codes))
}

# summary_endo_herb <- endo_herb %>% 
#   mutate(Sample_id_temp = Sample_id) %>% 
#   separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
#   # filter(seed_scored>0) %>% 
#   filter(month<=12&month>0) %>% 
#   group_by(species) %>% 
#   dplyr::summarize(n(),
#             avg_seed = mean(seed_score, na.rm = T),
#             avg_month = mode(as.numeric(month)))


##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
coords <- cbind(endo_herb$lon, endo_herb$lat)


non_convex_bdry <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(100, 100))

plot(non_convex_bdry)

sf::sf_use_s2(FALSE)
bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
                           mutate(lon = V1,  lat = V2) %>% 
                           st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                           summarise(geometry = st_combine(geometry)) %>% 
                           st_cast("POLYGON"))


coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = 4326))
#plot(coastline)

st_crs(coastline) <- 4326


bdry <- st_intersection(coastline$geom, bdry_st)

bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")

plot(bdry_polygon)

max.edge = diff(range(coords[,2]))/10
# mesh10 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
#                        boundary = bdry_polygon,
#                        offset = c(1,4),
#                        cutoff = max.edge/(5))
# plot(mesh10)
# 
# 
# mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
#                       boundary = bdry_polygon,
#                       offset = c(1,4),
#                       cutoff = max.edge/(10))
# plot(mesh1)


mesh <- fmesher::fm_mesh_2d_inla(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                     boundary = bdry_polygon,
                     offset = c(1,4),
                     cutoff = max.edge/(10))
# plot(mesh1)


ggplot() +
  gg(data = mesh) +
  geom_point(data = endo_herb, aes(x = lon, y = lat, col = species), size = 1) +
  coord_sf()+
  theme_bw() +
  labs(x = "", y = "")


# make spde (stochastic partial differential equation)

# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- .5

# The priors from online tutorials are :   # P(practic.range < 0.05) = 0.01 # P(sigma > 1) = 0.01
# For ESA presentation, I used the following which at least "converged" but seem sensitive to choices
# for AGHY =  P(practic.range < 0.1) = 0.01 # P(sigma > 1) = 0.01
# for AGPE =  P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
# for ELVI = P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
spde <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.1),
  
  prior.sigma = c(prior_sigma, 0.1)
)

# inlabru makes making spatial effects simpler compared to "INLA" because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# This is the model formula with a spatial effect (spatially varying intercept). To this, we can add predictor variables

# formula
# version for each species separately
# s_components <- ~ Intercept(1) +
#   year +
#   space_int(coords, model = spde)

# version with all species in one model. Note that we remove the intercept, and then we have to specify that the species is a factor 
s_components <- ~ Intercept(1) +
   PercentUrban + year + X_mean  + PercentAg +
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

fit$dic$dic
fit$mode$mode.status # a 0 or low value indicates "convergence"

fit$summary.fixed
fit$summary.random

saveRDS(fit, file = "fit_AGHY.rds")

fit <- read_rds(file = "fit_AGHY.rds")

saveRDS(fit, file = "fit_AGPE.rds")

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
  X_mean + yearEffect(year, model = "linear") + yearXnitEffect(year*X_mean, model = "linear") +
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


endo_herb <- endo_herb %>% mutate_at(c('PercentUrban', 'PercentAg', 'X_mean', "year"), ~(scale(.) %>% as.vector))

min_ag <- min(endo_herb$PercentAg, na.rm = TRUE)
max_ag <- max(endo_herb$PercentAg, na.rm = TRUE)
min_nit <- min(endo_herb$X_mean, na.rm = TRUE)
max_nit <- max(endo_herb$X_mean, na.rm = TRUE)

all_preddata <- expand.grid(PercentUrban = seq(min_urb, max_urb), PercentAg = mean(endo_herb$PercentAg, na.rm = TRUE), X_mean = c(min_nit, max_nit), year = NA)

#all_preddata <- data.frame(PercentUrban = seq(min_urb, max_urb), PercentAg = seq(min_ag, max_ag, length.out = 100), X_mean = seq(min_nit, max_nit, length.out =100), = seq(min_year, max_year, length.out = 100))

all.pred <- predict(
  fit,
  newdata = all_preddata,
  formula = ~ invlogit(Intercept +PercentUrban + X_mean + year))


ggplot(all.pred) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
  # geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal), shape = "|")+
  geom_line(aes(PercentUrban, mean, color = X_mean, group = X_mean)) +
 # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(PercentUrban, ymin = mean - 1 * sd, ymax = mean + 1 * sd, group = X_mean, fill = X_mean), alpha = 0.2) 
  lims(y = c(0,1))

  

#######################################################################################################################
################################################## straightening this up ##############################################
#######################################################################################################################
#creating full model
#and then model with x_mean and ag together


#4.22.24, fixed effect thing josh sent me  
  #removed: + PercentAg*X_mean,
s_components <- ~ my_fixed_effects(
  main =  std_Ag * std_year + std_Urb*std_year,
  model = "fixed")+
  space_int(coords, model = spde) 


# formula, with "." meaning "add all the model components":
s_formula <- Endo_status_liberal ~ .


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
min_urb <- min(endo_herb$std_Urb, na.rm = TRUE)
max_urb <- max(endo_herb$std_Urb, na.rm = TRUE)

min_ag <- min(endo_herb$std_Ag, na.rm = TRUE)
max_ag <- max(endo_herb$std_Ag, na.rm = TRUE)

min_nit <- min(endo_herb$std_NO3, na.rm = TRUE)
max_nit <- max(endo_herb$std_NO3, na.rm = TRUE)

#prediction dataframe. removed: X_mean = seq(min_nit, max_nit)
allpreddata <- expand.grid(PercentAg = seq(min_ag, max_ag), PercentUrban= c(min_urb, max_urb), year = c(1900,2000))
#which goes into this:
allpreddata <- predict(
  fit,
  newdata = allpreddata,
  formula = ~ invlogit(Intercept + my_fixed_effects)) 
# allpreddata <- allpreddata %>%
#   mutate(std_year = as.factor(std_year))


#now plot:
# ggplot(allpreddata) +
#   # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
#  # geom_point(data =endo_herb, aes(x = PercentAg, y = Endo_statu), shape = "|")+
#   geom_line(data = allpreddata, aes(PercentAg, color = as.factor(year), group = as.factor(year))) +
#   # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
#   geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975, fill = as.factor(year), group = as.factor(year), alpha = 0.2))
#  # facet_wrap(~ year)
# #lims(y = c(0,1))
  

########################################################### loop attempt #############################################

# making separate dataframes for each species

#creating list for for loop


species_codes <- c("AGHY", "AGPE", "ELVI")

#dataframe for output
pred.list <- list()


#explicitly writing out model inputs
ModMat_endo_herb <- endo_herb %>% 
  mutate(std_Ag_std_NO3 = std_Ag*std_NO3,
         std_Ag_std_year = std_Ag*std_year,
         std_NO3_std_year = NO3_mean*std_year,
         std_Ag_std_year_std_NO3 = std_Ag*std_year*std_NO3)

#removing rows with nas for the interaction terms
endo_herb <- ModMat_endo_herb
endo_herb <- endo_herb[complete.cases(endo_herb$std_Ag_std_NO3), ]
endo_herb <- endo_herb[complete.cases(endo_herb$std_Ag.std_year),]
endo_herb <- endo_herb[complete.cases(endo_herb$std_Ag_std_year),]
endo_herb <- endo_herb[complete.cases(endo_herb$std_Ag_std_year_std_NO3),]

#creating data frames for specific species
endo_herb_AGHY <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year))

endo_herb_AGPE <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lat>=20) # dropping one plant in mexico

endo_herb_ELVI <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(Spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 

endo_herb_list <- list()
endo_herb_list[[1]] <- endo_herb_AGHY
endo_herb_list[[2]] <- endo_herb_AGPE
endo_herb_list[[3]] <- endo_herb_ELVI

# pred.list[[1]] <- 
# pred.list[[2]]
# pred.list[[3]]
#changed out s_components
s_components <- ~ Intercept(1) +
  std_Ag + std_year + std_NO3 + std_Ag_std_year+std_Ag_std_NO3 +std_NO3_std_year +
  space_int(coords, model = spde) 

s_formula <- Endo_status_liberal ~ .


for (s in 1:3){
  data <- endo_herb_list[[1]]
  
  # Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
  coords <- cbind(data$lon, data$lat)
  
  
  non_convex_bdry <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(100, 100))
  
  plot(non_convex_bdry)
  
  sf::sf_use_s2(FALSE)
  bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
                             mutate(lon = V1,  lat = V2) %>% 
                             st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                             summarise(geometry = st_combine(geometry)) %>% 
                             st_cast("POLYGON"))
  
  
  coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = 4326))
  #plot(coastline)
  
  st_crs(coastline) <- 4326
  
  
  bdry <- st_intersection(coastline$geom, bdry_st)
  
  bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
    as("Spatial")
  
  plot(bdry_polygon)
  
  max.edge = diff(range(coords[,2]))/10

  
  mesh <- fmesher::fm_mesh_2d_inla(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                                   boundary = bdry_polygon,
                                   offset = c(1,4),
                                   cutoff = max.edge/(10))
  
  min_urb <- min(data$std_Urb, na.rm = TRUE)
  max_urb <- max(data$std_Urb, na.rm = TRUE)
  
  min_ag <- min(data$std_Ag, na.rm = TRUE)
  max_ag <- max(data$std_Ag, na.rm = TRUE)
  
  min_nit <- min(data$std_NO3, na.rm = TRUE)
  max_nit <- max(data$std_NO3, na.rm = TRUE)
  allpreddata <- expand.grid(std_Ag = c(min_ag, max_ag), std_NO3= seq(min_nit, max_nit, length.out = 10),std_Urb= c(min_urb,max_urb),Spp_code = species_codes[s], std_year = c(1900,2000))%>%
    mutate(std_Ag_std_NO3 = std_Ag*std_NO3,
           std_Ag_std_year = std_Ag*std_year,
           std_NO3_std_year = std_NO3*std_year,
           std_Ag_std_year_std_NO3 = std_Ag*std_year*std_NO3)

 

  # gennerating predictions and back-transforming the standardized year variable


 # Now run the model
  fit <- bru(
    s_components,
    like(
      formula = s_formula,
      family = "binomial",
      Ntrials = 1,
      data = data
    ),
    options = list(
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
      control.inla = list(int.strategy = "eb"),
      verbose = TRUE))
    pred.list <- list()
    pred.list[[1]] <- predict(
      fit,
      newdata = allpreddata,
      formula = ~ invlogit(Intercept + std_Ag + std_year + std_NO3 + std_Ag_std_year+std_Ag_std_NO3 +std_NO3_std_year),
      transform = TRUE)
    
}

ggplot(pred.list) +
  # geom_point(data = endo_herb_binned, aes(x = mean_year, y = mean_endo, size = sample))+
 # geom_point(data =endo_herb, aes(x = PercentAg, y = Endo_statu), shape = "|")+
  geom_line(data = pred.list,aes(x = std_NO3, y= mean, color = as.factor(std_Ag), group = as.factor(std_Ag))) +
  # geom_ribbon(aes(PercentAg, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  #geom_ribbon(aes(std_Ag, ymin = q0.025, ymax = q0.975, fill = as.factor(std_year), group = as.factor(std_year), alpha = 0.2))
  facet_wrap(~ std_year)
#lims(y = c(0,1))



coefvalues <- list(
  mean = pred.list[[1]]$summary.fixed$mean,
  sd = pred.list[[1]]$summary.fixed$sd
)

# Convert the extracted data to a data frame
coefvalues <- data.frame(coefvalues)

# Compute upper and lower bounds for error bars
coefvalues$upper_bound <- coefvalues$mean + coefvalues$sd
coefvalues$lower_bound <- coefvalues$mean - coefvalues$sd

#plot

coefplot <- barplot(coefvalues$mean, 
                     main = "Coefficient Plot with Error Bars", 
                     xlab = "Coefficient Value", ylab = "Effect Size", 
                     names.arg = c("Int", "Ag", "yr", "NO3", "Agyr","AgNO3",
                                           "NO3yr"), 
                     xlim = range(-coefvalues$lower_bound, coefvalues$upper_bound), 
                     col = "purple",
                     horiz = TRUE,
                     width = 0.1,
)
segments(coefplot, x0= coefvalues$mean + coefvalues$sd * 2, coefplot,
       x1 = coefvalues$mean - coefvalues$sd * 2, lwd = 1.5)


#instead
ggplot(data = coefvalues, aes(x = mean, y = c("Int", "Ag", "yr", "NO3", "Agyr","AgNO3","NO3") +
  geom_point() +
  geom_errorbarh(aes(xmin = coefvalues$lower, xmax = coefvalues$upper), height = 0.2) +
  theme_minimal() +
  labs(title = "coefplot", x = "Estimate", y = "Term")))

############################################################################################################
############################################ model w/o NO3 ########################################
##########################################################################################################
# making separate dataframes for each species

#creating list for for loop

allpreddata_b <- expand.grid(std_Ag = c(min_ag, max_ag),std_Urb= c(min_urb,max_urb),Spp_code = species_codes[s], std_year = c(1900,2000))%>%
  mutate(std_Ag_std_urb = std_Ag*std_urb,
         std_Ag_std_year = std_Ag*std_year,
         std_urb_std_year = std_urb*std_year,
         std_Ag_std_year_std_urb = std_Ag*std_year*std_urb)

#changed out s_components
b_components <- ~ Intercept(1) +
  std_Ag + std_year + std_urb + std_Ag_std_year+std_Ag_std_urb +std_urb_std_year +
  space_int(coords, model = spde) 

b_formula <- Endo_status_liberal ~ .


for (s in 1:3){
  data <- endo_herb_list[[s]]
  
  # Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
  coords <- cbind(data$lon, data$lat)
  
  
  non_convex_bdry <- fmesher::fm_nonconvex_hull_inla(coords, -0.03, -0.05, resolution = c(100, 100))
  
  plot(non_convex_bdry)
  
  sf::sf_use_s2(FALSE)
  bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
                             mutate(lon = V1,  lat = V2) %>% 
                             st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                             summarise(geometry = st_combine(geometry)) %>% 
                             st_cast("POLYGON"))
  
  
  coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = 4326))
  #plot(coastline)
  
  st_crs(coastline) <- 4326
  
  
  bdry <- st_intersection(coastline$geom, bdry_st)
  
  bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
    as("Spatial")
  
  plot(bdry_polygon)
  
  max.edge = diff(range(coords[,2]))/10
  
  
  mesh <- fmesher::fm_mesh_2d_inla(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                                   boundary = bdry_polygon,
                                   offset = c(1,4),
                                   cutoff = max.edge/(10))
  
  min_urb <- min(data$std_Urb, na.rm = TRUE)
  max_urb <- max(data$std_Urb, na.rm = TRUE)
  
  min_ag <- min(data$std_Ag, na.rm = TRUE)
  max_ag <- max(data$std_Ag, na.rm = TRUE)
  
  # min_nit <- min(data$std_NO3, na.rm = TRUE)
  # max_nit <- max(data$std_NO3, na.rm = TRUE)
  allpreddata_b <- expand.grid(std_Ag = c(min_ag, max_ag), std_urb= c(min_urb,max_urb),Spp_code = species_codes[s], std_year = c(1900,2000))%>%
    mutate(std_Ag_std_urb = std_Ag*std_urb,
           std_Ag_std_year = std_Ag*std_year,
           std_urb_std_year = std_urb*std_year,
           std_Ag_std_year_std_urb = std_Ag*std_year*std_urb)
  
  
  # gennerating predictions and back-transforming the standardized year variable
  
  
  # Now run the model
  b_fit <- bru(
    b_components,
    like(
      formula = s_formula,
      family = "binomial",
      Ntrials = 1,
      data = data
    ),
    options = list(
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
      control.inla = list(int.strategy = "eb"),
      verbose = TRUE))
  
  pred.list2[[s]] <- predict(
    b_fit,
    newdata = allpreddata_b,
    formula = ~ invlogit(Intercept + std_Ag + std_year + std_urb + std_Ag_std_year+std_Ag_std_urb +std_urb_std_year),
    transform = TRUE)
  
}


coefvalues2 <-  pred.list2(
  mean = b_fit$summary.fixed$mean,
  sd = b_fit$summary.fixed$sd
)

# Convert the extracted data to a data frame
coefvalues2 <- data.frame(coefvalues2)

# Compute upper and lower bounds for error bars
coefvalues2$upper_bound <- coefvalues2$mean + coefvalues2$sd
coefvalues2$lower_bound <- coefvalues2$mean - coefvalues2$sd

#plot

coefplot2 <- barplot(coefvalues2$mean, 
                    main = "Coefficient Plot with Error Bars", 
                    xlab = "Coefficient Value", ylab = "Effect Size", 
                    names.arg = c("Int", "Ag", "yr", "urb", "Agyr","AgUrb",
                                  "Urbyr"), 
                    xlim = range(-coefvalues2$lower_bound, coefvalues2$upper_bound), 
                    col = "purple",
                    horiz = TRUE,
                    width = 0.1,
)
segments(coefplot2, x0= coefvalues2$mean + coefvalues2$sd * 2, coefplot2,
         x1 = coefvalues2$mean - coefvalues2$sd * 2, lwd = 1.5)




# Calculate lower and upper bounds
coefvalues2$lower <- coefvalues2$mean - coefvalues$sd
coefvalues2$upper <- coefvalues2$mean + coefvalues$sd

# Create dot-and-whisker plot using ggplot2
ggplot(data = coefvalues2, aes(x = mean, y = c("Int", "Ag", "yr", "urb", "Agyr","AgUrb",
                                               "Urbyr"))) +
  geom_point() +
  geom_errorbarh(aes(xmin = coefvalues2$lower, xmax = coefvalues2$upper), height = 0.2) +
  theme_minimal() +
  labs(title = "coefplot", x = "Estimate", y = "Term")

