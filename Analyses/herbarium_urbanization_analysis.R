# Purpose: Fits spatial model of change in Epichloe endophyte prevalence in herbarium specimens and assesses effect of different anthropogenic land-use on trends.
# Authors: Mallory Tucker, and Joshua Fowler
# Updated: Feb 13, 2024

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(inlabru)
library(fmesher)

library(sf)
# library(rmapshaper)
library(terra)
library(tidyterra)
# library(ggpolypath)



library(patchwork)
library(ggmap)
library(ROCR)

invlogit<-function(x){exp(x)/(1+exp(x))}
species_colors <- c("#1b9e77","#d95f02","#7570b3")
endophyte_colors <- c("#fdedd3","#f3c8a8", "#5a727b", "#4986c7", "#181914",  "#163381")


################################################################################
############ Read in the herbarium dataset ############################### 
################################################################################

endo_herb_georef <- read_csv(file = "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/DigitizedHerbariumRecords/endo_herb_georef.csv") %>% 
  # filter(Country != "Canada") %>%
  mutate(sample_temp = Sample_id) %>% 
  separate(sample_temp, into = c("Herb_code", "spp_code", "specimen_code", "tissue_code")) %>% 
  mutate(species_index = as.factor(case_when(spp_code == "AGHY" ~ "1",
                                             spp_code == "AGPE" ~ "2",
                                             spp_code == "ELVI" ~ "3"))) %>% 
  mutate(species = case_when(spp_code == "AGHY" ~ "A. hyemalis",
                             spp_code == "AGPE" ~ "A. perennans",
                             spp_code == "ELVI" ~ "E. virginicus")) %>% 
  mutate(decade_bin = floor(year/10)*10) %>% 
  mutate(std_year = (year-mean(year, na.rm = T))/sd(year, na.rm = T))



endo_herb <- endo_herb_georef %>% 
  filter(!is.na(Endo_status_liberal)) %>%
  filter(!is.na(spp_code)) %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lon>-110 ) %>% 
  filter(Country != "Canada" ) %>% 
  mutate(year_bin = case_when(year<1970 ~ "pre-1970",
                              year>=1970 ~ "post-1970")) %>% 
  mutate(scorer_id = case_when(scorer_id == "BellaGuttierez" ~ "BellaGutierrez",
                               TRUE ~ scorer_id)) %>% 
  mutate(endo_status_text = case_when(Endo_status_liberal == 0 ~ "E-",
                                      Endo_status_liberal == 1 ~ "E+"))  %>% 
  filter(spp_code %in% c("AGHY", "AGPE", "ELVI")) 

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



# Loading in contemporary survey data for AGHY and ELVI, which we will use for model validation
# contemp_surveys <- read_csv(file = "contemp_surveys.csv") %>% 
#   mutate(decade_bin = floor(Year/10)*10) %>% 
#   filter(SpeciesID == "AGHY") %>% 
#   st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>% 
#   st_transform(epsg6703km) %>% 
#   mutate(
#     easting = st_coordinates(.)[, 1],
#     northing = st_coordinates(.)[, 2]
#   ) 
# 
# 
# contemp_random_sample <- read_csv(file = "contemp_random_sample.csv")


# register_google(key = "")
# map <- ggmap::get_map(zoom = 3, source = "google", maptype = c("satellite"))

outline_map <- map_data("world")
states_shape <- map_data("state")
counties_shape <- map_data("county")
# ggplot()+
#   geom_map(data = counties_shape, map = counties_shape, aes(long, lat, map_id = region))

collections_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes(map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = spp_code), alpha = .7, size = 1.2)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = species_colors)+
  theme_light()+
  theme(legend.text = element_text(face = "italic"))+
  labs(x = "Longitude", y = "Latitude", color = "Host Species")

collections_map
# ggsave(collections_map, filename = "collections_map.png", width = 7, height = 4)

endo_status_map <- ggplot()+
  geom_map(data = outline_map, map = outline_map, aes( map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
  geom_map(data = states_shape, map = states_shape, aes(map_id = region), color = "grey", linewidth = .1, fill = NA)+
  geom_point(data = endo_herb, aes(x = lon, y = lat, color = endo_status_text), alpha = .7, size = .8)+
  facet_wrap(~factor(year_bin, levels = c("pre-1970", "post-1970"))+species)+
  coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
  scale_color_manual(values = c(endophyte_colors[2],endophyte_colors[6]))+
  theme_light()+
  theme(strip.background = element_blank(),
        strip.text  = element_text(face = "italic", color = "black"))+
  labs(x = "Longitude", y = "Latitude", color = "Endophyte Status")
endo_status_map

# ggsave(endo_status_map, filename = "endo_status_map.png", width = 10, height = 5)



# testing if the trends are still are present without oldest samples
# endo_herb <- endo_herb_georef %>% 
#   filter(!is.na(Endo_status_liberal_1)) %>%
#   filter(!is.na(Spp_code)) %>% 
#   filter(!is.na(lon) & !is.na(year)) %>% 
#   filter(lon>-110) %>% 
#   filter(year>1910)

# summary_endo_herb <- endo_herb %>% 
#   mutate(score_match = case_when(Endo_status_conservative == Endo_status_liberal ~ 1,
#                                  Endo_status_conservative != Endo_status_liberal ~ 0)) %>% 
#   group_by(Spp_code) %>%
#   summarize(n  = n(),
#             n_match = sum(score_match),
#             percent_match = n_match/n)
# 
mode <- function(codes){
  which.max(tabulate(codes))
}

summary_endo_herb <- endo_herb %>% 
  mutate(Sample_id_temp = Sample_id) %>% 
  separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
  # filter(seed_scored>0) %>% 
  filter(month<=12&month>0) %>% 
  group_by(species) %>% 
  summarize(n(),
            avg_seed = mean(seed_scored, na.rm = T),
            avg_month = mode(as.numeric(month)))



# making separate dataframes for each species

endo_herb_AGHY <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "AGHY") %>% 
  filter(!is.na(lon) & !is.na(year))

endo_herb_AGPE <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "AGPE") %>% 
  filter(!is.na(lon) & !is.na(year)) %>% 
  filter(lat>=20) # dropping one plant in mexico


endo_herb_ELVI <- endo_herb %>% 
  filter(!is.na(Endo_status_liberal)) %>% 
  filter(spp_code == "ELVI") %>% 
  filter(!is.na(lon) & !is.na(year)) 


endo_herb <- endo_herb_AGHY
# endo_herb <- endo_herb_AGPE 
# endo_herb <- endo_herb_ELVI



##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
coords <- cbind(endo_herb$easting, endo_herb$northing)
# AGHY_coords <- cbind(endo_herb_AGHY$lon, endo_herb_AGHY$lat)
# AGPE_coords <- cbind(endo_herb_AGPE$lon, endo_herb_AGPE$lat)
# ELVI_coords <- cbind(endo_herb_ELVI$lon, endo_herb_ELVI$lat)

# non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.1, resolution = c(100, 100))


non_convex_bdry <- fm_extensions(
  coords,
  convex = c(250, 500),
  concave = c(250, 500)
)

# plot(non_convex_bdry[[1]])


non_convex_bdry[[1]] <- st_sfc(non_convex_bdry[[1]], crs = epsg6703km) 
non_convex_bdry[[2]] <- st_sfc(non_convex_bdry[[2]], crs = epsg6703km) 




coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE))) %>% st_transform(epsg6703km)
# plot(coastline)




bdry <- st_intersection(coastline$geom, non_convex_bdry[[1]])

# plot(bdry)

bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
  as("Spatial")

non_convex_bdry[[1]] <- bdry_polygon



max.edge = diff(range(coords[,1]))/(50)


mesh <- fm_mesh_2d_inla(
  # loc = coords,
  boundary = non_convex_bdry, max.edge = c(max.edge, max.edge*2), # km inside and outside
  cutoff = 20,
  crs = fm_crs(endo_herb)
  # crs=CRS(proj4string(bdry_polygon))
) # cutoff is min edge
# plot it



# 
# # AGHY_non_convex_bdry <- inla.nonconvex.hull(AGHY_coords, -0.03, -0.05, resolution = c(100, 100))
# # AGPE_non_convex_bdry <- inla.nonconvex.hull(AGPE_coords, -0.04, -0.05, resolution = c(100, 100))
# # ELVI_non_convex_bdry <- inla.nonconvex.hull(ELVI_coords, -0.03, -0.05, resolution = c(100, 100))
# 
# sf::sf_use_s2(FALSE)
# bdry_st <- st_multipolygon(non_convex_bdry[[1]])
#                            mutate(easting = V1,  northing = V2) %>%
#                            st_as_sf(coords = c("easting", "northing"), crs = fm_crs(endo_herb)) %>%
#                            summarise(geometry = st_combine(geometry)) %>%
#                            st_cast("POLYGON"))
# # 
# # AGHY_bdry_st <- st_make_valid(as_tibble(AGHY_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# # AGPE_bdry_st <- st_make_valid(as_tibble(AGPE_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# # ELVI_bdry_st <- st_make_valid(as_tibble(ELVI_non_convex_bdry$loc)  %>% 
# #                                 mutate(lon = V1,  lat = V2) %>% 
# #                                 st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
# #                                 summarise(geometry = st_combine(geometry)) %>% 
# #                                 st_cast("POLYGON"))
# 
# coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = fm_crs(endo_herb)))
# # plot(coastline)
# 
# sf::st_transform(coastline, fm_crs(endo_herb))
# 
# 
# bdry <- st_intersection(coastline$geom, bdry_st)
# # AGHY_bdry <- st_intersection(coastline$geom, AGHY_bdry_st)
# # AGPE_bdry <- st_intersection(coastline$geom, AGPE_bdry_st)
# # ELVI_bdry <- st_intersection(coastline$geom, ELVI_bdry_st)
# 
# bdry_polygon <- st_cast(st_sf(bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
#   as("Spatial")
# # AGHY_bdry_polygon <- st_cast(st_sf(AGHY_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # AGPE_bdry_polygon <- st_cast(st_sf(AGPE_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # ELVI_bdry_polygon <- st_cast(st_sf(ELVI_bdry), "MULTIPOLYGON", group_or_split = TRUE) %>% st_union() %>% 
# #   as("Spatial")
# # plot(bdry_polygon)
# 
# max.edge = diff(range(coords[,2]))/(10)
# mesh10 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*max.edge,
#                        boundary = bdry_polygon,
#                        offset = c(1,4),
#                        cutoff = max.edge/(5))
# plot(mesh10)
# mesh5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2),
#                       boundary = bdry_polygon,
#                       offset = c(1,4),
#                       cutoff = max.edge/(5))
# plot(mesh5)
# # # 
# # # mesh2.5 <- inla.mesh.2d(loc = coords, max.edge = c(1,2)*(max.edge/2/2),
# # #                         boundary = bdry_polygon,
# # #                         offset = c(1,4),
# # #                         cutoff = max.edge/(10))
# # # # plot(mesh2.5)
# # 
# mesh1 <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/4),
#                       boundary = bdry_polygon,
#                       offset = c(1,4),
#                       cutoff = max.edge/(10))
# plot(mesh1)
# # 
# mesh <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/4),
#                      boundary = bdry_polygon,
#                      offset = c(1,4),
#                      cutoff = max.edge/(15),
#                      crs = 4326)

ggplot() +
  gg(data = mesh) +
  geom_point(data = endo_herb, aes(x = easting, y = northing, col = species), size = 1) +
  coord_sf()+
  theme_bw() +
  labs(x = "", y = "")


# make spde
# The priors from online tutorials are :   # P(practic.range < 0.05) = 0.01 # P(sigma > 1) = 0.01
# For ESA presentation, I used the following which at least "converged" but seem sensitive to choices
# for AGHY =  P(practic.range < 0.1) = 0.01 # P(sigma > 1) = 0.01
# for AGPE =  P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01
# for ELVI = P(practic.range < 1) = 0.01 # P(sigma > 1) = 0.01



# In general, the prior on the range of the spde should be bigger than the max edge of the mesh
prior_range <- max.edge*3
# the prior for the SPDE standard deviation is a bit trickier to explain, but since our data is binomial, I'm setting it to .5
prior_sigma <- 1

spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(prior_range, 0.05),
  
  prior.sigma = c(prior_sigma, 0.05),
  # constr = TRUE
)
# inlabru makes making weights spatial effects simpler because we don't have to make projector matrices for each effect. i.e we don't have to make an A-matrix for each spatially varying effect.
# this means we can go strat to making the components of the model


# This is the model formula with a spatial effect (spatially varying intercept) as well a spatially varying time slopes


# setting the random effects prior
pc_prec <- list(prior = "pcprec", param = c(1, 0.1))



# components

svc_components <- ~ Intercept(1)+ std_year +
  space_int(coords, model = spde) +   # can add group argument to make specific to each species
  collector(collector_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))+
scorer(scorer_factor, model = "iid", constr = TRUE, hyper = list(pc_prec))




# formula, with "." meaning "add all the model components":
svc_formula <- Endo_status_liberal ~ .



# Now run the model

fit <- bru(svc_components,
           like(
             formula = svc_formula,
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

# make mesh projector to get model summaries from the mesh to the mapping grid
mesh_proj <- fm_evaluator(
  mesh,
  xlim = xlim, ylim = ylim, dims = grd_dims
)



space_int <- data.frame(
  median = invlogit(fit$summary.random$space_int$"0.5quant"),
  range95 = invlogit(fit$summary.random$space_int$"0.975quant") -
    invlogit(fit$summary.random$space_int$"0.025quant")
)



# loop to get estimates on a mapping grid
pred_grids <- lapply(
  list(space_int = space_int),
  function(x) as.matrix(fm_evaluate(mesh_proj, x))
)


# make a terra raster stack with the posterior median and range95
out_stk <- rast()
for (j in 1:1) {
  mean_j <- cbind(expand.grid(x = mesh_proj$x, y = mesh_proj$y),
                  Z = c(matrix(pred_grids[[j]][, 1], grd_dims[1]))
  )
  mean_j <- rast(mean_j)
  range95_j <- cbind(expand.grid(X = mesh_proj$x, Y = mesh_proj$y),
                     Z = c(matrix(pred_grids[[j]][, 2], grd_dims[1]))
  )
  range95_j <- rast(range95_j)
  out_j <- c(mean_j, range95_j)
  terra::add(out_stk) <- out_j
}
names(out_stk) <- c(
  "space_median", "space_range95"
)
#Masking the raster to our boundary
out_stk <- terra::mask(out_stk,bdry_st)

make_plot_field <- function(data_stk, scale_label) {
  ggplot() +
    geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = "#FAF9F6")+
    geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "grey", linewidth = .1, fill = NA)+
    geom_sf(fill = NA) +
    coord_sf(xlim = c(-109,-68), ylim = c(21,49))+
    geom_spatraster(data = data_stk) +
    labs(x = "", y = "") +
    scale_fill_viridis_c(name = scale_label, option = "turbo", na.value = "transparent")+
    theme(text = element_text(size = 2))+
    theme_bw()
  }



# This is the plot of the spatial intercept
ps <- make_plot_field(
  data_stk = out_stk[["space_median"]],
  scale_label = "Spatial\nPosterior Mean"
)
ps

ps_r <- make_plot_field(
  data_stk = out_stk[["space_range95"]],
  scale_label = "Spatial\nPosterior 95 CI\n"
)
ps_r

# Taking alot of material from this blog post: https://inlabru-org.github.io/inlabru/articles/2d_lgcp_covars.html

# plotting the spde range and sd posteriors
spde.range <- spde.posterior(fit, "space_int", what = "range")
spde.logvar <- spde.posterior(fit, "space_int", what = "log.variance")
spde.var <- spde.posterior(fit, "space_int", what = "variance")

range.plot <- plot(spde.range)
var.plot <- plot(spde.var)


range.plot
var.plot

# and plot the matern covariance (our spatial decay effect)

cov.plot <- plot(spde.posterior(fit, "space_int", what = "matern.covariance"))
cov.plot



# Making a plot of the marginal posterior of the year slope. We can see the posterior has a small positive slope for AGHY
flist <- vector("list", NROW(fit$summary.fixed))
for (i in seq_along(flist)) {
  flist[[i]] <- plot(fit, rownames(fit$summary.fixed)[i])
}
multiplot(plotlist = flist, cols = 2)



###### Getting and plotting prediction from model #####
# Here I am showing just the trend across year, but you could presumably show different interactions
min_year <- min(endo_herb$std_year)
max_year <- max(endo_herb$std_year)
preddata <- data.frame(std_year = seq(min_year, max_year), x_mean = c(1,1,1,1,1,1,5))


preddata <- data.frame(std_year = seq(min_year, max_year))
preddata <- expand.grid(std_year = seq(min_year, max_year), X_mean = c(1,5))

year.pred <- predict(
  fit,
  newdata = preddata,
  formula = ~(invlogit( Intercept + std_year))) %>%  
    mutate(year = std_year*(sd(endo_herb$year)) + mean(endo_herb$year))


ggplot(year.pred) +
  # geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal))+
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  lims(y = c(0,1))

