# Purpose: Fits spatial model of change in Epichloe endophyte prevalence in herbarium specimens and assesses effect of different anthropogenic land-use on trends.
# Authors: Mallory Tucker, and Joshua Fowler
# Updated: Feb 13, 2024

library(tidyverse) # for data manipulation and ggplot
library(INLA) # for fitting integrated nested Laplace approximation models
library(INLAutils) # supposedly has a function to plot residuals
library(inlabru)
library(fmesher)

library(sf)
library(rmapshaper)
library(terra)
library(tidyterra)


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


Mallorypath <- ""
Joshpath <- "~/Dropbox/Josh&Tom - shared/Endo_Herbarium/"
path<-Joshpath


endo_herb_georef <- read_csv(file = paste0(path, "DigitizedHerbariumRecords/endo_herb_georef.csv")) %>%
  # filter(Country != "Canada") %>%
  mutate(Spp_code = case_when(grepl("AGHY", Sample_id) ~ "AGHY",
                              grepl("ELVI", Sample_id) ~ "ELVI",
                              grepl("AGPE", Sample_id) ~ "AGPE")) %>% 
  mutate(species_index = as.factor(case_when(Spp_code == "AGHY" ~ "1",
                             Spp_code == "AGPE" ~ "2",
                             Spp_code == "ELVI" ~ "3"))) %>% 
  mutate(species = case_when(Spp_code == "AGHY" ~ "A. hyemalis",
                             Spp_code == "AGPE" ~ "A. perennans",
                             Spp_code == "ELVI" ~ "E. virginicus")) %>% 
  mutate(decade = floor(year/10)*10)

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



collector_levels <- levels(as.factor(endo_herb$collector_string))
collector_no <- paste0("Collector",1:nlevels(as.factor(endo_herb$collector_string)))

endo_herb$collector_id <- collector_no[match(as.factor(endo_herb$collector_string), collector_levels)]




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

summary_endo_herb <- endo_herb %>% 
  mutate(Sample_id_temp = Sample_id) %>% 
  separate(Sample_id_temp, sep = "_", into = c("herbarium", "spp_code", "plant_no")) %>% select(-spp_code, -plant_no) %>% 
  # filter(seed_scored>0) %>% 
  filter(month<=12&month>0) %>% 
  group_by(species) %>% 
  summarize(n(),
            avg_seed = mean(seed_scored, na.rm = T),
            avg_month = mode(as.numeric(month)))


##########################################################################################
############ Setting up and running INLA model with inlabru ############################### 
##########################################################################################

##### Building a spatial mesh #####

# Build the spatial mesh from the coords for each species and a boundary around each species predicted distribution (eventually from Jacob's work ev)
coords <- cbind(endo_herb$lon, endo_herb$lat)


non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))

plot(non_convex_bdry)

sf::sf_use_s2(FALSE)
bdry_st <- st_make_valid(as_tibble(non_convex_bdry$loc)  %>% 
                           mutate(lon = V1,  lat = V2) %>% 
                           st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
                           summarise(geometry = st_combine(geometry)) %>% 
                           st_cast("POLYGON"))


coastline <- st_make_valid(sf::st_as_sf(maps::map("world", regions = c("usa", "canada", "mexico"), plot = FALSE, fill = TRUE), crs = 4326))
# plot(coastline)

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


mesh <- inla.mesh.2d(loc = coords, max.edge = c(.5,2)*(max.edge/2/2),
                     boundary = bdry_polygon,
                     offset = c(1,4),
                     cutoff = max.edge/(15))
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
spde <- inla.spde2.pcmatern(
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
s_components <- ~ -1 +
  year + species(endo_herb$species_index, model = "factor_full")+
  space_int(coords, group = species_index, model = spde)


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
# Here I am showing u
min_year <- min(endo_herb$year)
max_year <- max(endo_herb$year)

year.pred <- predict(
  fit,
  data.frame(year = seq(min_year, max_year)),
  formula = ~ invlogit(Intercept + year))


ggplot(year.pred) +
  geom_point(data =endo_herb, aes(x = year, y = Endo_status_liberal))+
  geom_line(aes(year, mean)) +
  geom_ribbon(aes(year, ymin = q0.025, ymax = q0.975), alpha = 0.2) +
  geom_ribbon(aes(year, ymin = mean - 1 * sd, ymax = mean + 1 * sd), alpha = 0.2) +
  lims(y = c(0,1))

