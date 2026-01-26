# Purpose: Visualizes covariates and correlations between covariates
# Authors: Mallory Tucker and Joshua Fowler
# Updated: Jan 16, 2025
library(dplyr)
library(tidyverse) # for data manipulation and ggplot

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
library(maps)




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



#########################################################################################
####################### Urban Map Figure ################################################
########################################################################################

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
  labs(x = "Longitude", y = "Latitude", fill = "% Urban")+
  theme_light()+
  theme(legend.text = element_text(face = "italic",size = 14),
        legend.title = element_text(size = 10)
  )+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
endo_urb_map

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
  theme(legend.text = element_text(face = "italic", size = 14),
        legend.title = element_text(size = 10)
  ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  labs(x = "Longitude", y = "Latitude", fill = "% Agricultural")
# endo_ag_map

#################################################################################
################################## nitrogen map #################################
################################################################################

#Join map outline data with endo_herb data
tiny_df_nit <- data.frame(endo_herb$region, endo_herb$mean_TIN_10km)
colnames(tiny_df_nit) <- c("region", "mean_TIN_10km")

counties_data_nit <- counties %>%
  left_join(tiny_df_nit, by = "region")

counties_data_nit <- na.omit(counties_data_nit)

# Create the choropleth map
endo_nit_map <- ggplot(data = counties_data_nit) +
  geom_map(data = outline_map, map = outline_map, aes(long, lat, map_id = region), color = "darkgray", linewidth = .3, fill = "#E8E7E7")+
  geom_map(data = states_shape, map = states_shape, aes(long, lat, map_id = region), color = "gray", linewidth = .2, fill = "#E8E7E7")+
  geom_polygon(aes(x = long, y = lat, group = group, fill = mean_TIN_10km), color = "darkgray", linewidth = .1) +
  # coord_fixed(ratio = 1) +
  coord_cartesian(xlim = c(-109, -68), ylim = c(25, 48))+
  scale_fill_gradient(low = "lightgray", high = "#BF00A0", na.value = NA)+
  theme_light() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.text = element_text(face = "italic",size = 14),
        legend.title = element_text(size = 10)
  ) +
  labs(x = "Longitude", y = "Latitude", fill = "Kg N/sqkm/year")
endo_nit_map
#Compile maps into one panel
#mapfig <-  endo_ag_map +endo_urb_map + endo_nit_map + plot_layout(ncol = 3) + plot_annotation(tag_levels = "A")
mapfig <-  endo_ag_map /endo_urb_map / endo_nit_map + plot_layout(ncol = 1) + plot_annotation(tag_levels = "A")

#Save map file
ggsave(mapfig, file = "Map_Figure.png", width = 5, height = 9)



################################################################################
############ looking at correlations between covariates ############################### 
################################################################################

predictors <- data.frame("PercentAg" = endo_herb$PercentAg,"PercentUrb" = endo_herb$PercentUrban, "NitrogenDeposition" = endo_herb$mean_TIN_10km,
                         "MeanTemp" = endo_herb$tmean_10km, "AnnualPrecip" = endo_herb$ppt_10km) %>% 
  na.omit()
pearson_correlations <- as_tibble(cor(predictors, method = "pearson"), rownames = "rownames") %>% 
  pivot_longer(cols = -rownames, names_to = "colnames", values_to = "correlation") %>% mutate(type = "Pearson Correlation")
spearman_correlations <- as_tibble(cor(predictors, method = "spearman"), rownames = "rownames") %>% 
  pivot_longer(cols = -rownames, names_to = "colnames", values_to = "correlation")%>% mutate(type = "Spearman Correlation")

correlations <- pearson_correlations %>% 
  bind_rows(spearman_correlations) %>% 
  mutate(correlation = case_when(rownames == "PercentAg" & colnames == "PercentUrb" | 
                                   rownames == "NitrogenDeposition" & colnames == "PercentUrb" |rownames == "NitrogenDeposition" & colnames == "PercentAg" | 
                                      rownames == "MeanTemp" & colnames == "NitrogenDeposition" | rownames == "MeanTemp" & colnames == "PercentAg" | rownames == "MeanTemp" & colnames == "PercentUrb" | 
                                        rownames == "AnnualPrecip" & colnames == "MeanTemp" | rownames == "AnnualPrecip" & colnames == "NitrogenDeposition" | rownames == "AnnualPrecip" & colnames == "PercentAg" | rownames == "AnnualPrecip" & colnames == "PercentUrb"   ~ NA,
                                 TRUE ~ correlation) )


predictor_correlations_pearson <- ggplot(filter(correlations, type == "Pearson Correlation"))+
  geom_tile(aes(x = colnames, y = rownames, fill = correlation))+
  geom_text(aes(x = colnames, y = rownames, label = round(correlation, 3)))+
  scale_fill_gradientn( colors = c("#998ec3",  "grey95", "#f1a340"), limits = c(-1,1), na.value = "white")+
  facet_wrap(~type)+
  labs(x = "", y = "", fill = "Correlation\n Coefficient")+
  coord_equal()+
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
predictor_correlations_spearman <- ggplot(filter(correlations, type == "Spearman Correlation"))+
  geom_tile(aes(x = colnames, y = rownames, fill = correlation))+
  geom_text(aes(x = colnames, y = rownames, label = round(correlation, 3)))+
  scale_fill_gradientn( colors = c("#998ec3",  "grey95", "#f1a340"), limits = c(-1,1), na.value = "white")+
  facet_wrap(~type)+
  labs(x = "", y = "", fill = "Correlation\n Coefficient")+
  coord_equal()+
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
predictor_correlations <- predictor_correlations_pearson + predictor_correlations_spearman + plot_annotation(tag_levels =  "A") + plot_layout(nrow = 2, guides = "collect")
ggsave(predictor_correlations, filename = "Plots/predictor_correlations.png", width = 7, height = 8)

# calculating variance inflation factors
vif.m <- lm(Endo_status_liberal~PercentAg +  PercentUrban + mean_TIN_10km + tmean_10km + ppt_10km + year, data = endo_herb)

vif <- car::vif(vif.m, type="terms")
vif_values <- tibble("vif" = vif, "term" = names(vif)) %>% 
  mutate(Predictor = case_when(term == "PercentAg" ~ "Percent Agr.",
                               term == "PercentUrban" ~ "Percent Urb.",
                               term == "mean_TIN_10km" ~ "Total Nit.",
                               term == "tmean_10km" ~ "Temp.",
                               term == "ppt_10km" ~ "PPT",
                               term == "year" ~ "Year"))

vif_plot <- ggplot(vif_values)+
  geom_point(aes(x = Predictor, y = vif))+
  geom_line(aes(x = Predictor, y = vif, group = 1))+
  ylim(0,2) + labs(y = "Variance Inflation Factor", x = "Predictor")+
  theme_minimal()

vif_plot
ggsave(vif_plot, filename = "Plots/vif_plot.png", width = 4, height = 4)




# now plotting "pairs" plots between covariates

ag_urb <- ggplot(endo_herb)+
  geom_point(aes(x = PercentAg, y = PercentUrban), alpha = .2)+
  labs(x = "Agr. Cover (%)", y = "Urban Cover (%)")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
ag_nit <- ggplot(endo_herb)+
  geom_point(aes(x = PercentAg, y = mean_TIN_10km), alpha = .2)+
  labs(x = "Agr. Cover (%)", y = "Nit. Dep (kg N/sqkm/year)")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
ag_temp <- ggplot(endo_herb)+
  geom_point(aes(x = PercentAg, y = tmean_10km), alpha = .2)+
  labs(x = "Agr. Cover (%)", y = "Mean Annual Temp.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
ag_ppt <- ggplot(endo_herb)+
  geom_point(aes(x = PercentAg, y = ppt_10km), alpha = .2)+
  labs(x = "Agr. Cover (%)", y = "Annual Precip.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))




urb_nit <- ggplot(endo_herb)+
  geom_point(aes(x = PercentUrban, y = mean_TIN_10km), alpha = .2)+
  labs(x = "Urban Cover (%)", y = "Nit. Dep (kg N/sqkm/year)")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
urb_temp<- ggplot(endo_herb)+
  geom_point(aes(x = PercentUrban, y = tmean_10km), alpha = .2)+
  labs(x = "Urban Cover (%)", y = "Mean Annual Temp.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
urb_ppt <- ggplot(endo_herb)+
  geom_point(aes(x = PercentUrban, y = ppt_10km), alpha = .2)+
  labs(x = "Urban Cover (%)", y = "Annual Precip.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))



nit_temp <- ggplot(endo_herb)+
  geom_point(aes(x = mean_TIN_10km, y = tmean_10km), alpha = .2)+
  labs(x =  "Nit. Dep (kg N/sqkm/year)", y = "Mean Annual Temp.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
nit_ppt <- ggplot(endo_herb)+
  geom_point(aes(x = mean_TIN_10km, y = ppt_10km), alpha = .2)+
  labs(x =  "Nit. Dep (kg N/sqkm/year)", y = "Annual Precip.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))


tmean_ppt <- ggplot(endo_herb)+
  geom_point(aes(x = tmean_10km, y = ppt_10km), alpha = .2)+
  labs(x =  "Mean Annual Temp.", y = "Annual Precip.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))



row_1 <- wrap_elements(ag_urb |ag_nit | ag_temp | ag_ppt )

row_2 <- wrap_elements(urb_nit | urb_temp | urb_ppt| plot_spacer() )
row_3 <-wrap_elements(nit_temp | nit_ppt| plot_spacer() |plot_spacer() )
row_4 <- wrap_elements(tmean_ppt |plot_spacer() |plot_spacer() |plot_spacer())
pairs_plot <- row_1 / row_2 / row_3 /row_4 + plot_annotation(tag_levels = "A")
ggsave(pairs_plot, filename = "Plots/pairs_plot.png", width = 8, height = 10)




# now plotting covariates against collection year

year_temp <- ggplot(endo_herb)+
  geom_point(aes(x = year, y = tmean_10km), alpha = .2)+
  geom_smooth(aes(x = year, y = tmean_10km))+
  labs(x =  "Collection Year", y = "Mean Annual Temp.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
year_ppt <- ggplot(endo_herb)+
  geom_point(aes(x = year, y = ppt_10km), alpha = .2)+
  geom_smooth(aes(x = year, y = ppt_10km))+
  labs(x =  "Collection Year", y = "Annual Precip.")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))
year_nit <- ggplot(endo_herb)+
  geom_point(aes(x = year, y = mean_TIN_10km), alpha = .2)+
  labs(x =  "Collection Year", y = "Nit. Dep (kg N/sqkm/year)")+
  theme_classic() + theme(panel.grid.major = element_line(color = "gray"))

year_plot <- year_temp + year_ppt + year_nit
ggsave(year_plot, filename = "Plots/year_pairs_plot.png", width = 8, height = 4)


