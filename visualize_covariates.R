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
