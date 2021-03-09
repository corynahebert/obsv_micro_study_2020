## Line Point Data for observational study 

setwd("C:/Users/Coryna/Documents/Wetlan_USU/2020_observational_study/data/raw")

## import raw data and look at srtucture 

line_point <- read.csv("line_point_data.csv")

str(line_point)

micro <- read.csv("micro_measurements.csv")

str(micro)

## load packages 

library(tidyverse)
library("RColorBrewer")

## calculate roughness and relief 

micro_measurements <- micro %>%
  dplyr::summarise(roughness = (pin_profile_length_cm/profilometer_width_cm),
                   relief = (max_pin_height_cm-min_pin_height_cm))
micro_measurements

## put them in the same table?? 
