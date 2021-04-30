## NMDS attempt ---------

#install packages

install.packages("Car")
install.packages("emmeans")
install.packages("DHARMa")
install.packages("MuMIn")
install.packages("magrittr")
install.packages("xts")
install.packages("plyr")
install.packages("BiodiversityR")
install.packages("MASS")
install.packages("labdsv")
install.packages("grid")
install.packages("ggrepel")


# load packages ----
library(glmmTMB)
library(car)
library(emmeans)
library(DHARMa)
library(MuMIn)
library(tidyverse)
library(magrittr)
library(xts)
library(plyr)
library(dplyr)
library(vegan)
library(BiodiversityR) 
library(MASS)
library(ecodist)
library(labdsv)
library(ggplot2)
library(grid)
library(ggrepel)

#I do not want R to use exponential notation
options(scipen=999)

## Importing raw data, the 'fileencoding argument gets ride of weird chartacters 
## in my column names 

plant_data <- read.csv("line_point_data.csv",
                       fileEncoding="UTF-8-BOM",
                       stringsAsFactors = FALSE)
micro_data <- read.csv("micro_measurements.csv",
                       fileEncoding="UTF-8-BOM",
                       stringsAsFactors = FALSE)

glimpse(plant_data)

## calculating relief and roughness ------
## average values and create separate dataframe

topo_data <- micro_data %>% mutate(relief = max_pin_height_cm-min_pin_height_cm,
                                   roughness = pin_profile_length_cm/profilometer_width_cm)

rough_relief <- subset(topo_data, select=c("site_id", "roughness", "relief"))

topo_stats <- rough_relief %>%
  group_by(site_id) %>%
  dplyr::summarise(mean_rough = mean(roughness, na.rm = TRUE),
                   SD_rough = sd(roughness, na.rm = TRUE),
                   nsize_rough = sum(!is.na(roughness)),
                   SE_rough = SD_rough/sqrt(nsize_rough),
                   mean_relief = mean(relief, na.rm = TRUE),
                   SD_relief = sd(relief, na.rm = TRUE),
                   nsize_relief = sum(!is.na(relief)),
                   SE_relief = SD_relief/sqrt(nsize_relief))



glimpse(topo_data)
view(rough_relief)


## GGplots of relief and roughness by site ----

reief_by_site <- ggplot(topo_stats, aes(site_id, mean_relief)) +
  geom_bar(stat = "identity")

roughness_by_site <- ggplot(topo_stats, aes(site_id, mean_rough)) +
  geom_bar(stat= "identity")



## Need to create a site by species matrix to perform NMDS, I have to transform 
## my sites to percent cover by species 

head(plant_data)

species <- unique(plant_data$top_layer)
            
site_by_species <- mutate()

?mutate
