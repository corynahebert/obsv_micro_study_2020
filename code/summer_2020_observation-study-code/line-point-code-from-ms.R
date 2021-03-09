### How are community composition and microtopography related?
# just trying a PCA
# points (rows) will be sites
# loadings (columns) will have all species and topo variables

setwd("C:/Users/Coryna/Documents/Wetlan_USU/2020_observational_study/data/raw")

install.packages("data.table")
install.packages("visreg")

library(data.table)
library(ggplot2)
library(ggfortify)
library(visreg)
library(RColorBrewer)

plant_data <- fread("line_point_data.csv", stringsAsFactors = FALSE)
topo_data <- fread("micro_measurements.csv", stringsAsFactors = FALSE)

# calculating topo variables
topo_data[, relief := max_pin_height_cm-min_pin_height_cm]
topo_data[, roughness := pin_profile_length_cm/profilometer_width_cm]

hist(topo_data$relief) # you probably want to invesitgate that negative relief value
hist(topo_data$roughness) # and those roughness values below 1

# cleaning plant data
plant_data[top_layer %in% c("NONE", "NOE", "NOINE", "NON"), "top_layer"] <- NA
plant_data[top_layer %in% c("DISIP", "DISPNONE", "DSIP"), "top_layer"] <- "DISP"

# site and species lists
sites <- sort(unique(plant_data$site_id))
species <- sort(unique(plant_data$top_layer))

plant_comm <- data.table(site = sites) # initializing up a master table

# summarizing species occurence data into a tidy "plant_comm" table with one column per species
for(sp in species){
  #sp <- species[1]
  sp_records <- plant_data[top_layer == sp, ] # pull out data for just one species - I used only the top_layer data
  sp_props <- table(sp_records$site_id)/150 # get relative proportion... there's 150 total plant measurements per site, right?
  site_props <- data.table(site = sites,
                           sp_prop = 0) # just setting up a temporary site & species proportion data table to keep things organized
  site_props[which(site_props$site %in% names(sp_props)), "sp_prop"] <- sp_props # inserting the species proportions
  
  plant_comm[,sp] <- site_props$sp_prop # adding in a new column into the "plant_comm" table with the calculated proportions
}

# function to pull out relief and roughness data by site
get.topo <- function(site, var){
  #site <- sites[1]
  if(var == "relief") relevant_data <- topo_data[site_id == site, relief]
  if(var == "roughness") relevant_data <- topo_data[site_id == site, roughness]
  return(median(relevant_data, na.rm=T)) # median so that single outliers don't derail average
}

# get median relief and roughness at each site
plant_comm$relief <- mapply(get.topo, site=sites, var="relief")
plant_comm$roughness <- mapply(get.topo, site=sites, var="roughness")

# getting rid of sites with NAs in topo variables
plant_comm <- na.omit(plant_comm)

# does it look about right?
plant_comm


# PCA with all species percentages and topography
pca_all <- prcomp(plant_comm[, -1])
autoplot(pca_all,  loadings=T, loadings.label=T,
         colour = "relief") # relief explains most of the variation, but nothing else strongly covaries with it
summary(pca_all)

# PCA with just species data
pca_community <- prcomp(plant_comm[, -c("site", "relief", "roughness")])
autoplot(pca_community,   loadings=T, loadings.label=T) # There's one weird site with lots of SCAM, most of the rest of the variation is explained by DISP
summary(pca_all)

### GGplots ------

relief_hist <- ggplot(data = topo_data, aes(x = relief, fill = ..x..))+
  geom_histogram(binwidth = 1)+
  xlab("Relief (cm)")+
  ylab("Count")+
  scale_fill_gradient(low = "blue", high = "green")
  

relief_hist_2 <- ggplot(data = topo_data, aes(relief)) +
  geom_bar(color="blue", fill="blue") +
  scale_x_binned() +
  xlab("Maximum Relief (cm)") +
  ylab("Count")


relief_by_site <- ggplot(data = plant_comm, aes(site, relief))+
  geom_point()

roughness_by_site <- ggplot(data = plant_comm, aes(site, roughness))+
  geom_point()

# Sequential color scheme
hp+scale_fill_gradient(low="blue", high="red")

roughness_hist <- ggplot(data = topo_data, aes(roughness))+ 
  geom_bar(color = "orange", fill = "orange")+
  scale_x_binned()+
  xlab("Roughness (unitless)")+
  ylab("Count")
  

head(topo_data)

## Linear models ----------

### much simpler quick analyses looking at if the proportion of 
### phragmites or DISP (the two most abundant species, it looks like) 
### predicts relief or roughness

relief_model <- lm(relief ~ PHRAG + DISP + SCAM + SCAC, data=plant_comm)
par(mfrow=c(2,2))
visreg(relief_model, ylab= "Relief (cm)", 
       xlab = "Phragmites", "Saltgrass", "Threesquare bulrush", 
       "Hardstem bulrush")
summary(relief_model) # looks like there might be relationships, but they're not significant at least when modeled like this
# this explains 19% of variation

plant_comm

?visreg

roughness_model <- lm(roughness ~ PHRAG + DISP + SCAM + SCAC, data=plant_comm)
par(mfrow=c(2,2))
visreg(roughness_model,
       xlab = "Pragmites" + "Saltgrass" + "Threesquare bulrush" + "Hardstem bulrush",
       ylab = "Roughness")
summary(roughness_model) # DISP actually significantly predicts roughness! (at alpha=0.05)
# this explains 26% of the variation

## NMDS attempt ---------

#install packages

install.packages("Car")
install.packages("emmeans")
install.packages("DHARMa")
install.packages("MuMIn")
install.packages("magrittr")
install.packages("xts")
install.packages("")


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

