### How are community composition and microtopography related?
# just trying a PCA
# points (rows) will be sites
# loadings (columns) will have all species and topo variables

setwd("C:/Users/Coryna/Documents/Wetlan_USU/2020_observational_study/data/raw")


library(data.table)
library(ggplot2)
library(ggfortify)
library(visreg)
library(RColorBrewer)
library(tidyverse)
library(vegan)
library(permute)
library(lattice)

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

## make another matix where we drop rare species 

plant_comm2 <- plant_comm %>% 
  group_by(site) 

## NMDS of plant_comm -------

plotdatamatrix<- as.matrix(plant_comm[, -c("site", "relief", "roughness")])

env <- data.frame(plant_comm$relief, plant_comm$roughness)

plotdatamatrix2<-as.matrix(plant_comm[, -c("site", "relief", "roughness", 
                                           "AP Z", "BOME", "HOMA", "LASE", 
                                           "MEAL", "POLYGONUM", "PONO",
                                           "SCPU", "SPMA", "SUCA", "TRFR", 
                                           "VEAN AQ", "UNK BROME", "UNK FORB",
                                           "UNK GRASS")])


#NMDS

obsv_NMDS = metaMDS(plotdatamatrix, k=3) 


## stress is minimized with 3 dimensions

en = envfit(obsv_NMDS, env, permutations = 999, na.rm = TRUE)


obsv_NMDS

plot(obsv_NMDS)
plot(en)
stressplot(obsv_NMDS)

#### NMDS dropped rare species 

NMDS_2 = metaMDS(plotdatamatrix2, k=2)
en = envfit(NMDS_2, env, permutations = 999, na.rm = TRUE)
plot(NMDS_2)
plot(en)
stressplot(NMDS_2)

#extract scores and make sure side_id columns added
data.scores <- as.data.frame(scores(NMDS_2))#Using the scores function to extract the site scores and convert to a data.frame
data.scores$site_id <- plant_comm$site  
data.scores$roughness <- plant_comm$roughness
data.scores$relief <- plant_comm$relief


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#species scores
species.scores <- as.data.frame(scores(NMDS_2, "species"))  
#Using the scores function from vegan to extract the species scores and convert\
#to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores


## plotting ------ 

gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = site_id), size = 5, alpha = 0.5) +
  geom_point(data = species.scores, size = 1, alpha = 0.5)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 2, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) +
  geom_text(data = species.scores, aes(label=species), hjust=0, vjust=0) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Site")

gg2 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = site_id), size = 5, alpha = 0.5) +
  geom_point(data = species.scores, size = 1, alpha = 0.5)+
  geom_segment(aes(x = NMDS1, y = NMDS2), 
               data = en_coord_cont, size =1, alpha = 2, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) +
  geom_text(data = species.scores, aes(label=species), hjust=0, vjust=0) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Site")

?geom_point

NMDS_Data <- data.frame(MDS1=obsv_NMDS$points[,1], MDS2=obsv_NMDS$points[,2], group=data.scores$site_id)
ord <- ordiellipse(obsv_NMDS, data.scores$site_id, display = "sites", kind= "se", conf = 0.95, label = T, draw="lines", col="red") 

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

nmds_plot <- ggplot(data = data.scores, aes(NMDS1, NMDS2)) + 
  geom_point(aes(color = site_id))

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

## NMDS attempt


# load packages
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

## calculating relief and roughness 
## average values and create separate dataframe

topo_data <- micro_data %>% mutate(relief = max_pin_height_cm-min_pin_height_cm,
                roughness = pin_profile_length_cm/profilometer_width_cm)

rough_relief <- subset(topo_data, select=c("site_id", "roughness", "relief"))

topo_data_stats <- rough_relief %>% group_by(site_id)
                

                

glimpse(topo_data)
view(rough_relief)




## Need to create a site by species matrix to perform NMDS 
