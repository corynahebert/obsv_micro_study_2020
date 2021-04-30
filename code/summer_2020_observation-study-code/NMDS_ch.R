### How are community composition and microtopography related?


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

hist(topo_data$relief) 
hist(topo_data$roughness) 

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
  sp_props <- table(sp_records$site_id)/150 # get relative proportion, there's 150 total plant measurements per site
  site_props <- data.table(site = sites,
                           sp_prop = 0) # just setting up a temporary site & species proportion data table 
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

plant_comm




## NMDS of plant_comm -------

plotdatamatrix<- as.matrix(plant_comm[, -c("site", "relief", "roughness")])

env <- data.frame(micro_measurements$relief, plant_comm$roughness)

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

