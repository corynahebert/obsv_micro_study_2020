## Line Point Data for observational study 

## load packages 

library(tidyverse)
library("RColorBrewer")
library("vegan")
library("permute")
library("lattice")
library(viridisLite)
library(viridis)
library(patchwork)


## import raw data using tidyverse function and look at structure 

line_point <- read_csv("line_point_data.csv")

micro <- read_csv("micro_measurements.csv")
str(micro)


## calculate roughness and relief and summary stats ---------

micro_measurements <- micro %>%
  mutate(roughness = (pin_profile_length_cm/profilometer_width_cm),
         relief = (max_pin_height_cm-min_pin_height_cm)) %>% 
  group_by(site_id) %>% 
  summarise(mean_roughness = mean(roughness, na.rm = TRUE),
            sd_roughness = sd(roughness, na.rm = TRUE),
            nsize_rough = sum(!is.na(roughness)),
            se_roughness = sd_roughness/sqrt(nsize_rough),
            mean_relief = mean(relief, na.rm = TRUE),
            sd_relief = sd(relief, na.rm = TRUE), 
            nsize_relief = sum(!is.na(relief)), 
            se_relief = sd_relief/sqrt(nsize_relief))

##ggplots of microtopography by site --------

roughness_by_site <- ggplot(micro_measurements, 
                         aes(x = site_id, y = mean_roughness, 
                             fill = site_id)) +
  geom_col() + 
  scale_fill_viridis(option = "plasma", discrete = TRUE) +
  labs(x = "Site", 
       y = "Mean roughness (unitless)") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 75)) 

relief_by_site <- ggplot(micro_measurements, 
                         aes(x = site_id, y = mean_relief, 
                             fill = site_id)) + 
  geom_col() +
  scale_fill_viridis(option = "plasma", discrete = TRUE) +
  labs(x = "Site", 
       y = "Mean relief (cm)") + 
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 75)) 

##patchwork them together for comparison 
roughness_by_site | relief_by_site

## Change line point data into percent cover of species by site
## Then this data can be used to conduct an NMDS

tl <- line_point %>% select(c(2,6)) %>% rename(species = top_layer)
c1 <- line_point %>% select(c(2,7)) %>% rename(species = code_1)
c2 <- line_point %>% select(c(2,8)) %>% rename(species = code_2)
c3 <- line_point %>% select(c(2,9)) %>% rename(species = code_3)

new_df <- rbind(tl,c1,c2,c3)

percent_cover <- new_df %>% 
  group_by(site_id) %>% 
  count(species) %>% 
  filter(!is.na(species)) %>% 
  mutate(n_hits = n,
         percent_cover = n/150) %>% 
  select(site_id,
         species,
         percent_cover) %>%
  pivot_wider(names_from = species,
              values_from = percent_cover) %>% 
  replace(is.na(.), 0) %>% 
  filter(site_id != "FABA-01") %>% 
  filter(site_id != "PSG-01") %>% 
  filter(site_id != "NPDC-01")




## NMDS ----

##convert percent cover data into matrix

df <- as.data.frame(percent_cover)

df2 <- subset(df, select = -site_id)

plotdatamatrix<- as.matrix(df2)

str(plotdatamatrix)

##create a dataframe of environmental variables 

env <- data.frame(micro_measurements$mean_roughness, 
                  micro_measurements$mean_relief)


#NMDS
obsv_NMDS = metaMDS(plotdatamatrix, k=3) ## stress is minimized with 3 dimensions

en = envfit(obsv_NMDS, env, permutations = 999, na.rm = TRUE)


obsv_NMDS

plot(obsv_NMDS)
plot(en)
stressplot(obsv_NMDS)


#extract scores and make sure side_id columns added
data.scores <- as.data.frame(scores(obsv_NMDS))

#Using the scores function to extract the site scores and convert to a data.frame

data.scores$site_id <- df$site_id  
data.scores$Roughness <- micro_measurements$mean_roughness
data.scores$Relief <- micro_measurements$mean_relief


en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

library(data.table)
(setattr(en_coord_cont, "row.names", c("Roughness", "Relief")))



#species scores
species.scores <- as.data.frame(scores(obsv_NMDS, "species"))  
#Using the scores function from vegan to extract the species scores and convert\
#to a data.frame
species.scores$species <- rownames(species.scores)  
# create a column of species, from the rownames of species.scores


## want to remove the uncommon species for plotting   
species.scores_common <- species.scores %>% 
  slice(1:3, 5, 6, 9, 11, 14, 18, 20:22) 

species.scores_common$species <- rownames(species.scores_common) 


## plotting ------ 

NMDS_1 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = site_id), size = 5, alpha = 0.5) +
  geom_point(data = species.scores_common, size = 0, alpha = 0.5)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size =1, alpha = 2, colour = "grey30") +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) +
  geom_text(data = species.scores_common, aes(label=species), hjust=0, vjust=0) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Site") + 
  theme_minimal()

ggsave("NMDS_1.png", width = 30, height = 25, units = "cm")





