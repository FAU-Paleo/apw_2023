# *************************************************
#
#   Analytical Paleobiology Workshop 2022
#
#   Module 2: Paleodiversity analyses
#   Day 4: Exploring fossil occurrence data  
#             and fossil record biases
#
#   Thursday, August 25th 2022
#   Emma Dunne (emma.dunne@fau.de)
# _________________________________________________
#
#   04: Modern and palaeo- maps
# 
# *************************************************


## Load package(s):
library(tidyverse) # for data manipulation functions and ggplot
library(rgplates) # for palaeomaps (and other cool stuff)
library(ncdf4) # helps rgplates and chronosphere
library(ggpubr) # for arranging plots into panels



# Modern world map --------------------------------------------------------

## First, let's pear down or occurrence data to only keep the info we need for making the map
locality_info <- occ_data %>% 
  select(collection_name, lat, lng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit()


## Grab a world map for ggplot to work with:
world_map <- map_data("world")
ggplot() + geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region)) 

## Let's make it pretty and add our data
modern_map <- ggplot() + 
  geom_map(data = world_map, map = world_map, aes(long, lat, map_id = region), 
           color = "grey80", fill = "grey90", size = 0.1) +
  geom_point(data = locality_info, aes(lng, lat), alpha = 0.7, colour = "#23AAB8") +
  theme_void() + theme(legend.position = "none")
modern_map

## And save as a .pdf
ggsave(plot = modern_map,
       width = 8, height = 5, dpi = 600, 
       filename = "./plots/Modern_map.pdf", useDingbats=FALSE)




# Palaeogeographic map ---------------------------------------------------

## In this section, we're going to plot the PBDB collections onto palaeogeographic maps
## Side note: You should always cite the R packages you use!
citation("rgplates")

## First, let's explore the rgplates package a little

## Fetch the PaleoDEMs (Scotese and Wright, 2018) - this might take a while to load
dems <- chronosphere::fetch(dat = "paleomap", var="dem")
dems # returns a RasterArray class object

## Example plot for the ~late Silurian (based on age in Ma)
plot(dems["210"], col = gradinv(210)) 
showPal() # Color palettes usually used for plotting environmental data


## Now, let's pear down or occurrence data to only keep the info we need for making the map
locality_info_paleo <- occ_data %>% 
  select(collection_name, paleolat, paleolng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit()

## Begin by setting a theme (I like to keep things pretty minimal):
palaeomap_theme <- theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                                           axis.title.y=element_blank(), axis.text.y=element_blank(),
                                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                           legend.title=element_blank())

## Now let's plot the each of the maps!
#£ First, create an empty list to populate through the loop below:
palaeomap_list <- list()

## Now, let's loop through each interval and plot a map for it
## This might take a minute or two to run as its doing a lot of work!
for(i in 1:nrow(intervals)){
  
  ## First, find the correct palaeomap based on the interval midpoint and reconstruct the landmasses:
  this_map <- reconstruct("plates", age = round(intervals$mid_ma[i])) 
  ## Filter the occurrence data to contain localities from this interval only:
  these_occs <- filter(locality_info_paleo, max_ma <= intervals$max_ma[i] & min_ma >= intervals$min_ma[i])
  
  ## Set up your ggplot layers:
  palaeomap_list[[i]] <-  ggplot() +
    ## Add the landmasses for each interval:
    geom_polygon(data = this_map, aes(x = long, y = lat, group = group), fill = "grey75", color = "grey75") +
    ## Use the Mollweide projection:
    coord_map("mollweide") +
    ## Add the occurrence data (and set your colour!):
    geom_point(data = these_occs, aes(x = paleolng, y = paleolat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
    ## Add lines from the x and y axes
    scale_y_continuous(breaks = seq(from = -90, to = 90, by = 30), limits = c(-90,90)) + 
    scale_x_continuous(breaks = seq(from = -180, to = 180, by = 30), limits = c(-180,180)) + 
    ## Add the interval name to the title of each map 
    ggtitle(paste(" ", intervals$interval_name[i], sep = "")) +
    ## Finally, add the custom theme
    palaeomap_theme
}

## To make a panel plot containing all the map, we'll need to arrange them first:
## Here, there are 7 intervals, so 7 map plots:
palaeomaps_panel <- ggarrange(palaeomap_list[[1]], palaeomap_list[[2]],
                              palaeomap_list[[3]], palaeomap_list[[4]],
                              palaeomap_list[[5]], palaeomap_list[[6]],
                              palaeomap_list[[7]],
                              ncol = 2, nrow = 4) # no. of rows and columns

palaeomaps_panel # Call the panel to the plot window to have a peep

## And finally, save as a .pdf
ggsave(plot = palaeomaps_panel,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomaps_panel.pdf", useDingbats=FALSE)

