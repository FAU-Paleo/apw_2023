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

## First, let's pear down or occurrence data to only keep the info we need for making the map
locality_info_paleo <- occ_data %>% 
  select(collection_name, lat, lng, early_interval, late_interval, max_ma, min_ma) %>% 
  distinct(collection_name, .keep_all = TRUE) %>% 
  na.omit()



# Paleocoordinates --------------------------------------------------------

## To accurately reconstruct paleomaps, we should make sure our paleocoordinates from our
## occurrence data match the same scheme as the paleogeographies we want to use on our map 
## This can be a little tedious to begin with, but once you have your time intervals 
## standardised, the reconstruct() function in rgplates can reconstruct paleocordinates 
## from the modern coordinates supplied by the PBDB. 

## Let's re-do our Late Triassic - Early Jurassic maps, but with reconstructed palaeocoordinates

## First, assign an interval to every collection 
## Make an empty column:
locality_info_paleo$interval <- NA

for(i in 1:nrow(intervals)){
  ## Filter the occurrence data to contain localities from this interval only:
	
  # which collections correspond to this interval
	intervalIndex <- which(
		locality_info_paleo$max_ma <= intervals$max_ma[i] & 
		  locality_info_paleo$min_ma >= intervals$min_ma[i]
	)
	# which interval does this belong to 
	locality_info_paleo$interval[intervalIndex] <- intervals$interval_name[i]
}

## There are gaps here due to the mismatch of the timescales
## We have to rely on the primary input for binning!
locality_info_paleo$interval[c(11,76:80, 391, 408, 409,414, 417 )] <- "Hettangian"
locality_info_paleo$interval[c(59)] <- "Carnian"
locality_info_paleo$interval[c(62)] <- "Norian"
locality_info_paleo$interval[c(65, 123, 124,125,349:358, 367 , 369, 371, 395, 408, 431, 450:452)] <- "Toarcian"


## Set ggplot theme
palaeomap_theme <- theme_minimal() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
                                           axis.title.y=element_blank(), axis.text.y=element_blank(),
                                           axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
                                           legend.title=element_blank())
## Create an empty list to populate through the loop below:
palaeomap_list2 <- list()

for(i in 1:nrow(intervals)){
  
  ## Find the correct palaeomap based on the interval midpoint and reconstruct the plates:
  this_map <- reconstruct("plates", age = round(intervals$mid_ma[i]))
  ## Filter the occurrence data to contain localities from this interval only:
  these_occs <- locality_info_paleo[which(locality_info_paleo$interval == intervals$interval_name[i]), ]
  ## coordinates reconstruction
  new_coords <- reconstruct(these_occs[,c("lng", "lat")], age= round(intervals$mid_ma[i]))
  colnames(new_coords) <- c("plng", "plat") # rename columns
  ## column names: paleolong paleolat
  these_occs<- cbind(these_occs, new_coords)
  
  ## Set up your ggplot layers:
  palaeomap_list[[i]] <-  ggplot() +
    ## Add the landmasses for each interval:
    geom_polygon(data = this_map, aes(x = long, y = lat, group = group), fill = "grey75", color = "grey75") +
    ## Use the Mollweide projection:
    coord_map("mollweide") +
    ## Add the occurrence data (and set your colour!):
    geom_point(data = these_occs, aes(x = plng, y = plat), color = "#0DA69B", size = 4,  alpha = 0.8) + 
    ## Add lines from the x and y axes
    scale_y_continuous(breaks = seq(from = -90, to = 90, by = 30), limits = c(-90,90)) + 
    scale_x_continuous(breaks = seq(from = -180, to = 180, by = 30), limits = c(-180,180)) + 
    ## Add the interval name to the title of each map 
    ggtitle(paste(" ", intervals$interval_name[i], sep = "")) +
    ## Finally, add the custom theme
    palaeomap_theme

	cat(i, "\r")
	flush.console()
}

## To make a panel plot containing all the map, we'll need to arrange them first:
## Here, there are 7 intervals, so 7 map plots:
palaeomaps_panel2 <- ggarrange(palaeomap_list2[[1]], palaeomap_list2[[2]],
                               palaeomap_list2[[3]], palaeomap_list2[[4]],
                               palaeomap_list2[[5]], palaeomap_list2[[6]],
                               palaeomap_list2[[7]],
                               ncol = 2, nrow = 4) # no. of rows and columns

palaeomaps_panel2 # Call the panel to the plot window to have a peep

## And finally, save as a .pdf
ggsave(plot = palaeomaps_panel2,
       width = 12, height = 10, dpi = 600, 
       filename = "./plots/Palaeomaps_panel2.pdf", useDingbats=FALSE)

