# Basic examples with rgplates 
# Erlangen, 2023
# Ádám T. Kocsis
# CC-BY (attribution)
# see more examples at : https://adamkocsis.github.io/rgplates/
setwd("3_rgplates")

library(rgplates)
# might need to install
# install.packages("geojsonsf") # for interacting with the GWS
# install.packages("sf") # general dependency 

################################################################################
# 1. The GPlates web service - requires internet connection
# https://gwsdoc.gplates.org/

# reconstructable features:
# - static_polygons: non-deformed polygons of plates
# - plate_polygons
# - coastlines: present-day location of coastlines 

################################################################################
# A. Polygons and coastlines

# The default model: Merdith et al. 2021
# https://gwsdoc.gplates.org/reconstruction-models#merdith2021
# present day
map0 <- reconstruct("static_polygons", age=0)

# plot via geometry
plot(map0$geometry)

# coastlines on top
coast0 <- reconstruct("coastlines", age=0)
plot(coast0$geometry, col="gray", add=TRUE, border="blue")


# Past reconstruction
map150 <- reconstruct("static_polygons", age=150)
coast150 <- reconstruct("coastlines", age=150)

# visualize
plot(map150$geometry)
plot(coast150$geometry, col="gray", add=TRUE, border="blue")


# making a nice plot
proj <- "ESRI:54009"
coast150Moll <- sf::st_transform(coast150, proj)
plot(coast150Moll$geometry)

# making reall ynice
edge <- mapedge()
edgeMoll <- sf::st_transform(edge, proj)
plot(edgeMoll, col="#0022AA")
plot(coast150Moll$geometry, col="white", border=NA, add=TRUE)


# contrast with another model (mantle reference frame)
# https://gwsdoc.gplates.org/reconstruction-models#muller2022
coast150MR <- reconstruct("coastlines", age=150, model="MULLER2022")
coast150MRMoll <- sf::st_transform(coast150MR, proj)

plot(coast150Moll$geometry, col="#DD000055", border=NA, add=TRUE)

# Deeper in time!
coast400 <- reconstruct("coastlines", age=400)
coast400MR <- reconstruct("coastlines", age=400, model="MULLER2022")
# transform to Mollweide
coast400Moll <- sf::st_transform(coast400, proj)
coast400MRMoll <- sf::st_transform(coast400MR, proj)

# and plot
plot(edgeMoll, col="#0022AA")
plot(coast400Moll$geometry, col="white", border=NA, add=TRUE)
plot(coast400MRMoll$geometry, col="#DD0000AA", border=NA, add=TRUE)

################################################################################
# B. Points
erlangen <- c(10.97, 49.58)
sydney <- c(151.19, -33.86)
both <- rbind(erlangen, sydney)
colnames(both) <- c("lng", "lat")

# reconstructing past positions, default model
both150 <- reconstruct(both, age=150)

# plotting
plot(edge, col="#0022AA")
plot(coast150$geometry, col="white", border=NA, add=TRUE)
points(both150, cex=2, pch=3, col="red")

# one Mollweide - project the points too!!
both150_sf <- sf::st_as_sf(as.data.frame(both150), coords=c("paleolong", "paleolat"))
# add CRS
sf::st_crs(both150_sf) <- "WGS84"
# finally project
both150Moll <- sf::st_transform(both150_sf, proj)

# and plot
plot(edgeMoll, col="#0022AA")
plot(coast150Moll$geometry, col="white", border=NA, add=TRUE)
# this is sf! needs plot() not points()!
plot(both150Moll$geometry, cex=2, pch=3, col="red", add=TRUE)

################################################################################

# compare Erlangen's past
ages <- seq(500,0, -50)
bothSeries <- reconstruct(both, age=ages)
bothSeriesMR <- reconstruct(both, age=ages, model="MULLER2022")
bothSeriesPM <- reconstruct(both, age=ages, model="PALEOMAP")

# extract Erlangen's latitude
erLat<- unlist(lapply(bothSeries, function(x) x["erlangen", "lat"]))
erLatMR<- unlist(lapply(bothSeriesMR, function(x) x["erlangen", "lat"]))
erLatPM <- unlist(lapply(bothSeriesPM, function(x) x["erlangen", "lat"]))

library(divDyn)
data(stages)
tsplot(stages, boxes="sys",ylim=c(-90,90))
abline(h=0, lty=2)
lines(ages, erLat, type="o", pch=16)
lines(ages, erLatMR, type="o", col="red", pch=16)
lines(ages, erLatPM, type="o", col="blue", pch=16)
legend("topleft", legend=c("MERDITH2021", "MULLER2022", "PALEOMAP"), col=c("black","red", "blue"), lwd=1, pch=16)


################################################################################
# 2. Offline reconstruction with the GPlates desktop application
# This only works with the GPlates Desktop application installed
# standard files: -> ?reconstruct

# download files from the chronosphere
library(chronosphere)

# The Merdith et al model
merdith2021<- fetch("EarthByte", "MERDITH2021")

# Why? 1. More elaborate, persistent model
# see reconstructable features (not all work in rgplates, yet!)
polygons65 <- reconstruct("static_polygons", age=65, model=merdith2021)
continents65 <- reconstruct("continents", age=65, model=merdith2021)
coastlines65 <- reconstruct("coastlines", age=65, model=merdith2021)
cratons65 <- reconstruct("cratons", age=65, model=merdith2021)

plot(polygons65$geometry, col="gray80")
plot(continents65$geometry, col="gray40", add=TRUE)
plot(cratons65$geometry, col="black",  add=TRUE)
plot(coastlines65$geometry, col=NA, border="red", add=TRUE)

# Why? 2. more data processing
# -> 4_paleorasters.R



