setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/4_GIS/")

library(sf) # actual spatial package
library(viridis) # for using a color scheme
library(divDyn) # for example occurrence data from the Paleobiology Database
library(chronosphere) # for map edge and shaper()
library(geosphere) # great circles


# basic plotting
# Natural Earth
ne <- sf::st_read("data/ne_land/ne_110m_land.shx")

# show
ne

# structure
str(ne)

# the geometry: where the things are
ne$geometry

#plotting just the available shape
plot(ne$geometry)

# regular plotting parameters
plot(ne$geometry, col="green", border="black", xlim=c(-45, 45), ylim=c(45,45))

# single points
erlangen <- c(10.97, 49.58)
sydney <- c(151.19, -33.86)
both <- rbind(erlangen, sydney)


plot(ne$geometry)
points(both, col="red", pch=3)

# great circle
gc <- gcIntermediate(erlangen, sydney)
points(gc, col="blue", pch=3)

# Haversine distance
distHaversine(erlangen, sydney)


# the fossil occurrences
data(corals)
points(corals$lng, corals$lat, col="red", pch=3)

# making spatial points
collections <- unique(corals[, c("collection_no", "lng", "lat")])

# the number of occurrences
nOccs <- table(corals$collection_no)
collections$nOccs <- nOccs[as.character(collections$collection_no)]

# omit missing values
collections <- collections[!is.na(collections$lng) & !is.na(collections$lat), ]

# as simple feature!
sfColls <- st_as_sf(collections, coords=c("lng", "lat"))
plot(sfColls["nOccs"],pch=16)

# Change axis to log
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5)

# preparing multilayer plots
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, reset=FALSE)
plot(ne$geometry, add=TRUE)

# using it a background
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, reset=FALSE, main="Number of occurrences")
plot(ne$geometry, add=TRUE, border=NA, col="gray")
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, add=TRUE)


# projection change
# The mollweide projection
st_crs(ne)

mollProj <- "ESRI:54009"
# same in mollweide
neMoll <- st_transform(ne, crs=mollProj)

# won't work yet
coralsMoll <- st_transform(sfColls, crs=mollProj)

# we need to assign a crs to transform from
st_crs(sfColls) <- st_crs(ne)
coralsMoll <- st_transform(sfColls, crs=mollProj)

# the same thing as before in Mollweide
plot(coralsMoll["nOccs"], pch=16, reset=FALSE, logz=TRUE)
plot(neMoll$geometry, col="gray", border=NA, add=TRUE)
plot(coralsMoll["nOccs"], pch=16, add=TRUE, logz=TRUE)


# the boundary of the map? 
# This calculates it as a "Spatial*"
edge <- mapedge()
# use this to transform sp-type objects to sf-type objects
edge <- as(edge, "sf")
edgeMoll <- st_transform(edge, crs=mollProj)

# the same with the background map
plot(coralsMoll["nOccs"], pch=16, reset=FALSE, logz=TRUE, main="The number of occurrences in a collection")
plot(edgeMoll$geometry, col="#9bbff4", add=TRUE)
plot(neMoll$geometry, col="white", border=NA, add=TRUE)
plot(coralsMoll["nOccs"], pch=16, add=TRUE, logz=TRUE)

################################################################################
# Countries
countries <- sf::st_read("data/world_admin/world-administrative-boundaries.shx")
plot(countries$geometry)

# all available attributesg
plot(countries)

# one specific attribute
plot(countries["region"])

# overlay
joined <- st_join(countries, sfColls, join=st_intersects)

# the number of points
nPoints <- table(joined$iso3)
countries$nColl <- nPoints[countries$iso3]

# plot a country map based on the number of collections in it:
plot(countries["nColl"], logz=TRUE)

# manual ramping
breaks <- c(0, 5, 10, 50, 100, 500, 1000 )
colors <- viridis(6)
plot(countries["nColl"], breaks=breaks, pal=colors)


################################################################################
# IUCN range data
chameleons <- sf::st_read("data/CHAMELEONS/CHAMELEONS.shx")

plot(chameleons$geometry, col="#00808066", border=NA, reset=FALSE)
plot(ne$geometry, add=TRUE, border="black")

# plotting a single species
# subset as a data.frame!
one <- chameleons[which(chameleons$binomial=="Furcifer rhinoceratus"), ]

# preliminary plot - use shaper()
plot(ne$geometry, reset=FALSE, col="gray", xlim=c(34.9, 56.85), ylim=c(-31.92, -8.13))
plot(one$geometry, add=TRUE, col="blue")

# The area of this shape
st_area(one)

# the area measured
area <- rep(NA, nrow(chameleons))
for(i in 1:nrow(chameleons)){
	# basic exception handling in R
	try(area[i]<- st_area(chameleons[i, ]))
	cat(i, "\r")
	flush.console()
}
chameleons$area <- area

# Actual area
plot(chameleons$SHAPE_Area, chameleons$area)

# The total area
cleaned <-st_make_valid(chameleons$geometry, 1e7)

# the 186th geometry has a bit of problems!
union <- st_union(cleaned[-186])

plot(union, col="blue")
st_area(st_make_valid(union))

# which are tropical?
# Define a polygon: 
poly <- matrix(c(
		 -180, 23.5,
		180, 23.5,
		180, -23.5,
		-180, -23.5, 
		-180, 23.5	
	),
	ncol=2,
	byrow=TRUE)
