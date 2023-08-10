# Basic vector data examples
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/4_GIS/")

# necessary packages
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

# plotting a png image
dir.create("export", showWarnings=FALSE)
png("export/worldmap.png", height=3000, width=2000)
plot(ne$geometry)
dev.off()

# regular plotting parameters
plot(ne$geometry, col="green", border="red", xlim=c(-45, 45), ylim=c(45,45))

# single points
erlangen <- c(10.97, 49.58)
sydney <- c(151.19, -33.86)
both <- rbind(erlangen, sydney)
colnames(both) <- c("lng", "lat")

# plotting them both
plot(ne$geometry)
points(both, col="red", pch=3)

# great circle
gc <- geosphere::gcIntermediate(erlangen, sydney, n=100)
points(gc, col="blue", pch=3)
lines(gc, col="blue", pch=3)


# dateline problem with basic plotting!
# locator(2)

# two points around the dateline
point1 <- c(135,48 )
point2 <- c(-118, 44.8)
datePoints <- rbind(point1, point2)
colnames(datePoints) <- c("lng", "lat")

# the wrong plotting! - be aware of the dateline! 
plot(ne$geometry)
points(datePoints, col="red", pch=3)
gc2<- geosphere::gcIntermediate(point1, point2, n=100)
lines(gc2, col="blue", pch=3)

# There should be a break at the breakline
gc2<- geosphere::gcIntermediate(point1, point2, n=100, breakAtDateLine 	
=TRUE)
plot(ne$geometry)
points(datePoints, col="red", pch=3)
lines(gc2[[1]], col="blue", pch=3)
lines(gc2[[2]], col="blue", pch=3)

# Haversine distance (GCD)
distHaversine(erlangen, sydney)

# the fossil occurrences from divDyn
data(corals)
plot(ne$geometry)
points(corals$lng, corals$lat, col="red", pch=3)

# making spatial points
collections <- unique(corals[, c("collection_no", "lng", "lat")])

# the number of occurrences
nOccs <- table(corals$collection_no)
collections$nOccs <- nOccs[as.character(collections$collection_no)]

# omit missing values
collections <- collections[!is.na(collections$lng) & !is.na(collections$lat), ]

# as simple feature!
sfColls <- sf::st_as_sf(collections, coords=c("lng", "lat"))

# plotting all attributes (collection_id is numeric - alhough meaningless)
plot(sfColls)

# Subsetting to one variable
plot(sfColls[,"nOccs"],pch=16)

# Change axis to log
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5)

# preparing multilayer plots - without reset there is misalignment!
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, reset=FALSE)
plot(ne$geometry, add=TRUE)

# just map as a background - no legend!
plot(ne$geometry, border=NA, col="gray")
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, add=TRUE)


# using it as a background, with legend
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, reset=FALSE, main="Number of occurrences")
plot(ne$geometry, add=TRUE, border=NA, col="gray")
plot(sfColls["nOccs"],pch=16, logz=TRUE, cex=0.5, add=TRUE)


# projection change
# The mollweide projection
sf::st_crs(ne)

mollProj <- "ESRI:54009"
# same in mollweide
neMoll <- st_transform(ne, crs=mollProj)
plot(neMoll)

# won't work yet
coralsMoll <- st_transform(sfColls, crs=mollProj)

# we need to assign a crs to transform from
st_crs(sfColls) <- st_crs(ne)
coralsMoll <- sf::st_transform(sfColls, crs=mollProj)

# the same thing as before in Mollweide
plot(coralsMoll["nOccs"], pch=16, reset=FALSE, logz=TRUE)
plot(neMoll$geometry, col="gray", border=NA, add=TRUE)
plot(coralsMoll["nOccs"], pch=16, add=TRUE, logz=TRUE)


# the boundary of the map? 
# This calculates it as a "Spatial*"
# This will be moved to rgplates!
edge <- chronosphere::mapedge()

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

# Germany 
germany <- countries[which(countries$name=="Germany"), ]
plot(germany$geometry)

# all available attributes
plot(countries)

# one specific attribute
plot(countries["region"])

# overlay
joined <- sf::st_join(countries, sfColls, join=st_intersects)

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

# preliminary plot 
plot(ne$geometry, reset=FALSE, col="gray", xlim=c(34.9, 56.85), ylim=c(-31.92, -8.13))
plot(one$geometry, add=TRUE, col="blue")

# The area of this shape
sf::st_area(one)

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

#################################################################################
# Creating Polygons

# polygon coordinates
poly <- matrix(c(
		 -180, 23.5,
		180, 23.5,
		180, -23.5,
		-180, -23.5, 
		-180, 23.5	
	),
	ncol=2,
	byrow=TRUE)

# Creating a new geometry collection
tropics<- st_geometry(st_polygon(list(poly)))
plot(ne$geometry)
plot(tropics, col="#55555555", add=TRUE)
