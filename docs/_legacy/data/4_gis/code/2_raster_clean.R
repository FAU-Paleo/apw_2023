# Basic raster data examples
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/4_GIS/")

# necessary packages
library(viridis) # color palette
library(terra) # actual
library(divDyn) # fossils
library(sf) # interactions between

# 0. primer: The volcano data
data(volcano)
str(volcano)
?volcano

# plotting such a data 
contour(volcano)
filled.contour(volcano, col=viridis(25))

################################################################################
# 1. Raster data - single layer
# reading in a file, instantiate a raster object. raster::raster()
etopo1 <- rast("data/ETOPO1/ETOPO1_ice_c_20110606_tiff_1.tif")
etopo1

# crs as WKT
cat(terra::crs(etopo1))

# basic attributes - raster::res(), raster::xres(), raster::yres()
terra::res(etopo1)
terra::xres(etopo1)
terra::yres(etopo1)

# the extent - raster::extent()
terra::ext(etopo1)


# plotting
plot(etopo1)

# fully custommized plot
plot(etopo1
	, xlim=c(-45, 45) # limits to plotting
	, ylim=c(-45, 45) # y limit to plotting
	, col=viridis(255) # color palette
	, axes=FALSE # an axes
	, main="This is a plot"
)

# Overplotting something
data(corals)
colls <- unique(corals[, c("collection_no", "lng", "lat")])

# as regular points
points(colls[, c("lng", "lat")]
	, pch=3
	, col="red"
)

# Cropping: Changing the extent of the raster 
ex <- terra::ext(
	c(
		xmin=-45,
		xmax=+45,
		ymin=-30, 
		ymax=30
	)
)

# the actual cropping raster::crop()
cropped<- terra::crop(etopo1, ex)
plot(cropped)

# this is different here
wrong<- etopo1
terra::ext(wrong) <- c(30, 40, 10, 67)
plot(wrong)

# Changing raster resolution
wrong<- etopo1
terra::res(wrong) <- c(0.5,0.5)

# 1A. Aggregation: decrease resolution by combining cells, raster::aggregate()
aggregated <- terra::aggregate(etopo1, 2)
aggregated <- terra::aggregate(etopo1, 10)
plot(aggregated)

# 1B. Disaggregation: increase resolution by splitting cells, raster::disaggregate()
disaggregated <- terra::disagg(etopo1, 4)
plot(disaggregated)

# 2A. Resampling
# increase resolution (smoothing)
# Step 1: define a raster object. 
target <- terra::rast(res=c(0.2, 0.2))

# Step 2: actual resampling
retopo1 <- terra::resample(etopo1, target)

# decrease resolution, same steps
target <- terra::rast(res=c(1.3, 1.3))
retopo1 <- terra::resample(etopo1, target)
plot(retopo1)

# Projection change
et1Moll <- terra::project(etopo1, "ESRI:54009")
et1Moll
plot(et1Moll)

# missing values outside the projection!
valsMoll<- terra::values(et1Moll)
sum(is.na(valsMoll))/terra::ncell(et1Moll)

# Rasterization: iteration based on raster structure.
# eg. Count the number of items in cells
grid <- terra::rast(res=c(10, 10))
rasted <- terra::rasterize(
	x=as.matrix(colls[, c("lng", "lat")]), 
	y= grid, 
	fun=length)

# the density of occurrences
plot(rasted, col=viridis::viridis(256))

# Read in the coutnries data 
countries <- sf::st_read("data/world_admin/world-administrative-boundaries.shx")
plot(countries$geometry, col=NA, border="black", add=TRUE)

# making spatial points
collections <- unique(corals[, c("collection_no", "lng", "lat")])

# the number of occurrences
nOccs <- table(corals$collection_no)
collections$nOccs <- nOccs[as.character(collections$collection_no)]

# omit missing values
collections <- collections[!is.na(collections$lng) & !is.na(collections$lat), ]

# as simple feature!
sfColls <- sf::st_as_sf(collections, coords=c("lng", "lat"))
plot(sfColls$geometry, add=TRUE, pch=3, col="#99000033")

# masking with countries
countryElev <- terra::mask(etopo1, countries)
plot(countryElev)

# one country's elevation
canada <- countries[countries$name=="Canada",]
plot(canada$geometry, col="red")

plot(etopo1)
plot(canada$geometry, add=TRUE)

oneElev <- terra::mask(etopo1, canada)
plot(oneElev)

# correct from implementation point of view, but latitudinal patterns might be a problem!
canadaVals <- terra::values(oneElev)  
mean(canadaVals, na.rm=T)

################################################################################
# Multi-level raster
wcFiles <- list.files("data/WorldClim")
wc <- terra::rast(file.path("data/WorldClim",wcFiles))
plot(wc)

# extracting values
erlangen <- matrix(c(10.97, 49.58),  ncol=2)
terra::extract(etopo1, erlangen)

# average values
terra::extract(wc, erlangen)

# single layer
bio1 <- wc[[1]]
bio1 <- wc[["wc2.1_10m_bio_1"]]

terra::extract(bio1, erlangen)
plot(bio1)


# Lower-level operations
# replacing values, only above 1000m level
high<- etopo1
terra::values(high)[terra::values(high)<1000] <- NA
plot(high)

# Corresponding climate
bio1Re <- terra::resample(bio1, high)
bio1High <- terra::mask(bio1Re, high)
plot(bio1High)

################################################################################
# Climate averaging
# Latitude averaging
oneBio1<- terra::mask(bio1, canada)
plot(oneBio1)

# doing this in long-lat (wrong, polar areas are overrepresented!)
canadaVals <- terra::values(oneBio1)  
mean(canadaVals, na.rm=T)

# in equal-area projection (correct)
oneBio1<- terra::project(oneBio1, "ESRI:54009")
canadaVals <- terra::values(oneBio1)  
mean(canadaVals, na.rm=T)

################################################################################
# Getting x and y coordinates
xy <- terra::xyFromCell(etopo1, 1:terra::ncell(etopo1))
xyz <- cbind(xy, terra::values(etopo1))


