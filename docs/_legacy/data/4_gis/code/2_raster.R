setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/4_GIS/")
# Raster data
library(viridis) # color palette
library(terra) # 
library(divDyn) # 
library(sf)

# 0. primer: The volcano data
data(volcano)
str(volcano)

# plotting such a data 
contour(volcano)
filled.contour(volcano, col=viridis(25))

################################################################################
# 1. Raster data - single layer
# reading in a file, instantiate a raster object. raster::raster()
etopo1 <- rast("data/ETOPO1/ETOPO1_ice_c_20110606_tiff_1.tif")
etopo1

# basic attributes - raster::res(), raster::xres(), raster::yres()
res(etopo1)
xres(etopo1)
yres(etopo1)

# the extent - raster::extent()
ext(etopo1)

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
ex <- ext(
	c(
		xmin=-45,
		xmax=+45,
		ymin=-30, 
		ymax=30
	)
)

# the actual cropping raster::crop()
cropped<- crop(etopo1, ex)
plot(cropped)

# Changing raster resolution
# 1A. Aggregation: decrease resolution by combining cells, raster::aggregate()
aggregated <- aggregate(etopo1, 2)
aggregated <- aggregate(etopo1, 3)
plot(agg)

# 1B. Disaggregation: increase resolution by splitting cells, raster::disaggregate()
disaggregated <- disagg(etopo1, 4)
plot(disaggregated)

# 2A. Resampling
# increase resolution (smoothing)
# Step 1: define a raster object. 
target <- rast(res=c(0.2, 0.2))

# Step 2: actual resampling
retopo1 <- resample(etopo1, target)

# decrease resolution, same steps
target <- rast(res=c(1.3, 1.3))
retopo1 <- resample(etopo1, target)
plot(retopo1)

# Projection change
et1Moll <- project(etopo1, "ESRI:54009")
et1Moll
plot(et1Moll)

# Rasterization: iteration based on raster structure.
# eg. Count the number of items in cells
grid <- rast(res=c(10, 10))
rasted <- rasterize(as.matrix(colls[, c("lng", "lat")]), grid, fun=length)

# the density of occurrences
plot(rasted, col=viridis(256))

# Read in the coutnries data 
countries <- sf::st_read("data/world_admin/world-administrative-boundaries.shx")
plot(countries$geometry, col=NA, border="black", add=TRUE)
plot(sfColls$geometry, add=TRUE, pch=3, col="#99000033")

# masking with countries
countryElev <- mask(etopo1, countries)

# one country's elevation
oneElev <- mask(etopo1, countries[countries$name=="Canada",] )
plot(oneElev)



# extracting values
erlangen <- matrix(c(10.97, 49.58),  ncol=2)
extract(etopo1, erlangen)

################################################################################
# Multi-level raster
wcFiles <- list.files("data/WorldClim")
wc <- rast(file.path("data/WorldClim",wcFiles))

# average values
extract(wc, erlangen)

# single layer
bio1 <- wc[[1]]
plot(bio1)


# Lower-level operations
# Accessing values
vals <- values(etopo1)

# replacing values, only above 1000m level
high<- etopo1
values(high)[values(high)<1000] <- NA
plot(high)

# Corresponding climate
bio1Re <- resample(bio1, high)
bio1High <- mask(bio1Re, high)
plot(bio1High)


################################################################################
# Getting x and y coordinates
xy <- xyFromCell(etopo1, 1:ncell(etopo1))
xyz <- cbind(xy, values(land))


