# For plotting the occurrences on paleogeographic maps
# Analytical Paleo Workshop
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
# Continuing
library(divDyn)
library(chronosphere)
library(rgplates)
library(sf)
library(terra)
library(viridis)

setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/5_divDyn/spatials")

# for the chronosphere data!
dataDirectory <- "data/chronosphere"
dir.create("data", showWarnings=FALSE)
dir.create(dataDirectory, showWarnings=FALSE )

# the DEMS
dems<- fetch("paleomap", "dem", datadir=dataDirectory)

# the tectonic model
model <- fetch("paleomap", "model", datadir=dataDirectory)

# the stages
data(stages)
stages$map <- matchtime(names(dems), stages$mid)

# the occurrences
datReconstructed <- readRDS(file="export/processed_data.rds")

###############################################################################
# 3. Plotting occurrences
###############################################################################
# loop through every stage -based on the stage ID

# an example stage
i <- 67

# data that come from this slice
stgDat <- datReconstructed[which(datReconstructed$stg==i), ]

# finding the appropriate map
age <- stages$map[i]

########################################----------------------------------------
# 0. Using the plates themselves (for this you could do the stage mid point as well, but the coordinates also should match, then!)
plates <- reconstruct("plates", model=model, age=age)
plates <- as(plates, "sf")
plot(plates$geometry, col="gray", border=FALSE)

# the occurrences
points(
	unique(stgDat[, c("stg_lng", "stg_lat")]),
	pch=16,
	col=stages$col[i]
)

########################################----------------------------------------
# 1. THE DEMS
# the age that corresponds to this, make it terra::SpatRaster from raster::RasterLayer
# for raster
# currentMap<- dems[age]
# for terra
currentMap<- rast(dems[age])
# plot the map
library(viridis)
plot(currentMap, col=viridis(255), axes=FALSE)

# putting the actual occurrences on top of this
points(
	unique(stgDat[, c("stg_lng", "stg_lat")]),
	pch=16,
	col=stages$col[i]
)

########################################----------------------------------------
# 2. Using The Paleomap paleoatlas
pa <- fetch("paleomap", "paleoatlas", datadir=dataDirectory)
pa
# note the ages
stages$map %in% rownames(pa) # to use this series you need to reassign the map ages of the stages!

# it works for the example, but would not for the whole series. 
mapplot(pa[age, ], rgb=TRUE)
points(
	unique(stgDat[, c("stg_lng", "stg_lat")]),
	pch=16,
	col=stages$col[i]
)

########################################----------------------------------------
# 4. Using The Paleomap Paleocoastlines
pc <- fetch("paleomap", "paleocoastlines", datadir=dataDirectory)

# This is derived from the DEMs, so matching of maps to stage midpoints works the same
# But this might not the case in all time series!
stages$map %in% rownames(pc) # there is no Ediacaran in this!

margin <- as(pc[age, "margin"], "sf")
coastlines <- as(pc[age, "coast"], "sf")

# Still not a perfect match!
plot(margin$geometry, col="#006994", reset=FALSE, border=FALSE)
plot(coastlines$geometry, col="#bc815f", add=TRUE, border=FALSE)

points(
	unique(stgDat[, c("stg_lng", "stg_lat")]),
	pch=16,
	col=stages$col[i]
)



##############################
vals <- extract(dems[age], stgDat[, c("stg_lng", "stg_lat")])
