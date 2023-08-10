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

setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/5_divDyn/spatials")

# for the chronosphere data!
dataDirectory <- "data/chronosphere"
dir.create(dataDirectory, showWarnings=FALSE )

# the DEMS
dems<- fetch("paleomap", "dem", datadir=dataDirectory)

# the tectonic model
model <- fetch("paleomap", "model", datadir=dataDirectory)

# th Paleocoastlines
pc <- fetch("paleomap", "paleocoastlines", datadir=dataDirectory)

# the stages
data(stages)
stages$map <- matchtime(names(dems), stages$mid)

# the occurrences
datReconstructed <- readRDS(file="export/processed_data.rds")

###############################################################################
# loop through every stage -based on the stage ID

# an example stage
i <- 81

# data that come from this slice
stgDat <- datReconstructed[which(datReconstructed$stg==i), ]

# finding the appropriate map
age <- stages$map[i]

# simple plotting
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




###############################################################################
# icosahedral griding
# install.packages("icosa")
library(icosa)
data(tessguide)
tessguide[1:20,]

# a 10-degree mean edge-length grid
tri <- trigrid(2, sp=TRUE)
plot(tri)

hexa <- hexagrid(c(4, 8), sp=TRUE)
hexa

# two dimensional representation
plot(hexa)
gridlabs(hexa)

# detailed plots
gridSF<- as(hexa@sp, "sf")
gridSF$cell <- rownames(centers(hexa))

# grid midpoints
centers(hexa)

# coordinate lookup:
stgDat$cell <- locate(hexa, stgDat[,c("stg_lng", "stg_lat")])

###############################################################################
# Once this is done, everything is easy!
# A. Sampling in space: occurrences in a cell
tab <- table(stgDat$cell)

# add it to the sf object!
gridSF$occs <- tab[gridSF$cell]

# just this info
plot(gridSF["occs"], reset=FALSE, border="white")
plot(coastlines$geometry, add=TRUE, border="black", col="#55555544")

points(
	unique(stgDat[, c("stg_lng", "stg_lat")]),
	pch=16,
	col="#FF000044"
)

########################################----------------------------------------
# B. diversity in a grid cell
richness <- tapply(INDEX=stgDat$cell, X=stgDat$clgen, function(x) length(unique(x)))

# add it to the sf for plotting
gridSF$richness <- richness[gridSF$cell]

# richness plot
plot(gridSF["richness"], reset=FALSE, border="white")
plot(coastlines$geometry, add=TRUE, border="black", col="#55555544")

########################################----------------------------------------
# C. Geographic range measurement
# e.g. occupancy
occupancies <- tapply(INDEX=stgDat$clgen, X=stgDat$cell, function(x) length(unique(x)))
hist(occupancies)

