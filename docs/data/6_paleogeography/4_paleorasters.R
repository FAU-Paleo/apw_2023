# Basic examples with Raster-Paleogeography 
# Erlangen, 2023
# Ádám T. Kocsis
# CC-BY (attribution)
library(rgplates)
library(chronosphere)
library(divDyn)
library(via)
library(sf)
library(icosa)
library(viridis)


# The PaleoMAP model - Chris Scotese's reconstruction
# Frequently used in for environmental reconstruction

# needed:
# install.packages("ncdf4")
# install.packages("terra")

# via series of rasters
# R-array interface to multi-dimensional rasters
paleodems <- chronosphere::fetch("paleomap", "dem")
paleodems

# present-day topography
plot(paleodems['0'])

# 200 Ma
plot(paleodems['200'])

# corresponding model
mod <- chronosphere::fetch("paleomap", "model")
mod

# the plates on top of the DEMs
plates200 <- rgplates::reconstruct("static_polygons", age=200, model=mod)
plot(plates200, col="#AA000077", add=TRUE)


################################################################################
# Some example data
########################################----------------------------------------
# Updated coastlines - based on DEMs
pc <- chronosphere::fetch("paleomap", "paleocoastlines")
pc
dimnames(pc)

# updated 
plot(pc["70","margin"], col="blue", border=NA)
plot(pc["70","coast"], col="brown", border=NA, add=TRUE)

########################################----------------------------------------
# BRIDGE climate model data (Many more coming soon!) - also based on PaleoMAP
# Interpolate to match the paleocoastlines!
sst <- chronosphere::fetch("SOM-kocsis-provinciality")
plot(sst["70"])
plot(pc["70","coast"], col="black", border=NA, add=TRUE)

########################################----------------------------------------
li <- chronosphere::fetch("SOM-Li-PhanerozoicClimate")
dimnames(li)

# january TAS (near surface air temperature at 70Ma
one <- li["70", "jan", "T"]
plot(terra::rotate(one))

# updates done in one year!
plot(pc["70","coast"], col="#00000055", border=NA, add=TRUE)

################################################################################
# II. Working with PaleoDB data

# collections are not contemporaneous! Discretization is necessary
# Assign collections to bins -> reconstruct of them to the same target date.

# same as earlier for divDyn
data(stages)

# assign an age - whichever is closest!
stages$map <- divDyn::matchtime(names(sst), stages$mid)

## # earlier version of the timesclae
## data(stages2018)
## stages$midOld <-stages2018$mid 

# Data download with explicit version
dat <- chronosphere::fetch("pbdb")

# filter records not identified at least to genus
dat <-dat[dat$accepted_rank %in% c("genus", "species"),]

# omit non-informative genus entries
dat <- dat[dat$genus!="", ]

marineNoPlant <- c("",
    "Agmata",
    "Annelida",
    "Bilateralomorpha",
    "Brachiopoda",
    "Bryozoa",
    "Calcispongea",
    "Chaetognatha",
    "Cnidaria",
    "Ctenophora",
    "Echinodermata",
    "Entoprocta",
    "Foraminifera",
    "Hemichordata",
    "Hyolitha",
    "Mollusca",
    "Nematoda",
    "Nematomorpha",
    "Nemertina",
    "Onychophora",
    "Petalonamae",
    "Phoronida",
    "Platyhelminthes",
    "Porifera",
    "Problematica",
    "Rhizopodea",
    "Rotifera",
    "Sarcomastigophora",
    "Sipuncula",
    "Uncertain",
    "Vetulicolia",
    ""
)

# which rows?
bByPhyla <- dat$phylum%in% marineNoPlant
#B. classes
#   levels(factor(noNeed$class))
needClass <- c(
    "Acanthodii",
    "Actinopteri",
    "Actinopterygii",
    "Agnatha",
    "Cephalaspidomorphi",
    "Chondrichthyes",
    "Cladistia",
    "Coelacanthimorpha",
    "Conodonta",
    "Galeaspida",
    "Myxini",
    "Osteichthyes",
    "Petromyzontida",
    "Plagiostomi",
    "Pteraspidomorphi",
    # here come the Arthropods
    "Artiopoda",
    "Branchiopoda",
    "Cephalocarida",
    "Copepoda",
    "Malacostraca",
    "Maxillopoda",
    "Megacheira",
    "Merostomoidea",
    "Ostracoda",
    "Paratrilobita",
    "Pycnogonida",
    "Remipedia",
    "Thylacocephala",
    "Trilobita",
    "Xiphosura"
)

# which rows?
bNeedClass <- dat$class %in% needClass
#C.  mammals
#   mammals <- dat[dat$class=="Mammalia",]
#   levels(factor(mammals$order))
needMammalOrd <- c("Cetacea", "Sirenia")

# which rows?
bMammalOrder <- dat$order %in% needMammalOrd

# the carnivores
#   carnivores <- dat[dat$order=="Carnivora",]
#   levels(factor(carnivores$family))
needFam <- c("Otariidae", "Phocidae", "Desmatophocidae")

# which rows?
bNeedMamFam<- dat$family%in%needFam

# D. Reptiles
#   reptiles <- dat[dat$class=="Reptilia",]
#   levels(factor(reptiles$order))
needReptOrd<-c(
    "Eosauropterygia",
    "Hupehsuchia",
    "Ichthyosauria",
    "Placodontia",
    "Sauropterygia",
    "Thalattosauria"
)

# which rows?
bRept <- dat$order%in%needReptOrd

# E. turtles 
#   turtles <- dat[dat$order=="Testudines",]
#   levels(factor(turtles$family))

# Chelonioidea turtles
needTurtleFam <- c(
    "Cheloniidae",
    "Protostegidae",
    "Dermochelyidae",
    "Dermochelyoidae",
    "Toxochelyidae",
    "Pancheloniidae"
)

# which rows?
bTurtle <- dat$family%in%needTurtleFam

# subset the data
dat <- dat[bByPhyla | bNeedClass | bMammalOrder | bNeedMamFam | bRept | bTurtle , ]
dat$clgen <- paste(dat$class, dat$genus)
omitEnv <- c(
    "\"floodplain\"",
    "alluvial fan",
    "cave",
    "\"channel\"",
    "channel lag" ,
    "coarse channel fill",
    "crater lake",
    "crevasse splay",
    "dry floodplain",
    "delta plain",
    "dune",
    "eolian indet.",
    "fine channel fill",
    "fissure fill",
    "fluvial indet.",
    "fluvial-lacustrine indet.",
    "fluvial-deltaic indet.",
    "glacial",
    "interdune",
    "karst indet.",
    "lacustrine - large",
    "lacustrine - small",
    "lacustrine delta front",
    "lacustrine delta plain",
    "lacustrine deltaic indet.",
    "lacustrine indet.",
    "lacustrine interdistributary bay",
    "lacustrine prodelta",
    "levee",
    "loess",
    "mire/swamp",
    "pond",
    "sinkhole",
    "spring",
    "tar",
    "terrestrial indet.",
    "wet floodplain"
)

# actual omission
dat<-dat[!dat$environment%in%omitEnv, ]
dat <- dat[dat$lithification1!="unlithified",]
nrow(dat)
data(keys)
# B. the stg entries (lookup)
stgMin <- divDyn::categorize(dat[,"early_interval"],keys$stgInt)
stgMax <- divDyn::categorize(dat[,"late_interval"],keys$stgInt)
stgMin <- as.numeric(stgMin)
stgMax <- as.numeric(stgMax)
# empty container
dat$stg <- rep(NA, nrow(dat))
# select entries, where
stgCondition <- c(
# the early and late interval fields indicate the same stg
    which(stgMax==stgMin),
# or the late_intervarl field is empty
    which(stgMax==-1))

dat$stg[stgCondition] <- stgMin[stgCondition]

# load data
load(url("https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/cambStrat.RData"))

# correct it with this function
source("https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2018-08-31/cambProcess.R")
# load data
load(url("https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/ordStrat.RData"))

# correct it with this function
source("https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2019-05-31/ordProcess.R")


########################################----------------------------------------
# reconstruct coordinates -
collections <- unique(dat[, c("collection_no", "lng", "lat", "paleolng", "paleolat", "stg")])

# assign values - in this case to match that of the SST dataset
collections$map <- stages$map[collections$stg]

# with rgplates -> enumeration by default (takes a couple of minutes!)
paleoCoords <- rgplates::reconstruct(collections[, c("lng", "lat")], age=collections$map, model=mod, enumerate=FALSE)
colnames(paleoCoords) <- c("plng", "plat")

# cross check with default downloads!
plot(collections[,"paleolng"], paleoCoords[, "plng"])
plot(collections[, "paleolat"], paleoCoords[, "plat"])

# merge back
colls <- cbind(collections, paleoCoords)
pbdb <- merge(dat, colls[, c("collection_no", "plng", "plat")], by="collection_no")


########################################----------------------------------------
# Extract values from rasters

# plot a random map
i <- 87

# the model data
age <- stages$map[i]

# character subscript!!
curMap <- sst[age]
plot(curMap, main=paste0(stages$stage[i], ", ", age, " Ma" ))

# current dataset
curDat <- colls[which(colls$stg==i), ]

# still not perfect, loads to do still
points(curDat[, c("plng", "plat")], pch=3)

# extract sst values
vals <- extract(curMap, curDat[, c("plng", "plat")])
curDat$sst <- vals[, names(curMap)]

# cross plot with latitude
plot(curDat$plat, curDat$sst, xlab="Latitude (°)", xlim=c(-90,90), ylab="Annual mean SST (°C)")


################################################################################
# Gridding data with icosa. More info at https://adamkocsis.github.io/icosa/

# original icosahedron
tg <- icosa::trigrid(1)
tg

# install rgl to make a 3d plot (might need restart R)
# install.packages("rgl")
icosa::plot3d(tg)

# tessellation
tg <- icosa::trigrid(3)
plot3d(tg)

# the inversion
hg <- icosa::hexagrid(1)
hg

# create an SF representation for plotting
hg <- icosa::hexagrid(4,sf=TRUE)

# the sf-plot
plot(hg)

# degree approximation
hg <- icosa::hexagrid(deg=2.5, sf=TRUE)

# Present-day sampling
ne <- chronosphere::fetch("NaturalEarth")
plot(ne$geometry, col="gray", border=NA)
plot(hg, add=TRUE)

# Present-day sampling
# actual error in present-day coords
collections$cells <- icosa::locate(hg, collections[, c("lng", "lat")])

# omit invlaid 
collections <- collections[which(abs(collections$lat)<=90),]

# try again
collections$cells <- icosa::locate(hg, collections[, c("lng", "lat")])

# Example use: 
# number of collections in cells
collTab <- table(collections$cells)

# sf-style plotting readily available
plot(hg, collTab)

# multi-layers
plot(hg, collTab, logz=TRUE, border="white", reset=FALSE)
plot(ne$geometry, add=TRUE, col="#00000033", border="black")

# manual ramping
breaks <- c(0, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
colors <- viridis::viridis(length(breaks)-1, alpha=0.5)

plot(hg, collTab, breaks=breaks, pal=colors, border="#ffffff80", reset=FALSE)
plot(ne$geometry, add=TRUE, col="gray", border="gray")
plot(hg, collTab, breaks=breaks, pal=colors, border="#ffffff80", add=TRUE)











