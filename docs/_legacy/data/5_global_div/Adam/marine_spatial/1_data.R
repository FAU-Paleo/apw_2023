# Basic code for the Phanerozoic-scale, stage-level data preparetion for marine Genera
# Analytical Paleo Workshop
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
################################################################################
# Necessary packages
## library(chronosphere) # data distribution
## library(divDyn) # temporal binning
## library(rgplates) # reconstruction
## library(sf) # for examples
## library(terra) # for plotting
## library(viridis) # for plotting
################################################################################
# A. Loading data
################################################################################
# read the data table
	library(chronosphere)
	dat <- fetch("pbdb", ver="20220510")

# 1. taxonomic filtering
	# filter records not identified at least to genus
	dat <-dat[dat$accepted_rank %in% c("genus", "species"),]

	# omit non-informative genus entries
	dat <- dat[dat$genus!="", ]

	#A. phyla
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

	# the other
		noNeed <- dat[!bByPhyla,]
		needPhylum <- dat[bByPhyla,]

	#B. classes
	#	levels(factor(noNeed$class))
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
	#	mammals <- dat[dat$class=="Mammalia",]
	#	levels(factor(mammals$order))

		needMammalOrd <- c("Cetacea", "Sirenia")

		# which rows?
		bMammalOrder <- dat$order %in% needMammalOrd

		# the carnivores
		carnivores <- dat[dat$order=="Carnivora",]
		levels(factor(carnivores$family))

		needFam <- c("Otariidae", "Phocidae", "Desmatophocidae")

		# which rows?
		bNeedMamFam<- dat$family%in%needFam

	# D. Reptiles
	#	reptiles <- dat[dat$class=="Reptilia",]
	#	levels(factor(reptiles$order))

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
	#	turtles <- dat[dat$order=="Testudines",]
	#	levels(factor(turtles$family))
	
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
	dat <- dat[
		bByPhyla |
		bNeedClass |
		bMammalOrder |
		bNeedMamFam |
		bRept | 
		bTurtle
		, ]


	# resolve the potential homonymy problem
	dat$clgen <- paste(dat$class, dat$genus)

################################################################################
# 2. filter by environment
	levels(factor((dat$environment)))

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

	dat<-dat[!dat$environment%in%omitEnv, ]


# finally omit unlithified sediments
	dat <- dat[dat$lithification1!="unlithified",]

	
# resolving remaining marine environmental variables (not for this, but for additional analyses) 
	library(divDyn)
 	#load from divDyn
 	data(keys)

 	# lithology
	dat$lith<-categorize(dat$lithology1,keys$lith)

	# batyhmetry
	dat$bath <- categorize(dat$environment,keys$bath) 

	# grain size
	dat$gra <- categorize(dat$lithology1,keys$grain) 

	# reef or not?
	dat$reef <- categorize(dat$environment, keys$reef) 
	dat$reef[dat$lith=="clastic" & dat$environment=="marine indet."] <- "non-reef" # reef or not?/2

	# onshore - offshore
	dat$depenv <- categorize(dat$environment,keys$depenv) 

################################################################################
# 3. stratigraphic resolution
	# time scales
	data(stages)
	data(tens)
	
	# rename some names so they are the same (Cosmetic changes, needed for the divDyn-analyses)
	colnames(tens)[colnames(tens)=="X10"]<-"name"
	colnames(stages)[colnames(stages)=="stage"]<-"name"

	
#-------------------------------------------------------------------------------
# A.  the bin entries

	# a. categorize interval names to bin numbers
		# categorize is the new function of the package
		tenMin<-categorize(dat[,"early_interval"],keys$tenInt)
		tenMax<-categorize(dat[,"late_interval"],keys$tenInt)

		# convert to numeric
		tenMin<-as.numeric(tenMin)
		tenMax<-as.numeric(tenMax)

	# b. resolve max-min interval uncertainty
	# final variable (empty)
		dat$ten <- rep(NA, nrow(dat))

	# use entries, where
		tenCondition <- c(
		# the early and late interval fields indicate the same bin
			which(tenMax==tenMin),
		# or the late_interval field is empty
			which(tenMax==-1))

	# in these entries, use the bin indicated by the early_interval
		dat$ten[tenCondition] <- tenMin[tenCondition]

 # Sampling
	sampTen <- binstat(dat, bin="ten", tax="clgen", coll="collection_no", 
		duplicates=FALSE)  

	# sort the number of occs/ten
	tenOccs <- sampTen$occs
	names(tenOccs) <- rownames(sampTen)
	sort(tenOccs)

	# plot them
	# time scale
	tsplot(stages, boxes="sys", shading="sys", xlim=4:95, ylim=c(0,35000), 
		ylab="number of entries", plot.args=list(cex.lab=1.5, cex.axis=1.8),
		labels.args=list(cex=1.5))

	# occurrences
	lines(tens$mid, sampTen$occs, lwd=2)

	# collections
	lines(tens$mid, sampTen$colls, lwd=2, col="blue")

	# legend
	legend("topright", bg="white", legend=c("occurrences", "collections"), 
		col=c("black", "blue"), lwd=2, inset=c(0.1,0.01), cex=1.3)
	
#####################################-------------------------------------------
# do the same for the stages
# B. the stg entries (lookup)
		stgMin<-categorize(dat[,"early_interval"],keys$stgInt)
		stgMax<-categorize(dat[,"late_interval"],keys$stgInt)

		# convert to numeric
		stgMin<-as.numeric(stgMin)
		stgMax<-as.numeric(stgMax)

	# empty container
		dat$stg <- rep(NA, nrow(dat))

	# select entries, where
		stgCondition <- c(
		# the early and late interval fields indicate the same stg
			which(stgMax==stgMin),
		# or the late_intervarl field is empty
			which(stgMax==-1))

	# in these entries, use the stg indicated by the early_interval
		dat$stg[stgCondition] <- stgMin[stgCondition]

	# terrible in the pre-Silurian
# using the online items
	# additional treatment required for Cambrian
		# load data
		load(url("https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/cambStrat.RData"))
		# correct it with this function
		source("https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2018-08-31/cambProcess.R")

	# and the Ordovician
		# load data
		load(url("https://github.com/divDyn/ddPhanero/raw/master/data/Stratigraphy/2018-08-31/ordStrat.RData"))
		# correct it with this function
		source("https://github.com/divDyn/ddPhanero/raw/master/scripts/strat/2019-05-31/ordProcess.R")

# Sampling
	sampStg<- binstat(dat, bin="stg", tax="clgen", coll="collection_no", 
		duplicates=FALSE)  

	# time scale
	tsplot(stages, boxes="sys", shading="sys", xlim=4:95, ylim=c(0,20000), 
		ylab="number of entries", plot.args=list(cex.lab=1.5, cex.axis=1.8),
		labels.args=list(cex=1.5))

	# occurrences
	lines(stages$mid, sampStg$occs, lwd=2)

	# collections
	lines(stages$mid, sampStg$colls, lwd=2, col="blue")

	# legend
	legend("top", bg="white", legend=c("occurrences", "collections"), 
		col=c("black", "blue"), lwd=2, inset=c(0.15,0.01), cex=1.3)
	
# example metrics
	ddStages<-divDyn(dat, bin="stg", tax="clgen")

	# basic T-diversity curve
	tsplot(stages, boxes="sys", shading="sys", xlim=4:95, ylim=c(0, 4000))
	lines(stages$mid, ddStages$divRT, col="black")	

################################################################################
# 2. Reconstruction of paleocoordinates to stage midpoints
################################################################################
# Reconstructions are implemented in the package 'rgplates'
library(rgplates)

# Although the PBDB has paleocoordinates, they might not be perfect for some applications.
# 
# Two things define the exact values of the reconstructed coordinates
# - 1. The model used for the reconstruction. This is primarily a rotation file, but also a set of polygons for the plates
# are also necessary, otherwise you don't know how to tie the rotation file to your coordinates. These two are usually used
# and distributed together. The model used to tie spatial data to collections has to match the reconstruction model. 
# - 2. The age of the collection. The collections moved with the plates, and the paleocoordinates represent the expected age of the collection.
# Spatial data tends to have a coarse temporal resolution in the past, for instance every 5 million years. Check out these

# The PaleoDEMs
dems<- fetch("paleomap", "dem")

# The PaleomapPaleoatlas
pa<- fetch("paleomap", "paleoatlas")

# When we put multiple things on the same map, we imply that their age is the same - even though it isn't.
# This is dependent on the temporal gridding. The easiest solution is to say that each collection that falls into the
# same bin is reconstructed to the same age, which matches a spatial reconstruction.

# The simplest way to match these up is to 1. bin the data (which we already have), and 2. assign a map to every bin.
# This we can do with the convenience function matchtime().

# Let's use a framework that is similar to the PaleoDEMs, every 5 million years
# the ages of the dems
names(dems)

#reorder these so they align with the bin numbers, based on whichever is closest to the stage mid age.
# This information is specific to the time scale that we are using, so let's save it there!
stages$map <- matchtime(names(dems), stages$mid)

# this information is characteristic to our analyses now: whenever we refer to particular temporal bin, we use a
# distinct spatial reconstruction with it

# this means that the collections have to be reconstructed to these ages, when we want them to interact with our spatial data
# Let's separate the collection-level information from the occurrences
collections <- dat[ , 
	c(NULL
		, "collection_no" # unique id
		, "lng" # current coordinates
		, "lat" # current coordinates
		, "stg" # thetemporal bin
		, "paleolat" # the original paleocoord (for comparison)
		, "paleolng" # the original paleocoord (for comparison)
	)
]

# make it collection-level
collections <- unique(collections)

# it takes a lot of time to reconstruct the paleocoordinates with the web application
# for this reason, rgplates has an offline method for reconstructions, which uses the console interface of GPlates
# However for that to work we need plate reconstruction model files. These can be acquired from chronosphere.
# For the PaleoMAP project this is:
model <- fetch("paleomap", "model")

# the reconstruct works the same way as it does with the offline method. For instance to plot Erlangen on the plate polygons you would
	# This is also a test to see if gplates is accessible!
	# the Plates
	plates150 <- reconstruct("plates", age=150, model=model)
	erlangen150 <- reconstruct(matrix(c(11.0, 49.6), ncol=2), age=150, model=model)

	# plot as an sf!
	library(sf)
	plates150 <- as(plates150, "sf")
	plot(plates150$geometry, col="gray", border=NA)
	points(erlangen150, col="red", pch=16)


# to effectively reconstruct the coordinates we need to put the target ages to our tables
# for every stg id we have a mapage, we just need to copy this from the time scale table to the collections
collections$stg_map<- stages$map[collections$stg] # based on the rule that stg 5 is in the 5th row!

# it is also a good idea to do a bit of cleaning
# 1. We want to focus on those data the have stage assignments!
collections <- collections[!is.na(collections$stg), ]

# and valid present-day coordinates
collections <- collections[!is.na(collections$lat) & !is.na(collections$lng), ]

# weird entries? 
range(collections$lng)
range(collections$lat)
collections <- collections[collections$lat<=90, ]

# now we can execute the reconstruction, this takes a couple of minutes (3-5). 
system.time(
paleoCoords <- reconstruct(
	collections[, c("lng", "lat")] # the coordinates to reconstruct
	, age = collections$stg_map # the target reconstruction age, every coordinate pair (row in matrix), has one value
	, model=model # the R representation of the the platemodel
	, enumerate=FALSE # The way how the coordinates are pair with the ages. Enumerate FALSE means the 'paired' approach, described for age
	, plateperiod=FALSE # Avoid missing values produced by plate model inaccuracies, See ?reconstruct and 'plateperiod'
)
)

# The order of the coordinates is the same, 
colnames(paleoCoords) <- c("stg_lng", "stg_lat")
collections <- cbind(collections, paleoCoords)

# now we can contrast our results with what we had from the PBDB
plot(collections$stg_lat, collections$paleolat)

# merge data back with the occurrence dataset!
datReconstructed <- merge(
	dat,
	collections[, c("collection_no", "stg_map", "stg_lng", "stg_lat")], 
	by="collection_no",
	all=TRUE
)


# Now is a good idea to save what we have done!
setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/5_divDyn/")
dir.create("export", showWarning=F)
saveRDS(datReconstructed, file="export/processed_data.rds")




