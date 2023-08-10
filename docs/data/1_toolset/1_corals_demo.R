# Script to demonstrate the basic capabilities of R
# Using the divDyn toolchain
# Ádám T. Kocsis, 2021-10-15, Erlangen
# CC BY 4.0

# check whether packages are present, otherwise install them
required <- c("divDyn", "chronosphere", "rgdal")
installed <- installed.packages()

# a for loop: iterate for all packages in required
for(i in 1:length(required)){
	# conditional statement
	if(!required[i]%in%rownames(installed)){
		install.packages(required[i])
	}
}

# setting a project directory
projDir <- NULL
if(!is.null(projDir)) setwd(projDir)

# messages
message("Your current working directory is:")
message(getwd())

# using external packages
library(divDyn)
library(chronosphere)
library(rgdal)

# using built-in data
data(stages) # stage-level timescale
data(corals) # Paleobiology Database download

# executing function calls from packages
dd <- divDyn::divDyn(corals, bin="stg", tax="genus")
stats <- divDyn::binstat(corals, bin="stg", tax="genus")

# calculating simple values and printing them to the console
nCorals <- length(unique(corals$genus))

# concatenating text
statement <-paste(
	"The total number of sampled coral genera is:",
	nCorals,".",  sep="")

# and printing them to the console
message(statement)

# plotting variables against each other (scatter plots)
plot(
	stats$occs, dd$divSIB,
	xlab="Number of occurrences",
	ylab="Sampled-in-bin diversity")

# calculating simple statistics and printing them out to the console
cor.test(stats$occs, dd$divSIB, method="spearman")
message("It is important to understand warning messages!")

# higher and lower level plotting
divDyn::tsplot(stages, boxes="sys", xlim=52:95, ylim=c(0,300), boxes.col="col")
# low level plotting
lines(stages$mid, dd$divSIB, col="black", lwd=2)
lines(stages$mid, dd$divRT, col="blue", lwd=2)

# legend
legend(
	"topleft",
	legend=c("Sampled-in-bin diversity", "Range-through diversity"),
	col=c("black", "blue"),
	lwd=2)

# saving outputs to a directory
# create new directory from R
dir.create("corals_export")

# subsetting tables (matrices)
saveThis <- dd[52:95,]

# writing (and reading tables)
write.csv( saveThis, file="corals_export/divDyn.csv", row.names=FALSE)

# also the original result
saveRDS(file="corals_export/divDyn_original.rds", dd)

# downloading and drawing a world map
ne <-fetch("NaturalEarth")
plot(ne, col="gray", main="Coral occurrences from the world")

# with the coral occurrences
coordinates <- unique(corals[,c("lng","lat")])
# points with RGB colors
points(coordinates, pch=3, col="#AA223355")

# definition of functions
plot_occs_from_stage <- function(x, stage, ts=stages, map=ne){
	# occurrences from the stage
	thisStage <- x[which(x$stg==stage),]

	# a plot
	plot(map, col="gray",
		main=paste("Coral occurrences from the", ts$stage[stage]))

	
	# with the coral occurrences
	coordinates <- unique(thisStage[,c("lng","lat")])
	points(coordinates, pch=3, col=ts$col[stage])
}


# saving plots to the hard drive
png("corals_export/one_stage_occs.png", width=2000, height=1000)
	# calling of functions
	plot_occs_from_stage(corals, 81)	
dev.off()

	
