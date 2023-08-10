# Basic regionalization
# Analytical Paleo Workshop
# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
# Continuing
library(divdyn) # basic species cleaning
library(chronosphere) # data download
library(rgplates) # reconstruction
library(sf) # vector spatial
library(icosa) # icosahedral gridding
library(igraph) # graphs/networks

########################################
# workdir
setwd("/mnt/sky/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/5_divDyn/spatials")

# for the chronosphere data!
dataDirectory <- "data/chronosphere"
dir.create(dataDirectory, showWarnings=FALSE )

# th Paleocoastlines
ne <- fetch("NaturalEarth", datadir=dataDirectory)
ne <- as(ne, "sf")

# the occurrences
datReconstructed <- readRDS(file="export/processed_data.rds")

# add the cells
hexa <- hexagrid(7, sp=TRUE)
# the sf representation
gridSF<- as(hexa@sp, "sf")
gridSF$cell <- rownames(centers(hexa))


# selecting the very last 10my bin!
lastBin <- datReconstructed[which(datReconstructed$ten==49),]

# assign cells based on the present day coordinates
lastBin$cell <- locate(hexa, lastBin[,c("lng", "lat")])

# species names - very crude, just for demonstration purposes!
lastBin$sp<- NA
lastBin$sp[lastBin$identified_rank=="species"] <- lastBin$identified_name[lastBin$identified_rank=="species"]
lastBin$sp[lastBin$accepted_rank=="species"] <- lastBin$accepted_name[lastBin$accepted_rank=="species"]

allSpecies <- unique(lastBin$sp)
diag <- divDyn::cleansp(allSpecies, debug=TRUE)
both <- merge(lastBin, diag[, c("original", "binomen")], by.x="sp", by.y="original")

# how much is not species level
sum(is.na(both$binomen))/nrow(both)

########################################
# select only benthic core taxa
need <- which(
	both$class=="Anthozoa"  |
	both$class=="Bivalvia" |
	both$class=="Gastropoda" |
	both$order =="Decapoda" |
	both$phylum=="Brachopoda" |
	both$phylum=="Bryozoa" |
	both$phylum=="Echinodermata"
)
benthic <- both[need,]

# A. network-based
contingency <- table(both$binomen, both$cell)
class(contingency) <- "matrix"

#incidence: presence absence
incidence <- contingency
incidence[incidence>1] <- 1


library(igraph)

# make it a bipartite graph
adj <- graph_from_incidence_matrix(incidence)

# project for locality
loc <- bipartite.projection(adj)[[2]]

# clustering
info <- cluster_infomap(loc)

# membership
mem<-membership(info)

# adding to sf for plotting
gridSF$infomap<- as.character(mem[gridSF$cell])
plot(ne$geometry, border="black", add=TRUE, col=NA)
plot(gridSF["infomap"], reset=FALSE, border="white")

# with the labels!
text(centers(hexa)[names(mem),], label=mem)







