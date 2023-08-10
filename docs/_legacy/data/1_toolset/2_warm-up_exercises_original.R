# 0. Read in the Palmer penguin dataset!
setwd("C:/Users/Adam/Desktop/penguins")
penguins <- read.csv("data/penguins.csv")
penguins <- read.csv("https://fau-paleo.github.io/apw_2022/data/1_toolset/penguins.csv")

# 1. Select the subset that corresponds to the Gentoo and the Adelie penguins!
# first
gentoo <- subset(penguins,species == "Gentoo")
adelie <- subset(penguins, species == "Adelie")
both <- rbind(gentoo, adelie)
Gen_Ad <- subset(peng, species %in% c("Adelie", "Gentoo"))

x <-c("a", "b", "r", "f")
lookup <- c("a", "b")

x %in% lookup

# measuring code performance
system.time({
  for(i in 1:100){
       mean(1:1000)

  }
})

#
Gen_Ad <- penguins
peng <- penguins
Gen_Ad2 <- peng[Gen_Ad$species == "Adelie" | Gen_Ad$species == "Gentoo", ]

Gen_Ad2 <- peng[which(peng$species == "Adelie" | peng$species == "Gentoo"), ]

# second
Data.Peng1 <- Data.Peng %>% filter(species == "Gentoo" | species == "Adelie")

# operator objects
`%>%`
`%in%`

Gen_Ad2 <- peng[which(Gen_Ad$species == "Adelie" | Gen_Ad$species == "Gentoo"), ]
penguins_Gen_Ade<-subset(penguins, colnames(penguins) == c("Adelie", "Gentoo"))

# NA at eh end does not change this result
which(c(Gen_Ad$species == "Adelie"
| Gen_Ad$species == "Gentoo", NA))

# 2. With a loop for every species, select the species-specific part of the dataset, create a directory
# for it and write it into the new directory in rds and csv format!
# - subset the data by species
# - every file contains data that belong to one species.
sp <- unique(penguins$species)
# sp <- levels(factor(penguins$species))

penguins$species[1] <- NA

for(i in 1:length(sp)){
#for(i in seq_along(sp))

    # in every
    # create a new directory the species
    currentSpecies <- sp[i]
    dir.create(currentSpecies)
    # get the subset that contains all data of the species
    speciesData <- penguins[which(penguins$species == currentSpecies),]

    # write it into the directories
    # write.csv
    write.csv(speciesData, file=
        paste0(currentSpecies, "/", currentSpecies, ".csv"),
        row.names=FALSE)

    # "Chinstrap/data.csv"
    # use this "Chinstrap/Chinstrap.csv"

    # saveRDS
    saveRDS(speciesData,
        file=paste0(currentSpecies, "/", currentSpecies, ".rds"))
}

adelie <- readRDS("Adelie/Adelie.rds")

# .rds
# .RData
save(adelie, file="Adelie/Adelie.RData")
rm(adelie)
load(file="Adelie/Adelie.RData")

# not to do!
# o <- load(file="Adelie/Adelie.RData")

# 3. Change the previous code to save the subsets in a list. Use the species names as the names of the list!
li <- list(a=1:10, b=c("a", "b", "c"))

length(li)
li[1]
li[[1]]
li$a

container <- list()
for(i in 1:length(sp)){

    # in every
    currentSpecies <- sp[i]

    # get the subset that contains all data of the species
    speciesData <- penguins[which(penguins$species == currentSpecies),]

    # write it into the directories
    container[[i]] <- speciesData

}

str(container)
names(container) <- sp

container$Chinstrap
container[["Chinstrap"]]
container[["Chinstrap"]][1, ]



# 4. Write a function that returns parts of the dataset based on column name and value,
# either greater or equal/lower or lower or equal to the value. function(x, column, value, greater)
#' @param x The data frame to subset
#' @param column Which column you are basing the subsetting on. (Numeric column)
#' @param value A numeric value, that separates the subset that you want to select from what you don't want.
#' @param greater Should a subset be selected that correponds to rows where the values are
#' greater than \code{value} (TRUE), or lower tan equal (FALSE).
