# Code from the second day of the 2022 APW
# Erlangen, 2022-08-23


# 0. Read in the Palmer penguin dataset!
# setwd("C:/Users/Adam/Desktop/penguins")

# either from the hard drive
# penguins <- read.csv("data/penguins.csv")

# or from the internet
penguins <- read.csv("https://fau-paleo.github.io/apw_2022/data/1_toolset/penguins.csv")

### NOTE: This is also published as a package:
## install.packages("palmerpenguins")
## library(palmerpenguins)
## data(penguins)


# 1. Select the subset that corresponds to the Gentoo and the Adelie penguins!
# Method 1: separate subsets and rejoin
gentoo <- subset(penguins,species == "Gentoo")
adelie <- subset(penguins, species == "Adelie")
both <- rbind(gentoo, adelie)

### Note: this orders your data. Records from Gentoo will come first and then will come those of Adelie.

# using the subset function and the match operator
Gen_Ad <- subset(penguins, species %in% c("Adelie", "Gentoo"))

# simplified example for the match operator
x <-c("a", "b", "r", "f")
lookup <- c("a", "b")

# result is a logical vector, the length is the same as the length of 
x %in% lookup
# result c(TRUE, TRUE, FALSE, FALSE) 

# how to measure performance:
# the first argument is a block of code
system.time({
  for(i in 1:1000){
       mean(1:1000)
  }
})


# using the subset operators and the logical OR operator
Gen_Ad2 <- penguins[penguins$species == "Adelie" | penguins$species == "Gentoo", ]

# The same with which
# this is better for defense against missing values
Gen_Ad2 <- penguins[which(penguins$species == "Adelie" | penguins$species == "Gentoo"), ]

# using the dplyr package
library(dpylr)
Data.Peng1 <- penguins %>% filter(species == "Gentoo" | species == "Adelie")

# Operator objects can be accessed with the backticks
`%>%`
?`%in%`

# NA at the end does not change this result
which(c(penguins$species == "Adelie"
| penguins$species == "Gentoo", NA))

########################################----------------------------------------

# 2. With a loop for every species:
# a. subset the data by species
# b. create a directory for it
# c. write it into the new directory in rds and csv format!

# every file contains data that belong to one species.
# Two equally good ways to do this.
sp <- unique(penguins$species)
# sp <- levels(factor(penguins$species))

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


# Why we use which? 
# Create a different version of the penguins dataset,
# where we have a missing value in the species column
pen <- penguins 
pen$species[1] <- NA

# compare the subsetting with and
withWhich <- penguins[which(penguins$species == currentSpecies),]
# compare the subsetting without 
withoutWhich <- penguins[penguins$species == currentSpecies,]


# Contrasing .rds and the .RData
# readin in the rds file
adelie <- readRDS("Adelie/Adelie.rds")

# .RData: keeps the information of the names of the object
# saving as an RData file
save(adelie, file="Adelie/Adelie.RData")

# remove the file
rm(adelie)

# it is no longer there
ls()

# loading
load(file="Adelie/Adelie.RData")

# the object is back!
ls()

# not to do: Only gives you the names of the data
# o <- load(file="Adelie/Adelie.RData")

#######################################----------------------------------------

# 3. Change the previous code to save the subsets in a list. Use the species names as the names of the list!
# List can be heterogeneous, can contain different kinds of data 
li <- list(a=1:10, b=c("a", "b", "c"))

# how many elements are there
length(li)

# accessing elements:
# this is still a list
li[1]

# this is the element only: a vector
li[[1]]
li$a

# saving the results in an object.
# This is a container
container <- list()

# iterate through every species
for(i in 1:length(sp)){

    # The current species
    currentSpecies <- sp[i]

    # get the subset that contains all data of the species
    speciesData <- penguins[which(penguins$species == currentSpecies),]

    # write it into the directories
    container[[i]] <- speciesData

}

# the structure of the final list
str(container)
# good idea to indicate which element pertain to what species
names(container) <- sp


# acessing items in this list:
# the dollar operator
container$Chinstrap
# is essentially the same as this:
container[["Chinstrap"]]

# you can get deeper subsets of this: the first row of the data.frame
container[["Chinstrap"]][1, ]



# 4. Write a function that returns parts of the dataset based on column name and value,
# either greater or equal/lower or lower or equal to the value. function(x, column, value, greater)
#' @param x The data frame to subset
#' @param column Which column you are basing the subsetting on. (Numeric column)
#' @param value A numeric value, that separates the subset that you want to select from what you don't want.
#' @param greater Should a subset be selected that correponds to rows where the values are
#' greater than \code{value} (TRUE), or lower tan equal (FALSE).
