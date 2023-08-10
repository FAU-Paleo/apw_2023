# *************************************************
#
#   Analytical Paleobiology Workshop 2022
#
#   Module 2: Paleodiversity analyses
#   Day 4: Exploring fossil occurrence data  
#             and fossil record biases
#
#   Thursday, August 25th 2022
#   Emma Dunne (emma.dunne@fau.de)
# _________________________________________________
#
#   01: Data preparation
# 
# *************************************************


## Package(s) used in this script:
library(tidyverse) # for data manipulation functions and ggplot

## If you've been using a lot of different packages, some function names might be duplicated.
## This step ensures that the function 'select' is coming from the dplyr package (part of tidyverse)
select <- dplyr::select




# Getting PBDB data -------------------------------------------------------


## Choose a taxonomic group and time interval
## (Note that the scripts for today are set up for the interval from the Late Triassic - Early Jurassic, so you might prefer to chose a taxon group in this interval)
taxon_group <- "Pseudosuchia" # Taxon group
start_interval <- "Carnian" # Interval to start at
stop_interval <- "Toarcian" # Interval to stop at


## Create an API request form the Paleobiology Database and store this URL as an object
## A list of API options can be found here: https://paleobiodb.org/data1.2/
URL <- paste0("https://paleobiodb.org/data1.2/occs/list.csv?base_name=", # occurrence data, as a .csv
              taxon_group, "&interval=", start_interval, ",", stop_interval, # use our inputs from above
              "&show=full&pres=regular") # any additional columns we want 

## Then use this to load the data into R:
occ_data_raw <- as_tibble(read.csv(URL, header = TRUE, stringsAsFactors = FALSE))

## Take a peep:
glimpse(occ_data_raw) # view columns
View(occ_data_raw) # open as new tab


# Cleaning the data -------------------------------------------------------


## Before we do any plotting or analyses, let's clean our data up a bit

## Remove 'super-generic' identifications, so that we only retain occurrences to species- and genus-level
occ_data_raw2 <- filter(occ_data_raw, (identified_rank %in% c("species","genus")))

## Remove occurrences with “aff.”, “ex. gr.”, “sensu lato”, “informal”, or quotation marks
occ_data_raw3 <- occ_data_raw2 %>% filter(!grepl("cf\\.|aff\\.|\\?|ex\\. gr\\.|sensu lato|informal|\\\"", identified_name)) 

## We've already filtered our data to 'form taxa' (i.e. body fossils) in the API query, but just in case any trace fossils crept through...
## Remove entries marked as 'trace' or 'soft', and those with no genus name
occ_data_raw4 <- occ_data_raw3[occ_data_raw3$pres_mode != "trace", ] # trace taxa
occ_data_raw5 <- occ_data_raw4[!grepl("soft",occ_data_raw4$pres_mode), ] # 'soft' preservation
occ_data_raw6 <- occ_data_raw5[occ_data_raw5$genus != "", ] # missing genus name

## Finally, filter the data so any duplicate taxon names or collection numbers are eliminated:
occ_data <- distinct(occ_data_raw6, accepted_name, collection_no, .keep_all = TRUE)

## Take a look:
## How much has our data been reduced by?
glimpse(occ_data)


## Sometimes, with PBDB data, you'll find some annoying trace terms still remain after you've 
## cleaned your data through these steps. Some researchers (myself included) like to keep lists
## of pesky taxa so we can remove them more explicitly. You'll find code for this here:
## https://github.com/emmadunne/pbdb_cleaning


## Save a copy as a .csv file for safe keeping
#dir.create("./datasets") # create new folder if one doesn't exist
write_csv(occ_data, "./data/occs_cleaned_pseudo.csv")




# Set up time interval info -----------------------------------------------


## Next, let's grab information for time intervals from the PBDB so
## that their ages in Ma match the same ages as the occurrences

## Download names and ages of time intervals from the PBDB:
intervals_all <- read.csv("http://paleobiodb.org/data1.1/intervals/list.txt?scale=all&limit=all")

## Make a vector of stage names that we will specifically be looking at 
## (we'll be using this vector later for plotting too)
interval_names <- c("Carnian", "Norian", "Rhaetian", # Late Triassic
                    "Hettangian", "Sinemurian", "Pliensbachian", "Toarcian") # Early Jurassic

## Select these intervals from the full PBDB intervals tibble:
intervals <- filter(intervals_all, interval_name %in% interval_names)

## Pear this down to just the columns we'll need:
intervals <- dplyr::select(intervals, interval_name, early_age, late_age)

## For ease of use later, let's rename the age columns to match the occurrence data:
intervals <- rename(intervals, "max_ma" = "early_age", "min_ma" = "late_age")

## And finally, calculate the midpoint for each interval
intervals$mid_ma <- (intervals$min_ma + intervals$max_ma)/2

## Take a peep:
View(intervals) # open as new tab

## Save a copy as a .csv file
write_csv(intervals, "./data/intervals_Car_Tor.csv")

