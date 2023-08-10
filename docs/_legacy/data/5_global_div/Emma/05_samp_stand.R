# *************************************************
#
#   Analytical Paleobiology Workshop 2022
#
#   Module 5: Global diversity dynamic
#   Day 8 (AM): Sampling standardisation 
#
#   Tuesday, August 30th 2022
#   Emma Dunne (emma.dunne@fau.de)
# _________________________________________________
#
#   Coverage-based rarefaction in R
# 
# *************************************************

## This script will walk you through coverge-based rarefaction, also called
## Shareholder Quorum Subsampling (SQS) via the package iNEXT, following the 
## methods outlined in Dunne et al. (2018)


## Packages used in this script:
library(tidyverse)
library(iNEXT)


## Coverage-based rarefaction, or Shareholder Quorum Subsampling (SQS) uses 
## rank-order abundance to estimate diversity through subsampling at different 
## degrees of sampling coverage, or 'quorum levels', and estimates diversity 
## using a metric called Good's u. 

## We will use the package iNEXT, which estimates diversity using Hill numbers  
## via  subsampling with the equations of Chao and Jost (2012), and also implements
## extrapolation, using the Chao1 estimator - for more info, see the main citation:
citation("iNEXT")
## or the paper describing the package: 
## https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12613


## Start off by dowloading the occurrence data for early tetrapods using an API request
## Let's call the object tet_data
## And clean it like we did last week

## To off, here is the key info we need:
taxon_group <- "Tetrapoda"
start_interval <- "Bashkirian"
stop_interval <- "Kungurian"


# YOUR AMAZING CODE HERE ;)


## Once that's all cleaned up, pear it down to just the columns we need, while
## creating a new data object so that we don't overwrite the original:
taxa_data <- subset(tet_data, select=c(accepted_name, occurrence_no, collection_no,
                                        early_interval, late_interval, min_ma, max_ma))


## To get your data in the right shape for iNEXT, you'll need:
##   1. a dataframe of interval names
##   2. total number of occurrences in each interval (i.e. sampling units)
##   3. genus incidence frequencies for each interval


## If you're feeling fancy, you can grab these using the method we used last week,
## Or you can simply load the dataset below:
intervals <- read.csv("./data/ET_intervals.csv")


## 2 + 3. Get the genus incidence frequencies for each interval
## In the output, we need the first entry of each list object to be the total
## number of sampling units (e.g. occurrences, species, etc.), followed by the 
## incidence frequencies (i.e. how often a unique one occurs)
## For example: Interval_A : 150 99 96 80 74 68 60 54 46 45
## = there are 150 occurrences in Interval_A, 99 in the first collection, 
## 96 in the second, etc.

freq_data <- lapply(1:nrow(intervals), function(i) {
  tmp <- taxa_data %>% filter(max_ma >= intervals[i,"max_ma"] & min_ma <= intervals[i,"min_ma"]) %>% 
    count(., accepted_name) %>% arrange(desc(n)) %>% 
    add_row(n = sum(.$n), .before = 1) %>%
    select(n)
  freq_raw <- as.numeric(tmp$n)
  freq_raw
})
names(freq_data) <- intervals$interval_name # give each list element its correct interval name
str(freq_data) # call up the list to see a summary of its structure
freq_data[[6]] # check one of the intervals


## Now that we've got our data in the right shape, let's do some estimates!

## Create a vector of quorum levels that we want to compute
## 0.4 is considered the 'standard', but the fashion now is to plot multiple quorum levels
quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)

## And create a new list to store output
estD_output <- list() 

## estimateD() is the main function for running SQS in iNEXT

for(i in 1:length(quorum_levels)) {
  # estimateD function in iNEXT, running over each quorum level
  estD_tmp <- estimateD(freq_data, datatype = "incidence_freq", base = "coverage", level = quorum_levels[i])
  # filter to the diversity estimates (order = 1):
  estD_tmp <- filter(estD_tmp, order == 0)
  # organise the output:
  estD_tmp$quorum_level <- quorum_levels[i]
  estD_tmp$mid_ma <- intervals$mid_ma
  # add the output to the newly created list
  estD_output[[i]] <- estD_tmp
}

## The output is a 'list' object so we'll need to convert it into a dataframe and clean it up before plotting
estD_plotting <- bind_rows(estD_output) # binds rows of a list

## Ensure that the quorum level column is being treated as a 'factor' to avoid errors while plotting:
estD_plotting$quorum_level <- as.factor(estD_plotting$quorum_level)

## Create a colour gradient for as many colours as you have quorum levels:
teal_gradient <- scales::seq_gradient_pal("turquoise", "darkslategrey", "Lab")(seq(0, 1, length.out = 4))

## Set interval boundaries for the plot's axis
int_boundaries <- c(rev(intervals$max_ma), 272.2)

iNEXT_plot <- ggplot(estD_plotting, aes(x = mid_ma, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  ## Each quorum level is called individually to be plotted:
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.3), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[1], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.4), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[2], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.5), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[3], alpha = 0.2) +
  geom_ribbon(data = subset(estD_plotting, quorum_level == 0.6), aes(x = mid_ma, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = teal_gradient[4], alpha = 0.2) +
  ## Set our line and point sizes (and shapes):
  geom_line(size = 1) +
  geom_point(aes(pch = method), size = 4.5) +
  scale_shape_manual(values=c(15, 16, 17)) +
  ## Add our colours, theme, and axes labels:
  scale_colour_manual(values = teal_gradient) +
  scale_x_reverse(breaks = int_boundaries) +
  labs(x = "Time (Ma)", y = "Coverage rarified richness") +
  theme_minimal()
iNEXT_plot # Call the plot to the plots tab


## Save a copy of the plot to the plots folder
ggsave("./plots/iNEXT_plot.pdf", plot = iNEXT_plot, 
       width = 30, height = 18, units = "cm")
