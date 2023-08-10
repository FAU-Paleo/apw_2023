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
#   02: Exploring patterns of sampling
# 
# *************************************************


## Package(s) used in this script:
library(tidyverse) # for data manipulation functions and ggplot
library(viridis) # for colour scales


# Sampling proxy counts ---------------------------------------------------


## Let's explore sampling patterns!
## First we'll calculate counts of sampling proxies and plot these alongside raw diversity

# Taxa per interval 
count_taxa <- vector("numeric") # create empty vector for the loop below to populate
for (i in 1:nrow(intervals)) { # for-loop to count each taxon that appears in each interval
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma) # uses our intervals dataframe
  count_taxa[i] <- (length(unique(out$accepted_name)))
  print(count_taxa[i])
}

# Collections per interval
count_colls <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma)
  count_colls[i] <- (length(unique(out$collection_no)))
  print(count_colls[i])
}

# Formations per interval
count_formations <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(occ_data, max_ma > intervals[i,]$min_ma & min_ma < intervals[i,]$max_ma)
  count_formations[i] <- (length(unique(out$formation)))
  print(count_formations[i])
}


## For equal-area gird cells, I would recommend the package 'dggridR'. However, it is unfortuately not available
## for newer versions of R. So, we would need to load an older version of dggridR (version 2.0.3)
## To do this in your own time, you can uncomment the following lines and follow the instructions in the link

#library(devtools) # load functions from devtools
#install_version("dggridR", version = "2.0.3", repos = "http://cran.us.r-project.org") 
# ## When prompted, enter: 3 (i.e. 'None') - you don't want to update any other packages at this time
#library(dggridR) # load dggridR
#install_version("dggridR", version = "2.0.3", repos = "http://cran.us.r-project.org")
# ## Check out the instructions for the package here: https://github.com/r-barnes/dggridR


## Gather the proxy information together in a new dataframe for plotting:
proxy_counts <- data.frame(intervals$interval_name, intervals$mid_ma, count_taxa, count_colls, count_formations)
## Rename the columns for ease:
proxy_counts <- rename(proxy_counts, 
                       "interval_name" = "intervals.interval_name", 
                       "mid_ma" = "intervals.mid_ma")

## Finally, convert all zero's to NAs for plotting 
## This means that the plots won't register zero and instead will leave gaps 
## where there is no sampling instead - this gives a more realistic picture
proxy_counts[proxy_counts == 0] <- NA 



# Sampling plots ----------------------------------------------------------

## Let's get plotting these patterns so we can see what they look like!

## Plotting using ggplot

## Set interval boundaries for the dotted lines on the plot
## We'll also use this vector again, so its handy to have :)
int_boundaries <- c(rev(intervals$max_ma), 174.1)

## Set up your ggplot layers (first 'layer' goes on the bottom, then the next 'layer', and so on):
proxy_plot <- ggplot() + 
  # Formations (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", size = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_formations), colour = "orangered3", size = 4, shape = 16) +
  # Collections (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", size = 1.2, linetype = "dashed")  +
  geom_point(data = proxy_counts, aes(mid_ma, count_colls), colour = "peru", size = 5, shape = 16) +
  # Taxa (as dots and a line):
  geom_line(data = proxy_counts, aes(mid_ma, count_taxa), colour = 'black', size = 1.2)  +
  geom_point(data = proxy_counts, aes(mid_ma, count_taxa), colour = "black", size = 4, shape = 16) +
  # Add a minimal theme - but you can make your own custom themes too!
  theme_minimal() + 
  labs(x = "Time (Ma)", y = "Sampling proxy counts") +
  # Make sure to reverse the x-axis to match geological time!
  scale_x_reverse(breaks = int_boundaries) +
  # And tidy up our y-axis with even breaks that match the totals in our dataframe:
  scale_y_continuous(breaks = seq(0, 320, 20))
## Call the finished plot to the RStudio plots tab:
proxy_plot

## Set dimensions and save plot (as pdf) to the plots folder
#dir.create("./plots") # create new folder if have you haven't created on manually 
ggsave(plot = proxy_plot,
       width = 20, height = 15, dpi = 500, units = "cm", 
       filename = "./plots/sampling_proxies.pdf", useDingbats=FALSE)



# Collections per latitude ------------------------------------------------


## Now let's creates a plot to see where the collections across time and paleolatitude
## We'll also colour our plot according to the number of taxa in each collection
## This is also a visualisation of alpha diversity or 'local richness'

## There is evidence to suggest that alpha diversity is not as strongly affected by sampling biases
## as gamma (or 'global') diversity. For a more sophisticated way to calculate alpha diversity by
## treating taxonomically indeterminate occurrences as valid, see the method described in 
## Close et al. (2019) - code available here: https://github.com/emmadunne/local_richness

## Let's get our data set up:
lat_data <- occ_data # rename object to keep the original separate

## Create new column for mid_ma
lat_data$mid_ma <- (lat_data$max_ma + lat_data$min_ma)/2 

## Next, we'll need to count the number of taxa per collection (i.e. their frequency):
taxa_freqs <- count(lat_data, collection_no)

## Subset lat_data to only the columns we need:
lat_data <- lat_data %>% 
  select(collection_no, paleolat, paleolng, mid_ma) %>% 
  distinct() %>% na.omit()

## Add add the frequency information:
lat_data <- left_join(taxa_freqs, lat_data, by = "collection_no")

## Before we plot, let's order the frequencies and remove any NAs that have crept in:
lat_data <- lat_data %>% arrange(n) %>% na.omit()

## Take a look:
View(lat_data)


## Set up our ggplot layers
lat_plot <- ggplot(data = lat_data, aes(x = mid_ma, y = paleolat, colour = n)) +
  geom_vline(xintercept = int_boundaries, lty = 2, col = "grey90") +
  geom_hline(yintercept = 0, colour = "grey10") +
  scale_color_viridis(trans = "log", breaks = c(1, 2, 12), direction = -1, option = "D") + # set the break= to match your richness data
  #scale_y_continuous(labels = function(x) format(x, width = 5), limits = c(-70, 70), breaks = seq(from = -60, to = 60, by = 20)) +
  scale_x_reverse(breaks = int_boundaries) + 
  theme_minimal() + 
  theme(legend.position = "none", legend.direction = "vertical", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Palaeolatitude (º)") +
  geom_point(size = 4, alpha = 0.5) # (alpha sets point transparency)
lat_plot # call to plot window


## Set dimensions and save plot (as pdf)
ggsave(plot = lat_plot,
       width = 20, height = 10, dpi = 500, units = "cm", 
       filename = "./plots/lat_locations.pdf", useDingbats=FALSE)

