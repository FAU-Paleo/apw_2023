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
#   03: Simple stats to look at sampling
# 
# *************************************************

## Load package(s):
library(tidyverse) # for data manipulation functions and ggplot
library(vegan) # for diversity stats


# Regression --------------------------------------------------------------

## Let's do some simple regression plots to see how raw diversity correlates
## with sampling proxies:

## Raw diversity vs. collections
ggplot(proxy_counts, aes(x=count_taxa, y=count_colls)) + 
  geom_point(shape=17, size = 6, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1")  +
  theme_minimal()

## Raw diversity vs. formations
ggplot(proxy_counts, aes(x=count_taxa, y=count_formations)) + 
  geom_point(shape=16, size = 5, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1") +
  theme_minimal()


## Let's quantify these relationships:
lm_colls = lm(count_colls ~ count_taxa, proxy_counts) # linear model
summary(lm_colls) # summary of results

lm_forms = lm(count_formations ~ count_taxa, proxy_counts)
summary(lm_forms)




# Collector's curve -------------------------------------------------------

## First, get our data into shape

## Table the number of each species per collection
abun_data <- table(occ_data$collection_no, occ_data$accepted_name) # abundance
#abun_data[abun_data > 0] <- 1 # if we want a presence/absence table, we can add this step

## Turn this into a matrix
abun_matrix <- matrix(abun_data, ncol = length(unique(occ_data$accepted_name)))
colnames(abun_matrix) <- colnames(abun_data) # add the column names back in for when we need to check anything
rownames(abun_matrix) <- rownames(abun_data) # same for the row names

## Using the vegan R package, we can make a species accumulation curve
## Check out the help file for specaccum() to find out more about the methods
sp_accum <- specaccum(abun_matrix, method = "collector")
sp_accum <- specaccum(abun_matrix, method = "random")
summary(sp_accum) # bring up the summary

## Plot the curve in base R - you can make this pretty in ggplot if you prefer!
plot(sp_accum, ci.type = "poly", col = "chartreuse4", lwd = 2, ci.lty = 0, ci.col = "chartreuse")

