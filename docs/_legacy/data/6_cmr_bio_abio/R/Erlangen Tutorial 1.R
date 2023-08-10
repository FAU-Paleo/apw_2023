#Lee Hsiang Tutorial for Erlangen Analytical Paleobiology Workshop 1.9.2022
#Tutorial 1

#It is important to look at raw data in many different ways and to think about 
#potential issues and biases that may affect general or specific inferences. 
#Below are some questions and suggestions to get us thinking.
#If I mark the questions with (v), they can be discussed (with a desk partner or in your own head), 
#if they are marked with (c) they should be explored with code.
#After you have worked through the questions I have suggested, you can come up with different ones (that could have bearing on 
#diversification or richness analyses using paleo data) that you can share with the class.
#Suggested solutions are given in "Erlangen Tutorial 1 - with suggested answers.R" you can
#look at that whenever, but I suggest working on your "own" for half an hour first.

rm(list=ls()) #start with a clean slate
#check out https://paleobiodb.org/#/ if you haven't already.
#Q(v) How are PBDB data obtained and compiled?
#Q(v) What are some of the potential differences among taxa, sampling environments, biogeographic areas, 
#and time intervals and their "interactions", when considering PBDB data?

#Download canidae data from PBDB
library(paleobioDB)
#check this for more specifics https://paleobiodb.org/data1.1/occs/list for the function pbdb_occurrences
canidae=pbdb_occurrences( limit="all", vocab= "pbdb", base_name="Canidae", show=c("phylo", "ident", "late_interval"))

#Q(c)how many canidae genera are represented in pbdb and what are they? 

#Q(c)how many species are represented?

#Q(c) plot the observed temporal ranges of the genera

#Q(c) plot the observed temporal ranges of the species


#Q(c) what are the mean and median number of observations for genera and species for the 
#time interval in which they are observed? and for species

#Q(c) what are the mean and median number of observations in temporal bins (try 5 and 10 million year bins) you define for genera and species?

#Q(c) table the number of observations per 1 million years for, say the genus and compare this with what you see for the 5 million year 
#binning - Q(v) what thoughts come to mind? (notice the "missing" bins, just means there was no data)

#Q(c) plot the temporal ranges of the species if you have time

#Q(v) how "gappy" is the canidae record as represented in the PBDB? Do you think the function xxx represents what is close 
#to the "true" range of candidae? do you think if we sampled intensively if we can fill the temporal gaps?
# What about space?

#Q(c) Choose any taxon you are interested in yourself (genus level upwards) and repeat these exercises/questions. 
#Q(v) Do the distributions/numbers looks similar or different from canidae? why?



