#Lee Hsiang Liow Tutorial for Erlangen
#Tutorial 1

#Applies to all tutorials - please let me know if something doesn't "work". Or if there is something
#that is weird or you don't understand. It will help me improve my notes!

#It is important to look at raw data in many different ways and to think about 
#potential issues and biases that may affect general or specific inferences. 
#Below are some questions and suggestions to get us thinking.
#If I mark the questions with (v), they can be discussed (with a desk partner or in your own head), 
#if they are marked with (c) they should be explored with code.
#After you have worked through the questions I have suggested, you can come up with different ones (that could have bearing on 
#diversification or richness analyses using paleo data) that you can share with the class.

rm(list=ls()) #start with a clean slate
#Q(v) How are PBDB data obtained and compiled?
#Q(v) What are some of the potential differences among taxa, sampling environments, biogeographic areas, 
#and time intervals and their "interactions", when considering PBDB data?

#Download canidae data from PBDB
library(paleobioDB)
#check this for more specifics https://paleobiodb.org/data1.1/occs/list for the function pbdb_occurrences
canidae=pbdb_occurrences( limit="all", vocab= "pbdb", base_name="Canidae", show=c("phylo", "ident", "late_interval"))

#Q(c) how many canidae genera are represented in pbdb and what are they? 
unique(canidae$genus)[order(unique(canidae$genus))] #good to check names, oops there is an "NA" beware!
length(unique(canidae$genus)) #this length is not right! the "NA" should be removed
which(is.na(canidae$genus)) #check, oh, quite a lot of entries are NAs 
canidae[which(is.na(canidae$genus)),] #oops beware of these, should be probably removed before (some) analyses!
good_genus_only=canidae[which(!is.na(canidae$genus)),] #we will use this data set that is a bit "cleaner"

#Q(c)how many species are represented?
#remove those with "sp." and "indet."  first
good_sp_only=canidae[which(canidae$species_name != "sp." & canidae$species_name != "indet."),]
unique(good_sp_only$taxon_name)[order(unique(good_sp_only$taxon_name))] # list of names, check if you like this
length(unique(good_sp_only$taxon_name)) #there are 447 species if you are ok with this

#plot the observed temporal ranges of the genera
genus.list=unique(good_genus_only$genus)[order(unique(good_genus_only$genus))] 

data=matrix(data=NA, nrow=length(genus.list), ncol=3) #lets make an empty matrix to put data in
rownames(data)=genus.list
colnames(data)=c("FIRST.APP", "LAST.APP", "N")
data[,1]=as.numeric(paste(tapply(good_genus_only$early_age, as.factor(good_genus_only$genus), max )))
data[,2]=as.numeric(paste(tapply(good_genus_only$late_age, as.factor(good_genus_only$genus), min )))
data[,3]=as.numeric(paste(tapply(good_genus_only$genus, as.factor(good_genus_only$genus), length)))
data=as.data.frame(data)
data#check if it looks nice

#plot the observed temporal ranges of the genera
plot(0, 0, ylim=c(0,length(genus.list)), xlim=c(-max(data[,1]), 10), type="n", xlab="Mya", ylab="genus")

for (i in 1:length(genus.list)){
  segments(-data[i,1], i, -data[i,2], i) #draw observed stratigraphic ranges
  text(8, i,  genus.list[i], cex = .5)#add genus names
  text(-55, i,  data[i,3], cex = .5) #add number of observations
}

#Q(c) What arethe mean and median number of observations for genera and species for the 
#time interval in which they are observed?
#for genera, since we already have that nice table
median(data[,3]) #16
mean(data[,3]) #48.4
#and for species 
median(as.numeric(paste(tapply(good_sp_only$taxon_name, as.factor(good_sp_only$taxon_name), length)))) #1
mean(as.numeric(paste(tapply(good_sp_only$taxon_name, as.factor(good_sp_only$taxon_name), length)))) #5.6

#Q(c) What are the mean and median number of observations in temporal bins you define for genera and species?
#ok let's decide on that first, let's see what the underlying data look like first
#there is something called "late_interval" in pbdb, but I am going to ignore that for now, as many cells are empty 
#so lets just focus on the early_interval, whose range brackets are given (in early_age and late_age)

unique(canidae$early_interval)#there are quite a lot of different systems and names in PBDB
#I am going to try equal sized time bins to make things easy, let's use the "clean" genus data
#table the number of observations per 5 million years for, say the genus
good_genus_only$mid_age=0.5*(as.numeric(good_genus_only$early_age)+as.numeric(good_genus_only$late_age))
breaks=seq(-60,0, 5 )
good_genus_only$bins=cut((-good_genus_only$mid_age), breaks=breaks) 
head(good_genus_only)#see if the conversions look right
good_genus_only$bins_n=good_genus_only$bins

x=as.character(good_genus_only$bins_n)
for (i in 1:12){
x=replace(x, x==as.character(levels(good_genus_only$bins)[i]), 13-i)
}
good_genus_only$bins_n=as.numeric(x) 
table(good_genus_only$genus, good_genus_only$bins_n)

#table the number of observations per 1 million years for, say the genus and compare this with what you see for the 5 million year 
#binning - Q(v) what thoughts come to mind? (notice the "missing" bins, just means there was no data)
breaks=seq(-60,0, 1 )
good_genus_only$bins=cut((-good_genus_only$mid_age), breaks=breaks) 
head(good_genus_only)#see if the conversions look right
good_genus_only$bins_n=good_genus_only$bins

x=as.character(good_genus_only$bins_n)
for (i in 1:60){
  x=replace(x, x==as.character(levels(good_genus_only$bins)[i]), 61-i)
}
good_genus_only$bins_n=as.numeric(x) 
table(good_genus_only$genus, good_genus_only$bins_n)

approx.time=-as.numeric(colnames(table(good_genus_only$genus, good_genus_only$bins_n)))+0.5
obs=colSums(table(good_genus_only$genus, good_genus_only$bins_n))     
plot(approx.time,obs, pch=20)
points(approx.time,obs,type="l")
for (i in 1: 61){
points(approx.time, table(good_genus_only$genus, good_genus_only$bins_n)[i,], type="l", lwd=0.5, col=i+1)
}


#Q(v) how gappy is the canidae record as represented in the PBDB? Do you think the function xxx represents what is close 
#to the "true" range of candidae? do you think if we sampled intensively if we can fill the temporal gaps?
#what about space?

#Q(c) Choose any taxon you are interested in yourself (genus level upwards) and repeat these exercises/questions. 
#Q(v) Do the distributions/numbers looks similar or different from canidae? why?



