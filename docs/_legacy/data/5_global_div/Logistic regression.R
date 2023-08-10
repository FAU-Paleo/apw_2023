
rm(list=ls())

##Packages ----
library(divDyn) #loads of tools for paleobiological analysis
library(lme4) #glms
library(modEvA) # a nice package to evaluate glms
library(car) #logit transformation
library(metafor)

###Some logistic regression basics ----

#    glm(formula, family = binomial,data)

#The link function for the logistic model is the logit function

(x<-seq(0,1,0.01))

#Logit() transforms binomial data into a normal distribution (minus infinity to infinity)
plot(x,logit(x)) #note the warning to avoid infinite outputs
hist(logit(x))
#Coin flips are TRUE/FALSE or 1 and 0 in reality but probability scales 0 to 1
#This probability scale is where our extinction risk predictions fall e.g. 0.3 = 30% chance of extinction

#in one time interval 100 genera get one trial each at 50:50 survive or not, random odds
(x<- rbinom(100, 1, 0.5)) 
hist(logit(x))# in reality a species either survives (0) or goes extinct (1)
#But as you've seen above, we can also give non-integer logit probabilities e.g. 0.5
logit(0.5)

#https://www.geo.fu-berlin.de/en/v/soga/Basics-of-statistics/Logistic-Regression/The-Logit-Function/index.html

#To find out where logit comes from, we need to introduce the odds-ratio
#If the probability of an event is a 0.5, or 50:50, the odds are one-to-one or even
0.5/(1-0.5) #=1:1 or 1/1
#If the probability is 1/3, the odds are one-to-two
(1/3)/(1-(1/3)) #=1:2 or 1/2

#Log-odds, the logarithm of the odds ratio
#Using the above probability values
log(0.5/(1-0.5)) #1:1 becomes 0
log((1/3)/(1-(1/3))) #1:2
#Negative logits represent probabilities below one half 
#Positive logits correspond to probabilities above one half.

#To return from log-odds to odds-ratio use exp()
exp(log(1))
exp(log(0.5))
#Then to probability, prob = odds / (1 + odds)


##Simulating a coin toss ----

#To identify a weighted coin vs unweighted coin (two groups_)
#100 observations
coindat<-data.frame(results=NA,coin=c(rep("a",100),rep("b",100)))
coindat[,1]<- c(rbinom(100, 1, 0.5),rbinom(100, 1, 0.6))

m1<-glm(results ~ coin, family="binomial",data=coindat) 
summary(m1)

#10000 observations
coindat<-data.frame(results=NA,coin=c(rep("a",10000),rep("b",10000)))
coindat[,1]<- c(rbinom(10000, 1, 0.5),rbinom(10000, 1, 0.6))

m1<-glm(results ~ coin, family="binomial",data=coindat) 
summary(m1)

#prob = odds / (1 + odds)
#Intercept should be ~0.5
exp(summary(m1)$coefficients[,1])/(1+exp(summary(m1)$coefficients[,1]))


##Extinction determined by geographical range and its trajectory? ----

#The following R code is amended from the publication Kiessling, W., & Kocsis, Á. T. (2016). Adding fossil occupancy trajectories to the assessment of modern extinction risk. Biology letters, 12(10), 20150813.
#https://doi.org/10.1098/rsbl.2015.0813

##Load in datasets 
#Set the working directory to where the files from the repository are saved
setwd("C:/Users/...")
setwd("C:/Users/reddin/Documents/Teaching")

# Load data of occupancy history of extinct species of the last 10 myr
dat.10 <- read.table("dat10.csv", sep=";", header=T)

str(dat.10) #species names are row names
cor(dat.10[,1:3],use = "pairwise.complete.obs")

# Test simple models
r.ocl <- glm(el ~ ocl, family=binomial(), data=dat.10) # check occupancy
summary(r.ocl)
r.chl <- glm(el ~ chl, family=binomial(), data=dat.10) # check change of occupancy
summary(r.chl)

#Can also use interactive terms and model selection
#occupancy (oc)
#change in occupancy (ch)
#change of occupancy over two intervals (ch2)	
#change of occupancy over three intervals (ch3)

# Create most simple and most complex models
# most complex with all interactions
lim.10 <- subset(dat.10, !is.na(chl3)) # set filter to longer time series (to have same sample size) 
model <-  glm(el ~ ocl * chl * chl2 *chl3, family=binomial(), data=lim.10) 
summary(model) 
#note that NAs reduce the complete dataset size to
nrow(dat.10) - sum(apply(dat.10[,c("ocl","chl","chl2","chl3")],1,anyNA))

# most simple model (everything extinct regardless of variables)
nothing <- glm(el ~ 1, family=binomial, data=lim.10)
summary(nothing)

# Model selection
forwards = step(nothing, scope=list(lower=formula(nothing), upper=formula(model)), direction="forward")
summary(forwards)
Dsquared(forwards)  # Check deviance explained my model

AIC(nothing,forwards,model) #lowest AIC is most parsimonious

plotGLM(forwards)
termplot(forwards,se=T,partial.resid=T) #cannot plot interactions but shows odds of extinction (on logit scale)
#are higher when occupany is low or change is negative

plot(forwards)



####Extinction selectivity by clades and traits ----

#The following R code is amended from the publication Reddin et al. 2019, Marine clade sensitivities to climate change conform across time scales
#https://doi.org/10.1038/s41558-020-0690-7

#Preprepared fossil data (from online or local)
#slc.dat<-read.csv("Reddin et al. Fossil selectivity supporting data.csv") 
slc.dat<-read.csv("https://zenodo.org/record/3465281/files/Reddin%20et%20al.%20Fossil%20selectivity%20supporting%20data.csv?download=1")
#or
slc.dat<-read.csv("Reddin et al. Fossil selectivity supporting data.csv",header=T)

data(stages) #from divDyn

##Fill in body size data from Payne & Heim 2020
gen1 <- read.csv("https://datadryad.org/stash/downloads/file_stream/495584", header=T)
#or
gen1 <- read.csv("genus.sizes.ranges.rev.csv", header=T)

#Run the next few lines as one
#amend one different class name
unique(gen1$class);table(gen1$class) 
levels(gen1$class)<-c(levels(gen1$class),"Actinopterygii")
gen1$class[gen1$class=="bony fish"]<-paste("Actinopterygii")

#Make clgen for Payne & Heim data
gen1$clgen<- paste(gen1$class, gen1$genus)

sum(duplicated(gen1$clgen))#Remove duplicates
gen1<-gen1[!duplicated(gen1$clgen),]

#Fill in body size as a new column for our dataset
slc.dat <-merge(slc.dat,gen1[,c("clgen","logvol")],by="clgen",all.x=T)
names(slc.dat)[names(slc.dat)=="logvol"]<-paste("size")
sum(is.na(slc.dat$size)==F)/length(slc.dat$size) #42% of the genera are filled


###Investigating extinction selectivity ----

#True time intervals are expected to have non-random or directional odds. 
#Usually, most organisms survive, with some more likely than others

hyperther<-c(51,58,61,74,76,83) #Bin numbers of the hyperthermal onsets: PT, TJ, PLiToa, Apt-Alb, Cen-Tur, PETM
stages[hyperther,]

i=83 #PETM = extinctions observed in Selandian-Thanetian stages
#Note, in this dataset, 1 or TRUE is extinction

#the intercept is the mean of logit(response), or the extinction intensity (no extinction selectivity yet)
m0<-glm(ext ~ 1, family="binomial",data=subset(slc.dat,slc.dat$Slc==i)) 
summary(m0)

#We can check this:
mean(subset(slc.dat,slc.dat$Slc==i)$ext) #binomial mean
logit(mean(subset(slc.dat,slc.dat$Slc==i)$ext)) #Compare against regression intercept


#size is a numerical variable so intercept remains partial mean (i.e. the mean extinction intensity after accounting for variation associated with size)
#extinction intensity variation associated with size is extinction selectivity (i.e. extinction rate of one group relative to another)
#We hypothesise a size selective extinction pressure at the PETM
m1<-glm(ext ~ size, family="binomial",data=subset(slc.dat,slc.dat$Slc==i)) 
summary(m1) #not significant? If it was, negative slope would mean larger size genera more likely to have response = 0, or less likely to go extinct

#Remember it's a good idea to check the diagnostic plots (though they look different to usual datasets)
plot(m1)

#What about higher clades and categorical regression?
m2<-glm(ext ~ Taxa, family="binomial",data=subset(slc.dat,slc.dat$Slc==i)) 
summary(m2)
#In this case the intercept becomes one of the groups, taken as a reference point, 
#The intercept coefficient is the extinction intensity of this group with a p-value testing if it is significantly difference from zero
#R uses the first group in alphabetical order as a reference
table(subset(slc.dat,slc.dat$Slc==i)$Taxa) #Clades in alphabetical order
#Means Actinopterygii is used as the intercept
#All other coefficients are relative to the reference group, their p-values are testing if 
#their extinction intensity is different to the reference extinction intensity i.e. extinction selectivity
table(slc.dat[slc.dat$Slc==i & slc.dat$Taxa=="Retaria",]$Class) # Forams removed from 'Retaria'
#=Radiolarians significantly more likely to go extinct than fish

stages[hyperther,]

#Not many extinctions at PETM. Try end-Permian
i=51 #end-Permian = extinctions observed in Changhsingian stage
m1<-glm(ext ~ Taxa, family="binomial",data=subset(slc.dat,slc.dat$Slc==i)) 
summary(m1) 

#Extinction selectivity by trait
m2<-glm(ext~motility,family="binomial",data=subset(slc.dat,slc.dat$Slc==i))
summary(m2)


##Averaging coefficients ----
#How does this hyperthermal extinction selectivity (intensity variation among groups) compare with other times?
#Using -1 removes the intercept term, which gives easy access* to the expected mean value and standard errors per group but DO NOT do this for a numerical slope as it forces the slope to go through the origin (unless that makes theoretical sense)
#Discussion on this practice here: https://stats.stackexchange.com/questions/7948/when-is-it-ok-to-remove-the-intercept-in-a-linear-regression-model
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==51)))$coefficients[,1]
#vs
earlier_coeff<-cbind(
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==50)))$coefficients[,1],
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==49)))$coefficients[,1],
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==48)))$coefficients[,1],
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==47)))$coefficients[,1])
apply(earlier_coeff,1,mean) #simple averaging
#The difference between PT coefficients and their previous coefficients by motility groups
summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==51)))$coefficients[,1]-apply(earlier_coeff,1,mean)

#Collect the standard errors
earlier_err<-cbind(
  summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==50)))$coefficients[,2],
  summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==49)))$coefficients[,2],
  summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==48)))$coefficients[,2],
  summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==47)))$coefficients[,2])

#Can use meta-analysis, rma(), to give a weighted average of the extinction coefficients
m.all.active<-rma(earlier_coeff[row.names(earlier_coeff)=="motilityactively mobile",],sei=earlier_err[row.names(earlier_err)=="motilityactively mobile",],method="REML")
forest(m.all.active,main="actively mobile",slab=stages$stage[50:47],xlab="Extiction risk, log-odds") #Look at contribution of each time bin
#Dashed line to mark the PT coefficient
abline(v=summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==51)))$coefficients["motilityactively mobile",1],lty="dashed")

m.all.statio<-rma(earlier_coeff[row.names(earlier_coeff)=="motilitystationary",],sei=earlier_err[row.names(earlier_err)=="motilitystationary",],method="REML")
forest(m.all.statio,main="stationary",slab=stages$stage[50:47],xlab="Extiction risk, log-odds")
#Dashed line to mark the PT coefficient
abline(v=summary(glm(ext~motility-1,family="binomial",data=subset(slc.dat,slc.dat$Slc==51)))$coefficients["motilitystationary",1],lty="dashed")


#Using random effects to account for nested factors ----
#e.g. motility groups are nested within clades. Fish, for example are nearly always actively mobile. Corals are nearly always stationary
#we can account for this nested information using random effects
#what is the effect of motility after accounting for the motility associated with clade differences?

x<-subset(slc.dat,slc.dat$Slc==i)
(betterTax<-table(x$Taxa)[table(x$Taxa)>10]) #remove the least well-represented clades

m2<-glmer(ext~motility+ (motility | Taxa),family="binomial",data=subset(slc.dat,slc.dat$Slc==i & slc.dat$Taxa %in% names(betterTax) ))
summary(m2) #look at the fixed effects
table(slc.dat[slc.dat$Slc==i & slc.dat$motility=="passively mobile",]$Class)

#Try this for size selectivity at the Olenekian
x<-subset(slc.dat,slc.dat$Slc==53)
(betterTax<-table(x$Taxa)[table(x$Taxa)>5]) #remove the least well-represented clades

m3<-glmer(ext~size+ (size | Taxa),family="binomial",data=subset(slc.dat,slc.dat$Slc==53 & slc.dat$Taxa %in% names(betterTax)))
summary(m3)
#From Payne & Heim 2020, "selective loss of large-bodied animals is the exception, rather than the rule, in the evolution of marine animals"
#Singularity "leads to random-effect variance estimates of (nearly) zero, or estimates of correlations that are (almost) exactly -1 or 1"



##Have a play? Choose a stage, choose a trait or clade grouping ----
#but first, tidying up this dataset (sorry, I've mixed code and datasets a bit above!)
slc.dat$diet[slc.dat$diet %in% c("grazer, carnivore","grazer, suspension feeder","grazer, photosymbiotic","browser, omnivore","grazer, browser")]<-paste("grazer")
slc.dat$composition[slc.dat$composition %in% c("low Mg calcite, agglutinated")]<-paste("low Mg calcite")
slc.dat$composition[slc.dat$composition %in% c("chitin","agglutinated","\"sclero-protein\"","no hard parts")]<-NA
slc.dat$life_habit[slc.dat$life_habit %in% c("epifaunal, solitary, polymorph, depth=surface")]<-paste("epifaunal")
slc.dat$life_habit[slc.dat$life_habit %in% c("boring","gregarious")]<-NA

stages[hyperther,] ##e.g. hyperthermals listed in Foster et al. 2018
names(slc.dat)

i=58 #Choose your favorite (but don't go p-fishing, decide a hypothesis first)
m2<-glm(ext~composition,family="binomial",data=subset(slc.dat,slc.dat$Slc==i))
summary(m2)




########Prob extinctions ----
#The code for converting a extinctions matrix (e.g. built by function divDyn::modeltab() ) to a 
#probabilistic extinctions matrix is provided at the bottom of the code but needs the original dataset 
#to obtain information

#needs spdat (the time-binned PBDB data) with a column called Taxa (the clades you are interested in)
sam3t<-divDyn(spdat,"clgen","slc")$samp3t #Best do this before removing data e.g. singletons

taxnams<-c("Actinopterygii","Ammonoidea","Anomalodesmata","Annelida","Arcida","Arthropoda","Asteroidea",     
           "Cardiida","Carditida","Chondrichthyes","Coleoidea","Crinoidea","Echinoidea","Gastropoda","Gymnolaemata", "Hippuritida",
           "Holothuroidea","Lucinida", "Megalodontida", "Modiomorphida","Myalinida","Mytilida","Nautiloidea","Nuculanida", "Nuculida","Ophiuroidea","Ostreida",
           "Pectinida","Pholadida","Porifera","Rhynchonelliformea" ,"Scaphopoda",
           "Scleractinia","Solemyida","Stenolaemata","Trigoniida","Venerida") 
taxnams<-unique(spdat$Taxa)

taxnamsall<-matrix(NA,ncol=2,nrow = length(taxnams))
for (i in 1:length(taxnams)){
  if(length(unique(spdat[spdat[,"Taxa"]==taxnams[i],"slc"]))>2){
    taxsam3t<-divDyn(spdat[spdat$Taxa==taxnams[i],],"clgen","slc")$samp3t
    taxnamsall[i,]<-c(taxnams[i],median(taxsam3t[!taxsam3t%in%c(1,0)],na.rm=T))# remove 1s, 0s and NAs for averaging
  }else{taxnamsall[i,1]<-taxnams[i]}}
taxnamsall[order(as.numeric(taxnamsall[,2])),]

#cycle through each of your time bins, here 43:94
subsdat2<-NA
for (i in 43:94){ 
  subsdat<-subset(slc.dat,slc.dat$Slc==i)
  
  #This uses the sampling completeness as the probability of true extinction.   
  for (j in 1:length(taxnams)){
    subsdat$ext[subsdat$Taxa==taxnamsall[j,1]]<-subsdat$ext[subsdat$Taxa==taxnamsall[j,1]] * ((sam3t[i+1] + as.numeric(taxnamsall[j,2])+0.09)/2)}  #+0.09 is the difference between the means of tax and time

  subsdat3<-subsdat #Store this as just the extinctions from the current bin, else the below line will keep accumulating the taxa as we move forward in time
  subsdat<-rbind(subsdat,subsdat2) #add subsdat2, probability remainder from previous bin, to subsdat
  
  #Do here whatever you want to do with these data – i.e. analyse or store them
  slc.dat2<-rbind(slc.dat2,subsdat)
  
  subsdat2<-subsdat3[subsdat3$ext>0,] #The remainder of the probablility of LAD (via sampling completeness) is then saved for the next bin 
  subsdat2$ext<- 1-subsdat2$ext 
  
  print(i)
}
slc.dat2<-slc.dat2[-1:-2,]
