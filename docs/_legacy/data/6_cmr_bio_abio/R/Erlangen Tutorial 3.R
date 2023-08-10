#Lee Hsiang Liow Tutorial for Erlangen
#Tutorial 3

#This tutorial has a bit of data massaging before and after for PBDB like data for openCR (making serc objects)
#analyses. We can use the data from Tutorial 1 downloaded from PBDB
#for some realism, rather than bird observation data
#When you are done with this or anything else you want to try with your own data/quetsions
#you can explore the R script from Lidgard et al. https://royalsocietypublishing.org/doi/full/10.1098/rspb.2021.1632

#useful for input data for openCR (historically coming from serc)
#https://www.otago.ac.nz/density/pdfs/secr-datainput.pdf
library(openCR) 
library(paleobioDB) #we need functions from paleobioBD to connect with the paleobiolog database
#check this for more specifics https://paleobiodb.org/data1.1/occs/list for the function pbdb_occurrences
data=pbdb_occurrences( limit="all", vocab= "pbdb", interval= "Cenozoic", base_name=c("Cheilostomata"), show=c("phylo", "ident", "late_interval"))
#getting data using this function took a couple of minutes for me so don't worry if you have to wait a bt
#N =5151
#work with data that is id-ed to genus at least
length(which(data$taxon_rank=="species" | data$taxon_rank=="genus")) #N=5028, lost a few
data=data[(which(data$taxon_rank=="species" | data$taxon_rank=="genus")),]

# we will work on genus level, let's check everything looks ok; check you know what the columns mean and that we
#are extracting what we want https://paleobiodb.org/data1.1/occs/list
unique(data$genus)[order(unique(data$genus))]

#let's reduce the dataset a bit
bryo=subset(data, select=c("family", "genus", "species_name","early_interval", "early_age", "late_age"))
head(bryo) #ah.. much nicer
length(unique(bryo$family)) #117 family names
length(unique(bryo$genus)) #416 genera
length(unique(paste(bryo$genus, bryo$species_name)) )
#there are 1862 species names there are going to be some "sp." etc, but now
#we don't bother to do any cleaning, since I plan to do analyses at genus level for now. But think about this one.
       
#let's explore the data a bit first and notice you might have to convert the data type
str(bryo)

bryo$family=as.factor(bryo$family)       
bryo$genus=as.factor(bryo$genus)     
bryo$early_interval=as.factor(bryo$early_interval) 
bryo$early_age=as.numeric(bryo$early_age) 
bryo$late_age=as.numeric(bryo$late_age) 
bryo$mid_age=(bryo$early_age+bryo$late_age)*0.5      
str(bryo) # better

# I am going to randomly add a covariate at genus level to this data so we can play with setting up different models
#imagine it to be e.g. encrusting versus erect cheilostome bryozoans, which can be real data, but I
# don't have time to check each genus for now :-)
genera=as.data.frame(unique(bryo$genus))
colnames(genera)=c("genera")
genera$covariate=runif(dim(genera)[1],min=0,max=1)
genera$covariate=replace(genera$covariate, which(genera$covariate>0.5), 1)
genera$covariate=replace(genera$covariate, which(genera$covariate<=0.5), 0)
one=genera$genera[which(genera$covariate==1)]

bryo$covariate=0
bryo$covariate=replace(bryo$covariate, which(bryo$genus %in% one), 1)

#quick a dirty look at the temporal distribution of the data we downloaded: 
hist(-bryo$mid) 
#you see that there is the "pull of the recent". There are way fewer records  before 20 mya. 
#is this "real" maybe there were really fewer cheilostomes between 50-40 my??

#let's assign the data to standard stages, get the data for stages in first.
stages=c("Danian","Selandian","Thanetian",
"Ypresian", "Lutetian", "Bartonian",  "Priabonian",
"Rupelian", "Chattian", "Aquitanian", "Burdigalian","Langhian", "Serravallian",
"Tortonian","Messinian","Zanclean","Piacenzian","Gelasian","Calabrian", "Almost.now")
max=c(66, 61.6,  59.2, 56, 47.8, 41.3, 38,33.9,28.1, 23.03,  20.44,
      15.97, 13.82, 11.62,7.246, 5.333, 3.6, 2.58,1.8,  0.781)
min=c(61.6, 59.2,56, 47.8, 41.3, 38, 33.9,  28.1, 23.03, 20.44,
       15.97,13.82, 11.62, 7.246, 5.333, 3.6, 2.58, 1.8, 0.781, 0)
num=seq(1,20,1)
stages=as.data.frame(cbind(stages, max, min, num))
stages$stages=as.factor(stages$stages)
stages$max=as.numeric(stages$max)
stages$min=as.numeric(stages$min)
stages$mid=-(stages$max+stages$min)*0.5
stages$num=as.numeric(stages$num)
breaks=c(-(stages$max),0 )
stages$bins = cut(stages$mid, breaks=breaks) 
stages$dur=stages$max-stages$min

bryo$bins = cut(-bryo$mid_age, breaks=breaks)

bryo=merge(bryo,stages, by="bins")

#tidy up
bryo=subset(bryo, select=c("family", "genus", "species_name", "stages","num", "mid",
                           "bins", "covariate"))

#make openCR format
trapXY.1=trapXY[1,] #all traps the same number, trapXY is from Open.CR, just lazy copying 
#but because we are not using spatial data, it doesn't matter. It's just that these cannot be missing
dummytraps <- read.traps(data = trapXY.1)

#now we will set it up like openCR wants it, it wants Sessions and IDs and trapID's minimally and the Occasion
#has the temporal sequence

detection=as.data.frame(rep(1, dim(bryo)[1]))
colnames(detection)=c("Session") #one session
detection$ID=bryo$genus
detection$Occasion=bryo$num#ordered bin
detection$trapID=rep(1, dim(bryo)[1])
detection$covariate=bryo$covariate

dnd.matrix=suppressWarnings(make.capthist (detection, dummytraps, fmt = "trapID", covnames="covariate"))
Erect=subset(dnd.matrix, covariates(dnd.matrix)$covariate==1) #dirty
Encrust=subset(dnd.matrix, covariates(dnd.matrix)$covariate==1)


JS.direct.est=JS.direct(dnd.matrix)
?JS.counts #will exaplan the parameters you see 
?JS.direct #will exaplan the parameters you see 
#e.g. n = detected, R released, m= previously marked, r = detected later, z alive but not marked
#phi is survival from one interval to the next, while p is sampling probability within interval

#converting from from se (e.g. sep, sepN etc) to CI
#The standard approach to calculating 95% confidence limits for some parameter theta is theta ±
#(1.96 × SE). However, to guarantee that the calculated 95% CI is [0, 1] bounded for parameters 
#(like phi or p) that are [0, 1] bounded,we can calculate the 95% CI on the logit scale, 
#before back-transforming to the real probability scale. However, because the logit transform is 
#not linear, the reconstituted 95% CI will not be symmetrical around the parameter estimate, 
#especially for parameters estimated near the [0, 1] boundaries.

#Note: JS.direct spits out the estimates on the "real" scale

par(mfrow=c(1,1), mar=c(5,4,4,2))
plot(stages$mid, JS.direct.est$N, pch=20, cex=1, axes=F, ylab="N", xlab="MYA",
     main="", cex.main=1, ylim=c(0,600))
points(stages$mid, JS.direct.est$N, type="l")
arrows(stages$mid, JS.direct.est$N+(1.96*JS.direct.est$seN), stages$mid, JS.direct.est$N-(1.96*JS.direct.est$seN), code=3, angle=90, length=0.1)
#grey lines are n + z which is range through 
points(stages$mid, JS.direct.est$n+JS.direct.est$z, pch=20, cex=1, col="darkgrey")
points(stages$mid, JS.direct.est$n+JS.direct.est$z, type="l", col="darkgrey")
axis(1)
axis(2)
box()

#Figure you just plotted above: it is a bit "odd" that it declines toward the Recent for both, 
#but that will be "remedied" by add e.g. genera we know are extant from WORMS, our own data and adding worms show a different picture
#sorry for a bit of self-advertising but you can check 
#https://royalsocietypublishing.org/doi/full/10.1098/rspb.2021.1632 fig. 1
# also, this tells us that we have thinking to do about estimating sampling (and the scope of the data)
#easier to discuss in class, remind me!

#many species in the Ypresian (4th dot from left), but maybe it is a longer stage? so maybe not so surprising? but 
#interesting this does not show when you do rangethroughs

plot(stages$mid,JS.direct.est$p, type="l", ylim=c(0,1), lty=2, axes=F, lwd=2, ylab="sampling", xlab="Mya")
points(stages$mid, JS.direct.est$p, type="l")
arrows(stages$mid, JS.direct.est$p+(1.96*JS.direct.est$sep), stages$mid, JS.direct.est$p-(1.96*JS.direct.est$sep), code=3, angle=90, length=0.1)
axis(1)
axis(2)
box()

#blue is "corrected" for the time interval, i.e. the instantaneous" sampling in that time interval
#so it's more comparable.
p=-log(1-JS.direct.est$p)/stages$dur #have to drop the penultimate one
pl=-log(1-(JS.direct.est$p-(1.96*JS.direct.est$sep)))/stages$dur #have to drop the penultimate one
pu=-log(1-(JS.direct.est$p+(1.96*JS.direct.est$sep)))/stages$dur #have to drop the penultimate one
points(stages$mid, p, type="l", col="blue", lty=2)
arrows(stages$mid, pl, stages$mid, pu, code=3, angle=90, length=0.1,col="blue",)

#JSSAl = POPAN

#Exercise, fit some models we tried in tutorial 2 and maybe some others? Do you expect the randomly 
#assigned covariates to have any bearing on sampling or some "biological" rate you image?

#try some other datasets and have fun and ask questions!



#these will take a while and wil give some warning messages which we can discuss
full.Pradelg <- openCR.fit(dnd.matrix, type='Pradelg', model= list(p~t, phi~t, gamma~t)) #=JSSAgCL
full.Pradel <- openCR.fit(dnd.matrix, type='Pradel', model= list(p~t, phi~t, lambda~t))#== JSSAlCL
full.JSSAg <- openCR.fit(dnd.matrix, type = 'JSSAg', model = c(phi~t, p~t, gamma~t))
full.JSSAl <- openCR.fit(dnd.matrix, type = 'JSSAl', model = c(phi~t, p~t, lambda~t))
full.JSSAgCL <- openCR.fit(dnd.matrix, type = 'JSSAgCL', model = c(phi~t, p~t, gamma~t))
full.JSSAlCL <- openCR.fit(dnd.matrix, type = 'JSSAlCL', model = c(phi~t, p~t, lambda~t))
const.JSSAg <- openCR.fit(dnd.matrix, type = 'JSSAg', model = c(phi~1, p~1, gamma~1))
const.JSSAl <- openCR.fit(dnd.matrix, type = 'JSSAl', model = c(phi~1, p~1, lambda~1))

#compare full.JSSAg$fit$estimate and full.Pradelg$fit$estimate
#compare  full.JSSAl$fit$estimate and full.Pradel$fit$estimate
#The non-spatial model types ‘Pradel’ and ‘Pradelg’ are implemented
#in openCR using sufficient statistics (Pradel 1996) and therefore fall outside the main 
#framework (Table 3). They correspond to ‘JSSAlCL’ and ‘JSSAgCL’ respectively, and estimate 
#the same parameters as those models. Estimates should coincide except when there are losses 
#on capture. ‘Pradel’ is parameterized in terms of population growth rate (lambda) and 
#‘Pradelg’ is parameterized in terms of seniority (gamma)


#Exercises
#fit models that 
#Convert the probabilities to rates so they are comparable across time intervals
#fit some models using the covariates I created
#fit some time models, where for instance, you specify differences between the Paleogene and the Neogene


#Try some simulations (do your estimates give you data that look plausible? can you get the estimates back? 
#What' the ucnertainty? how much data do you need?)
## using type = 'Pradelg' and extractfn = derived
## homogeneous p
fitarg <- list(type = 'Pradelg', model = list(p~t, phi~t, gamma~t))
turnover <- list(phi = c(rep(0.85, 5), rep(0.5,5)), lambda = 1.0, recrmodel = 'discrete')
cores=1
outg <- runsim.nonspatial(nrepl = 1, N = 500, ncores = cores, turnover = turnover,
                          p = 0.2, recapfactor = 4, nsessions = 10, noccasions = 1,
                          fitargs = fitarg, extractfn = derived)
outg[[1]]
apply(sapply(outg, function(x) x$estimates$lambda),1,mean)
