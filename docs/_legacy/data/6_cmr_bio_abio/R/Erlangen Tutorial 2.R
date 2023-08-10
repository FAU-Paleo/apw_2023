#Lee Hsiang Liow Tutorial for Erlangen
#Tutorial 2

#This tutorial introduces capture recapture data structure and *one* of the software available
#namely openCR. I will use openCR to illustrate a few things using an ecological dataset available 
#in openCR this also one of the standard examples used in many other different capture recapture software
#(which may have different "conventions" and functions), but if you understand this one, it is easier to
# use "new" ones. I am using openCR for this class because we can run everything I want to show you in 
#R (mac or windows) but be aware there are many other software that may be also useful or even more 
#appropriate/convenient for you such as MARK (with RMark), or you might even want to write your own in R 
#or Python.

#openCR and its "ancestor" serc (which you also need to have installed) are spatial capture-recapture 
#software, but the author (Murray Efford) has included non-spatial, which are the ones we will focus on.
#We will also only focus on "open models"

#openCR documentation and examples (which hopefully this course will help you "translate" in paleo questions)
#https://cran.r-project.org/web/packages/openCR/vignettes/openCR-vignette.pdf
#https://cran.r-project.org/web/packages/openCR/openCR.pdf
#https://www.otago.ac.nz/density/pdfs/openCR-examples.pdf

#To get the most of of this tutorial, look up definitions or parameters that are unclear or unknown to
#you as you go through stuff and read up in the links that come with this tutorial and ask questions. The
#code takes less than a minute to run through :-) It's easy to run capture recapture analyses ! But it takes
#thinking to figure out if the estimates make sense and if the models make sense etc.

rm(list=ls()) #start with a clean slate
library(openCR) # load openCR, this requires serc

#Let's first explore a dataset called dipperCH, so you can start thinking about the elements and structure of 
#this data and the similarity to paleo data, where individual animals are equivalent to individual taxa
#dipperCH is a capture recapture history (hence "CH) of the European
#Dipper (*Cinclus cinclus*) compiled in some locality from the years 1981–1987.This is a formally published dataset

#DipperCH in openCR is single-session secr capthist object. In other software, it may be formatted differently

#Q(v)Explore and describe the data verbally and discuss how it corresponds to paleo taxon data
# e.g. what is a "session"?

dipperCH #the data as a single-session secr capthist object
str(dipperCH) #explore how the data is structured

dipperCH[1:dim(dipperCH)[1],,] #the capture history (or encounter/detection history); each row is an individual bird
#rows will eventually be our individual taxa; and each column is the "time", which in this case is
#years, see attr(x = dipperCH, which = "sessionlabels") 
#or
unclass(dipperCH)[,,1]

#lets look at the attributes in turn
attr(x = dipperCH, which = "covariates") #these are discrete covariates associated with each individual
#you can imagine that if you are looking at say marine invertebrates, this might be infaunal versus epifaunal
#or if you are looking a vertebrate, this might be herbivore, omnivore and carnivore
#Q(v) Discuss with your desk partner what types of covariates might be interested in, we can make a class list during discussion
attr(x = dipperCH, which = "intervals") # ignore this, this is for more complex analyses with primary and secondary sessions 
#e.g. for robust design analyses see http://www.phidot.org/software/mark/docs/book/pdf/chap15.pdf
attr(x = dipperCH, which = "sessionlabels") #so these tell you which year the data are collected in

#summaries of the data, 
#Q(v)study n, R, n, r, and z
JS.counts(dipperCH)
?JS.counts

#ok let's fit some first models using this dataset so we don't have to bother with making a caphist
#object yet from real paleo data, but let's imagine how this data can be analogous to paleo data
#Q(v)Here are four different models below, what do each of them "do"?
#what do the things specified in the function openCR.fit mean?
dipper.const <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~1, p~1))
dipper.phi.t <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~1))
dipper.p.t <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~1, p~t))
dipper.phi.t.p.t <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~t))

#you got a warning messag fitting dipper.phi.t.p.t, check to see if you know
#what it means and what it means for inference. Is it ok to fit models wherer variance cal
#fails? What does it say about your model and your data ? How might you proceed (no one "correct)
#answer

#compare models using AIC
AIC(dipper.const,dipper.phi.t, dipper.p.t, dipper.phi.t.p.t)
#Q(v) which is the best model according to AIC?
#Q(v)what are the parameters of this model and how do the estimated parameters compare with those from the other ones?
#Q(v)Can you think of other models to fit (just a thought experiment), given what you know from looking at this data?

#let's plot the time-varying model and the model where both p and phi are constant through time
par(mfrow=c(2,1))
plot(dipper.phi.t.p.t, par = 'phi', ylim = c(0,1), pch = 16)
abline(h=predict(dipper.const)$phi[2,2],lwd=2) #this line added from the model with both p and phi constant
abline(h=predict(dipper.const)$phi[2,4],lwd=1, lty=2) #lower 95% CI for phi
abline(h=predict(dipper.const)$phi[2,5],lwd=1, lty=2) #upper 95% CI for phi

plot(dipper.phi.t.p.t, par = 'p', ylim = c(0,1), pch = 16)
abline(h=predict(dipper.const)$p[2,2],lwd=2) #this line added from the model with both p and phi constant
abline(h=predict(dipper.const)$p[2,4],lwd=1, lty=2) #lower 95% CI for p
abline(h=predict(dipper.const)$p[2,5],lwd=1, lty=2) #upper 95% CI for p
#Q(v) discuss why the phi and p plots are "offset" in time
#Q(c) can you plot the death (= extinction if working with paleodata) 
#probability instead of phi (survival probability?)

# if you feel like you don't really know what the best model *really* is you can do model averaging, although
#in this case, the dipper.constant model does "distinguish itself".
#Q(v) Discuss what model-averaging does
#Q(c) if you feel up to it, you can try some other ways that do not fully take into account sampling to count up the 
#number dead and compare to what we infer using these models. e.g. assume that the last time you see
#the animal is the last time that animal is alive, and then calculate some probability of death using this dataset
modelAverage(dipper.const , dipper.phi.t, dipper.p.t , dipper.phi.t.p.t )
#make sure you understand the output. Why is the first line of "p" all NAs? Why is the
#last row of "phi" all NAs?

#Let's now think about some other models that one might fit as a paleobiologist 
#while imagining that this dipper dataset is one of some marine invertebrate group over some
#geological time interval of which there are 7 equal "chunks"

#model with fixed time values
#e think of the sessions as time intervals, e.g. Paleozoic versus post-Paleozoic
sesscov <- data.frame(t.period = c(1,1,1,0,0,0,0)) #dummy values for splitting up the time periods
sampling.time <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~t.period), 
                                       sessioncov = sesscov)
survivorship.time  <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t.period, p~t), 
                                       sessioncov = sesscov)
both  <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t.period, p~t.period), 
               sessioncov = sesscov)
AIC(sampling.time, survivorship.time, both)

#plot the model called both -  did it do what you expected?
par(mfrow=c(2,1))
plot(both, par = 'phi', ylim = c(0,1), pch = 16)
plot(both, par = 'p', ylim = c(0,1), pch = 16)


#model with individual covariates (as factors): replace "sex" with say bivalve versus brachiopod in your mind!
surivorship.different <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~sex+t, p~t), 
                                            sessioncov = sesscov)
preservation.different <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~sex+t), 
                                    sessioncov = sesscov)
both.different <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~sex+t, p~sex+t), 
                                     sessioncov = sesscov)
more.complex <- openCR.fit(dipperCH, type = 'CJS', model = c(phi~sex*t, p~1), 
                             sessioncov = sesscov)
AIC(surivorship.different, preservation.different, both.different,more.complex)


#plot the best covariate model, i.e. preservation.different
par(mfrow=c(2,1))
plot(preservation.different, par = 'phi', ylim = c(0,1), pch = 16)#this shows for "male", 
plot(preservation.different, par = 'p', ylim = c(0,1), pch = 16)
# in preservation.different sampling is different for males and females. 
#How do we plot (or find) the value for females?
#Aha! it's through the link function beta estimates
#read http://www.phidot.org/software/mark/docs/book/pdf/chap6.pdf and see lecture notes

#this is what is being plotted (but only for males)
predict(preservation.different)

#test it
plot(preservation.different, par = 'p', ylim = c(0,1), pch = 16)#this shows for "male", 
points( seq(1, 6,1)+0.05,
        predict(preservation.different)$p$estimate[2:7])#just showing this so you know where the values come from

#let me first make a smaller data.frame of the numbers we need for me to demonstrate
est=as.data.frame( preservation.different$fit$estimate) #make a data frame of the things directly estimates (the betas)
rownames(est)= preservation.different$betanames #give them the names that correct so it is easy to see
est

#in the preservation.different model, we are saying that females have a different sampling probability, but
#the temporal dynamics are the same (as the males). In other words,there is temporal variability in the
#sampling of the population of males and females, but females as systematically different from males 
#(but this does not "interact" with time).... to model this, we have to specify a model that is 
#openCR.fit(dipperCH, type = 'CJS', model = c(phi~1, p~sex*t), sessioncov = sesscov))

#now I will show you how to "convert" or transform the beta values (see lecture) for the "standard" link function
#used which is a logit function (more about that in lectures) : what we aer using is something like this:
# logit (sampling) = beta1 + beta2(sex)+ b3(timeinterval)
# and the betas are what are in the data.frame we made, "est"
# so say I am interested in males in 1982 for this model
# then I need inverse of logit(1.0438558  + 0 + 0) = invlogit(1.043855) in R speak
# the sex is the dummy variable 0 and so is the first time interval 
# check it
predict(preservation.different)$p
#right?
# or next one
# I am interested in males in 1983 for this model
# the I need inverse logit(1.0438558  + 0 +1.6782084) = invlogit(1.0438558  + 0 +1.6782084) in R speak
#and so on! So this just a linear model.

#now we can calculate the female values and plot against the males:
female.p=preservation.different$fit$estimate[1]+preservation.different$fit$estimate[2]
for (i in 1:5){
female.p[i+1]=preservation.different$fit$estimate[1]+preservation.different$fit$estimate[2]+preservation.different$fit$estimate[i+2]
}
invlogit(female.p)

plot(preservation.different, par = 'p', ylim = c(0,1), pch = 16)#this shows for "male", 
points( seq(1, 6,1)+0.05,
        predict(preservation.different)$p$estimate[2:7])#just showing this so you know where the values come from
points( seq(1, 6,1)+0.05,
        invlogit(female.p), col="blue")#just showing this so you know where the values come from
#Q(v)why aren't the values the same distance from the males for each year?


#edge effects
###confounded parameters and identifiability see lecture notes
#you got a Warning message that reads:
##variance calculation failed for some beta parameters; confounding likely

#fixing values (or dropping these in onward analyses)

openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~t), details = list(fixedbeta=c(rep(NA,5),1e3,rep(NA,6))))
openCR.fit(dipperCH, type = 'CJS', model = c(phi~t, p~t), details = list(fixedbeta=c(rep(NA,11),1e3)))


#messsage from Murray Efford : I wouldn't call this an 'error' -
#just one of the joys of open population analysis. Non-identifiability happens. 
#More constructively - try another parameterization of xxx - y is often more 
#robust than z, and you can always feed that model as starting values to a z model.


#Q(v) we ran the standard CJS model, which is a closed population model that I used for illustration
#but we don't normally use this model in the standard diversification rates modeling 
#why? What does the CJS model assume?

#run (some) of the other models in Table 3 of the openCR-vignette (choose one type of model, e.g.
#c(phi~1, p~1)) and compare the output and Table 3. Why are the AIC's different ? 
#why are the p and phi different, even though we specified  c(phi~1, p~1))? (clue, look at the other parameters)

dipper.const.JSSAbCL <- openCR.fit(dipperCH, type = 'JSSAbCL', model = c(phi~1, p~1))
dipper.const.JSSAfCL <- openCR.fit(dipperCH, type = 'JSSAfCL', model = c(phi~1, p~1))
dipper.const.JSSAgCL <- openCR.fit(dipperCH, type = 'JSSAgCL', model = c(phi~1, p~1))
dipper.const.JSSAlCL <- openCR.fit(dipperCH, type = 'JSSAlCL', model = c(phi~1, p~1))

AIC(dipper.const.JSSAbCL,dipper.const.JSSAfCL, dipper.const.JSSAgCL,dipper.const.JSSAlCL )

#think a little about model comparison using AIC, what is being evaluated? Should you always let AIC "decide"
#which models are best to use for any given dataset? or should you decide (why and how?)

# which are the closed versus open formulations? what's the difference?
dipper.const.JSSAf <- openCR.fit(dipperCH, type = 'JSSAf', model = c(phi~1, p~1))
dipper.const.JSSAfCL <- openCR.fit(dipperCH, type = 'JSSAfCL', model = c(phi~1, p~1))
AIC(dipper.const, dipper.const.JSSAf, dipper.const.JSSAfCL)

#compare the open formulations
dipper.const.JSSAb <- openCR.fit(dipperCH, type = 'JSSAb', model = c(phi~1, p~1))
dipper.const.JSSAf <- openCR.fit(dipperCH, type = 'JSSAf', model = c(phi~1, p~1))
dipper.const.JSSAg <- openCR.fit(dipperCH, type = 'JSSAg', model = c(phi~1, p~1))
dipper.const.JSSAl <- openCR.fit(dipperCH, type = 'JSSAl', model = c(phi~1, p~1))

AIC(dipper.const.JSSAb,dipper.const.JSSAf, dipper.const.JSSAg,dipper.const.JSSAl )

#compare all of them
AIC(dipper.const , dipper.const.JSSAbCL,dipper.const.JSSAfCL, dipper.const.JSSAgCL,dipper.const.JSSAlCL, 
    dipper.const.JSSAb,dipper.const.JSSAf, dipper.const.JSSAg,dipper.const.JSSAl )

#I hope you are realizing that the models that can be built are many! We barely scratched the surface, but
#this is where your (paleo)biological knowledge and intuition comes to play! What models are useful to build?
#what do we learn?


####

#note that in many paleo applications, it's hard to have the sampling occasions short (an assumption built in) 
#or even in time as the geologic stage increases in amount of time, there will be a greater probability of speciation 
#or extinction or sampling! So all things being equal, if a "session" is longer than other, than the probability for 
# anything to happen will increase right?
#similarly (but not identically) as the "amount of time" increases, there might be more rock/sediment available
#there are different ways to deal with this as we will see in tutorial 4, but if you are using Mark/Rmark there is also 
#an easy way to incorporate. I will show you in the next tutorial how to overcome this issue by converting the probability
#to a rate after analyses in openCR, but note there are multiple ways to do this.

###END


