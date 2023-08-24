# Example script to demonstrate some basic statistical ideas.
# Ádám T. Kocsis
# 2023-08-22
# Erlangen, APW 2023

################################################################################
# A sample from the standard normal distribution:
# 100: sample size (n)
# 0: population mean (mu)
# 0: population standard deviation (sigma)
one <- rnorm(100, 0, 1)

# There are some irregularities of course
hist(one)

################################################################################
#1.  Estimating the population mean with the sample mean

# Our ability to do this depends on the sample size. The higher the sample size see the variance in the mean. 
trials <- 20

# display the sample mean with n=10, actual sampling of the population
for(i in 1:trials){
  one <- rnorm(10, 0, 1)
  message(mean(one))
}


# display the sample mean with n=100, actual sampling of the population
for(i in 1:trials){
  one <- rnorm(100, 0, 1)
  message(mean(one))
}

# display the sample mean with n=100, actual sampling of the population
for(i in 1:trials){
  one <- rnorm(200, 0, 1)
  message(mean(one))
}

# In order to assess this variance, we need to record the results

# number of iterations 
trials <- 1000

# create an empty container
container <- rep(NA, trials)

# A for loop
for(i in 1:trials){
  one <- rnorm(100, 0, 1)
  container[i] <- mean(one)
#   hist(one)
}

# how much does this vary? 
sd(container)

# This variance reflects our ability to estimate the true mean from the sample. This is so important that mathematical statisticians derived a way to estimate this standard devi ation analytically (with equations). It turns out that this can be approximated reasonably well, if we take the standard deviation of the sample, and divide it with square root of the sample size:
sd(one)/sqrt(length(one))

# this is the sampled based estimation of the Standard Error of the Mean. This is working very well if the data are normally distributed (i.e. Gaussian). But it might not work accuratetly if this cannot be assumed. 


################################################################################
# For these cases, we might want to try resampling, which does not assume anything about the distribution of the values - the only thing it assumes is that the sample is **representative** of the population, and has high large enough size that random combinations of the sample elements (with replacement) can be almost unique. Then we can use a simulation called **bootstrapping** to estimate the standard error. 

# Let's assume that this sample represents the population very well. 
one

# such resampling typically need high iteration size
trials <- 10000
bootContainer <- rep(NA, trials)

# for every bootstrapping trial
for(i in 1:trials){
	# create a replicate (resampled elements)
	# has the same size as the original sample
	# replacement is on: some elements will be present multiple times, creating a ~unique combination of elements
	bootstrapReplicate <- sample(one, length(one), replace=TRUE)

	# then we calculate from it the same thing as from the sample
	bootContainer[i] <- mean(bootstrapReplicate)
}

# container holds a sample from the 'resampling distribution'. Which will have a mean that is very close to the sample's actual mean, but not exactly the same.
mean(bootContainer) - mean(one)

# If you need to have the same exact mean, then you need to do a balanced-bootstrap, which means that overall, the sample elements need to be used during the generation of the replicates at an equal frequency. 

# The standard deviation of this resampling distribution sample will be your estimate of the standard error. 
sd(bootContainer)

################################################################################
# Comparison of sample central values

# In order to make the examples in the script reproduce exactly the same way on your computer, you can set a seed, which will ake the pseudorandom number generators to work in sync (!DO NOT USE THIS IN ANALYSES!)
set.seed(0)

# Let's assume that we have two beds where we have fossils coming from. We are interested in the body size measurements of these, which we will simulated with a single value. 

#  Let's assume that in reality, bed2 has slightly larger organisms (but not too much), compare the means 100 vs 103
# We take relatively small samples from these. 
bed1 <- c(rnorm(10, 100, 3))
bed2 <- c(rnorm(20, 103, 3))

# The means will also be slightly off from the population means because of the randomness
mean(bed1)
mean(bed2)

# If we can assume that the two samples come from normal distributions, we can use parametric statistical tests to compare these means. These are called  parametric, becuase we assume and implicitly estimate some parameters of the population's distribution, which enables us to derive some test conditions analytically. Since the normal distribution is among the simplest distributions and is very frequent, many parametric tests require normality. 

# are they normal?
hist(bed1)
hist(bed2)

# Maybe. They do look more or less normal, so let's try a two-sample t.test
# A two sample t.test(): parametric test on means
t.test(bed1, bed2)

# This test is giving you information on whether the difference between these two means is meaningful. Your p-value will give you the probability of making the mistake that the the difference is actually meaningless (null hypothesis), assuming that you accept that the difference is meaningful (alternative hypothesis). This is pretty high, way above the conventional 0.05, or 0.01 values, so it makes more sense to stay with the null hypothesis, so that we cannot distinguish the two sample from each other in terms of means.


# This can be confirmed with a Wilcoxon rank sum test which, is the non-parametric counterpart of the T-test. It compares medians, rather than means, which in the case of the normal distribution are very similar (Because we are working with a made-up example, we do not really know that the sample come from a normal distribution!) 

# wilcox.test() median
median(bed1)
median(bed2)

# executing the wilcoxon test:
wilcox.test(bed1, bed2)
# The argumentation is similar as above. The p-value will give you the chance that you incorrectly reject the null hypothesis that the true difference in medians (location shift) is just 0. Since this p-value is really high, this is quite likely, so we can stick with the null hypothesis, that the difference is meaningless.

# Typically these tests are less powerful than ideal, which is why some people also use resampling to do these tests.

# The idea is that we are explicitly simulating what kind of differences occurr between the samples if the two samples created from the same source (really there is no difference between them). If we see that the observed value occurrs quite frequently when the fossils are distributed between bed1 and bed2 randomly, then it is very likely that the difference is meaningless.

# We typically choose high sample sizes to do this. 

# resampling
trials <- 10000

# hypothetically merging the two populations
pseudopop <- c(bed1, bed2)
# container to store the results
results <- rep(NA, trials)
# iterate for every trial
for(i in 1:trials){
	# recreating sample structure for bed1 and bed2
	repBed1 <- sample(pseudopop, length(bed1), replace=TRUE)
	repBed2 <- sample(pseudopop, length(bed2), replace=TRUE)

	# record the randomly emerging difference in replicates
	results[i] <-  mean(repBed2)- mean(repBed1)
}

# the resulting sample from the resampling distribution is a hat-shaped distribution 
hist(results)

# with the mean of approximately 0
mean(results)

# Our results emerge from comparing this distribution to what we observed
observation <- mean(bed2) - mean(bed1)
observation

# showing that on the histogram
abline(v=observation, col="red")

# We can derive a percentile p-value from this. How frequenly do we see that the the randomly generated differences are as high as what we observed? - if this is very high, than it means that there is a decent chance that the difference that we see occurrs becasue of randomness
sum(results >= observation)/trials

# In our case, the true mean difference is just 3, which is in general too small to detect. But try increasing the mean of bed2! You will see how it is getting less and less likely to generate the observed difference randomly - making the detection of the mean difference easier and easier.

# If you are interested in resampling methods, I wholeheartedly recommend the course notes by Kowalewski and Novack-Gottshall (2010) Resampling methods in paleontology in Alroy and Hunt (eds.) Quantitative methods in paleobiology.  
