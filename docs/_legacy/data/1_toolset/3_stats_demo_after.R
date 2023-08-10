# Erlangen, 2022
# Ádám T. Kocsis
# CC-BY (attribution)
# 0. Download and load the file."bed1_length_sample.rds"
setwd("/media/adam/work/Dropbox/Teaching/Workshops/2022-08-22_APW/teach/1_toolset/3_stats_demo/")
samp <- readRDS("data/bed1_length_sample.rds")

########################################----------------------------------------

# 1. Analyze and plot the distribution from the file.
# - Plot the distribution.
hist(samp)
# - Look at summary statistics
median(samp)
summary(samp)
quantile(samp, c(0.3, 0.5))
# - Do a hypothesis test.
shapiro.test(samp)
ks.test(samp, "pnorm", mean(samp), sd(samp))

########################################----------------------------------------

# 2. Load "bed1_length.RData". This make-belief example represents our population.
# - Generate a single random sample of 30 specimens!
load("data/bed1_length.RData")
len <- length
set.seed(1)
sample(len, 30)
sample(len, 30)

########################################----------------------------------------

# 3. Create 3 differently sized samples containing 20, 50 and 100 elements! Bind them together in a list!
# - name the list, so the names contain how many elements the list has.
# - calculate the mean and standard deviation of the samples and put the results in a 3 by 2 matrix. 
samp20 <- sample(len,  20)
samp50 <- sample(len,  50)
samp100 <- sample(len,  100)

# adding the data to a list
samples <- list(samp20, samp50, samp100)
names(samples) <- c(20, 50, 100)

# results
mean_sd <- matrix(NA, ncol=2, nrow=3)
colnames(mean_sd) <- c("mean", "sd")
rownames(mean_sd) <- names(samples)

# manual fill: use character-subscript when possible!
mean_sd["20","mean"] <- mean(samples[["20"]])
mean_sd["50","mean"] <- mean(samples[["50"]])
mean_sd["100","mean"] <- mean(samples[["100"]])

mean_sd["20","sd"] <- sd(samples[["20"]])
mean_sd["50","sd"] <- sd(samples[["50"]])
mean_sd["100","sd"] <- sd(samples[["100"]])


########################################----------------------------------------

# 4. Repeat 3. for every sample size between 5 and 400 with a for loop.
# - find the highest and lowest realized value across the samples!
# - calculate the means in the samples, store them in a named vector (names are sample sizes). Plot the means as a function of sample size! 
# - calculate the minimum and maximum (range of) values it every sample.
# Draw the maximum and minimum as a function of sample size! 
# sample sizes
sizes <- 5:400
# second take
# sizes <- 5:1000

# container
multiSamples <- list()

# for all sample size 
for(i in 1:length(sizes)){

	# one random sample
	oneSample <- sample(len,  sizes[i])

	# save it 
	multiSamples[[i]] <- oneSample

}
names(multiSamples) <- sizes


# iterating the mean calculation
means <- lapply(multiSamples, mean) # produces list
means <- sapply(multiSamples, mean) # produces vector

# extremes of the samples
ranges <- sapply(multiSamples, range)
rownames(ranges) <- c("min", "max")

# find the most extreme values for plotting
vect <- unlist(multiSamples)
range(vect)

# plotting usually looks ok when the limits are a bit beyond this
plot(sizes, means, ylim=c(10, 80))
lines(sizes, ranges[1,])
lines(sizes, ranges[2,])


########################################----------------------------------------

# 5. Add the width of the specimens, "bed1_width.RData".
# - Load the file!
# - Create a data.frame and match the width with the length measurements!
# - adjust the Collect function so it works, when you use data.frames as an input!
load("data/bed1_width.RData")

# inspect structure!
str(width)

# the specimen names do not match that of `len`

# reordering things based on the names is a relatively safe operation - even
# when the sets are different
namesLen1 <- c(names(len[1:5]), "asdfasdf")

# produces missign values when you look for a name that does not exist!
width[namesLen1]

# reorder to match len
widthOrdered <- width[names(len)]

# create a simple data.frame
df <- data.frame(length=len, width=widthOrdered)

# Defining a function: 

# step 1. define inputs (arguments)
# df: the data frame
# 20: the number of entities

# function_name <- function(df, 20){
# 
# }

# step 2. define output:
# a data frame: result
#
# function_name <- function(df, 20){
# 	# Some calculation
# 
#	# return
#   return(result)
# }

# step 3. write prototype code, how you derive the result from the input
#
# function_name <- function(df, 20){
#
#	# index of the rows (random 20 rows)
#	randomIndex <- sample(1:nrow(df), 20, replace=FALSE)
#
#	# actual subsetting
#	result <- df[randomIndex, ]
#
#	# return
#   return(result)
# }

# step 4. rename the arguments, generalize variables: df->x, 20->n
# FYI: 20 or other numbers cannot be arguments, arguments have to be variable names!
#
# function_name <- function(x, n){
#
#	# index of the rows (random n rows)
#	randomIndex <- sample(1:nrow(x), n, replace=FALSE)
#
#	# actual subsetting
#	result <- x[randomIndex, ]
#
#	# return
#   return(result)
# }

# step 5: Write documentation. The standard styl is Roxygen2, which is what we use in R packages.
#' Collect random specimens (rows) from a data.frame.
#'
#' Select a given number of rows randomly from a data.frame class object. 
#' 
#' @param x Data frame that you want to sample
#' @param n Target sample size
#' @return The function returns a data.frame with n (randomly selected rows)
Collect <- function(x, n){

	# index of the rows (random n rows)
	randomIndex <- sample(1:nrow(x), n, replace=FALSE)

	# actual subsetting
	result <- x[randomIndex, ]
	
	# return the data frame
	return(result)
}

# step 6: testing, testing, testing...
oneSample <- Collect(x=df, n=20)
oneSample <- Collect(x=df, n=30)

# corner cases!
oneSample <- Collect(x=df, n=0)
oneSample <- Collect(x=df, n=1)

# Defense
# Does the function fail when it should? 
# Is the error making sense? 
oneSample <- Collect(x=df, n=nrow(df)+1)

########################################----------------------------------------

# 6. Create a 20 element sample and:
# - draw a scatterplot
# - calculate Pearson's covariance and correlation coefficient

# a simple scatterplot
oneSample <- Collect(x=df, n=20)
plot(oneSample)

# covariance (Pearson's)
cov(oneSample$length, oneSample$width)
cov(oneSample) # check out the sign!

# correlation (Pearson's)
cor(oneSample$length, oneSample$width)
cor(oneSample)

# sample size senzitivity!
oneSample <- Collect(x=df, n=100)
plot(oneSample)
cor(oneSample$length, oneSample$width)
cor.test(oneSample$length, oneSample$width)

########################################----------------------------------------
# Linear models
# dependent ~ independent
# predicted ~ predictors
oneSample <- Collect(x=df, n=20)
lm(oneSample$width ~ oneSample$length)
mod <- lm(width ~ length, data=oneSample)

plot(oneSample)
abline(mod)

################################################################################
# This was not covered in class:
# 7. Following from 6. Create 1000 different samples with 20 elements. Calculate the correlation coefficient from them, check out the distribution!

# number of iterations
trials <- 1000

# container to hold results
samples <- list()

# iterate for every trials
for(i in 1:trials){
	# select sample and store it a matrix
	samples[[i]] <- Collect(x=df, n=20)
}

# coefficients
correlations <- sapply(samples, function(x) {
	cor(x$length, x$width)
})

# the correlations can vary quite much solely because of random sampling
hist(correlations)

# this is why it is important to test for significance!
cor.test(samples[[1]]$length, samples[[1]]$width)

# note that the distribution above is the actual sampling distribution based on the population!
# The stats of cor.test() function are based on the sample!

########################################----------------------------------------

# 8. Make a linear model for one of these samples. Also fit a 2nd and a 10th order polynomial. Which is the best model and why?    
set.seed(2)
oneSample <- Collect(x=df, n=40)
plot(oneSample, pch=16, col="black")
linear <- lm(width ~ length, data=oneSample)
second <- lm(width ~ poly(length, 2), data=oneSample)
tenth <- lm(width ~ poly(length, 10), data=oneSample)

# prediction
# order to make sure plotting looks nice
orderedLength <- sort(oneSample$length)
linearPred <- predict(linear, newdata=data.frame(length=orderedLength)) 
secondPred <- predict(second, newdata=data.frame(length=orderedLength)) 
tenthPred <- predict(tenth, newdata=data.frame(length=orderedLength)) 

lines(orderedLength, linearPred, col="blue", pch=16)
lines(orderedLength, secondPred, col="red", pch=16)
lines(orderedLength, tenthPred, col="green", pch=16)

# check fit! this is expressed by the likelihood value of the model (you can get log-likelihoods)
logLik(linear)
logLik(second)
logLik(tenth)

# The higher the polynomial, the better the fit. But it comes at the cost of more parameters (watch for the degrees of freedom).
# Model selection is based on this compromise. The Akaike Information Criterion is calculated from these two things: model fit and parsimony.

# The lower the AIC, the better the model:
AIC(linear)
AIC(second)
AIC(tenth)

# Now try to repeat this with different seeds! 

########################################----------------------------------------

# 9. Download the all.rds file, and generate a 50 element sample from it. Which is a better predictor of specimen age, length or width? 
all <- readRDS("data/all.rds")

# random sample
samp <- Collect(all, 50)

# making a mutiple regression
multiple <- lm (age ~ width + length, data=samp)

# are both important? step tries to omit parameters and searches for minimum AIC of parameter combinations
selected <- step(multiple)
summary(selected)





# + 1. Check out the function definition in unknown.R 
# - without running it, what does this function do? 
# - Run it as: TheThing(4, 5, 1, da=0.1), do not plot it, but calculate the correlation coeffficient. 
# - are the two variables related? 
