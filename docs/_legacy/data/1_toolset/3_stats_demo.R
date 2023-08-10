# 0. Download and load the file."bed1_length_sample.rds"

# 1. Analyze and plot the distribution from the file.
# - Plot the distribution.
# - Look at summary statistics
# - Do a hypothesis test.

# 2. Load "bed1_length.RData". This make-belief example represents our population.
# - Generate a single random sample of 30 specimens!
# - Write a function to execute this sampling with a given sample size!


# 3. With our new function, create 3 differently sized samples containing 20, 50 and 100 elements! Bind them together in a list!
# - name the list, so the names contain how many elements the list has.
# - calculate the mean and standard deviation of the samples and put the results in a 3 by 2 matrix. 


# 4. Repeat 3 for every sample size between 5 and 400 with a for loop.
# - find the highest and lowest realized value across the samples!
# - calculate the means in the samples, store them in a named vector (names are sample sizes). Plot the means as a function of sample size! 
# - calculate the minimum and maximum (range of) values it every sample.
# Draw the maximum and minimum as a function of sample size! 

# 5. Add the width of the specimens, "bed1_width.RData".
# - Load the file!
# - Create a data.frame and match the width with the length measurements!
# - adjust the Collect function so it works, when you use data.frames as an input!

# 6. Create a 20 element sample and:
# - draw a scatterplot
# - calculate Pearson's covariance and correlation coefficient

# 8. Following from 6. Create 1000 different samples with 20 elements. Calculate the correlation coefficient from them, check out the distribution!

# 8. Make a linear model for one of these samples. Also fit a 2nd and a 1th order polynomial. Which is the best model and why?    

# 9. Download the all.rds file, and generate a 50 element sample from it. Which is a better predictor of specimen age, length or width? 

# + 1. Check out the function definition in unknown.R 
# - without running it, what does this function do? 
# - Run it as: TheThing(4, 5, 1, da=0.1), do not plot it, but calculate the correlation coeffficient. 
# - are the two variables related? 
