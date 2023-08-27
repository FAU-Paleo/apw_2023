---
title: ""
layout: "post" 
permalink: "phylogenetics/revbayes"
---

## Introduction to RevBayes and the Rev language

[RevBayes](https://revbayes.github.io) is a flexible program that has been developed for Bayesian phylogenetics.
It has it's own command line language, similar to the R programming language, that can be used to set up a model in a Bayesian framework. We won't run any analyses in this exercise, we'll just explore the Rev languages and look at the syntax used to set up different types of variables. 

First make sure you can [download](https://revbayes.github.io/download) and install the software.
Once you have it up and running, see if you can produce any sensible output using the RevBayes terminal. For example, it can be used like a calculator in the same way as R.

```
4 + 10
```

## Constant variables {#constant}

Constant variables have a fixed value and are specified using the left arrow `<-` assignment.

```
a <- 4
b <- 10
a + b
```

Note that in R the `=` and `<-` assignment can almost always be used interchangeably. This is not the case in RevBayes, where `=` can be used to define variables in your environment but these can't then be used as part of a model. We'll see some examples of when we might use `=` in RevBayes in future weeks.

To clear the value of a variable you can use `clear(a)`.

## Deterministic variables {#deterministic}

Deterministic variables are linked to the value of other variables and specified using the colon-equal `:=` assignment.

```
# specify constant variable a
a <- 4
# specify deterministic variable b
b := a
b
```

```
# now see what happens if you change the value of a
a <- 5
b
```

Deterministic variables can also be transformed using any given formula you provide.

```
b := 3 * a
```

> Explore what happens when you change the value of `a` and experiment with different formulas / operators.

## Stochastic variables {#stochastic}

Stochastic variables are described with respect to some underlying probability distribution and are specified using the tilde `~` assignment. For example, we might want to say that variable `c` is a normally distributed random variable.

```
# specify the mean and standard deviation parameters as constant variables
a <- 1
b <- 0.1

# specify stochastic variable using the normal distribution 
c ~ dnNormal(a, b)
c
```

Any distribution function in Rev always begins with `dn`. You can also bring up help pages using a question mark like we do in R, e.g. `?dnNormal` will bring up help pages for the normal distribution.

You can define stochastic variables in RevBayes using any common statistical distributions (as well as more complex phylogenetic distributions that we'll meet in future weeks). A list of available Rev functions can be found [here](https://revbayes.github.io/documentation/).

>Can you think of any other common distributions and try to set them up in RevBayes?

If you run the above code again, `c` will take on a different random value from the normal distribution. Without more context, it might not seem obvious why we'd want to define a variable in this way. You can think of `c` as a potential parameter with an unknown value that we're interested in estimating. When we go on to set up the full analysis in the next class, you'll see how we can update the value of `c` during our analysis. And because stochastic variables like `c` are linked to an underlying distribution, `c` inherits some of the properties of that distribution. We can call functions related to that distribution by using `.` followed by a function name. For example, `c.probability()` will tell us the probability of observing the current value of `c` given the underlying distribution.

> Redefine your stochastic variable and compare the probability for different values under the same distribution -- do these probabilities make sense to you?

## Functions

A substitution model is defined by an instantaneous rate matrix, where each character transition is associated with a rate of change (e.g., the rate of changing from A to G). We could set this up manually, by specifying an entry for each row/column in a matrix but this is a kerfuffle. Luckily, RevBayes has many inbuilt functions, which begin `fn`, including one for all the most commonly used substitution models.

The simplest substitution model for nucleotide data is the Jukes-Cantor model, which can be specified using the `fnJC` function.

```
# specify a constant matrix variable 
Q <- fnJC(4) # 4 here refers to the number of states (i.e. ATGC)
Q
```

In this model all the substitution rates are equal. If substitution rates vary, we need to estimate this variability from the data and `Q` would become a deterministic variable. We'll see some examples of this next time!

Last thing for today - let's compare the Q matrix and the transition probabilities over different intervals of time. We can start by creating a simpler 2-state Q matrix.

```
# specify a constant matrix variable 
Q <- fnJC(2)
Q
# calculate the transition probabilities after 0.1 units of time
Q.getTransitionProbabilities(0.1)
# calculate the transition probabilities after 0.1 units of time
Q.getTransitionProbabilities(100) # 0.1 units
```

### Further reading

Here's a more [comprehensive intro](https://revbayes.github.io/tutorials/intro/rev) to the Rev language.
