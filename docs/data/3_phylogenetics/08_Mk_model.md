---
title: ""
layout: "post" 
permalink: "phylogenetics/morpho"
---

## Tree Building using Mk substitution model

In this exercise we'll estimate trees using different extensions of the Mk substitution model. Your choice of substitution model is any combination of the standard Mk model, with extensions for partitioning the data, allowing for unequal transition rates (+gamma), and correcting for collection biases (MkV).
There is then a total of 8 models to choose from:

```
Mk
Mk+G
MkV
Mk+multistate
MkV+G
MkV+multistate
Mk+G+multistate
MkV+G+multistate

```
 
For this exercise you will need a number of R packages.

Note: Start with the RevBayes analysis and download these R packages when the MCMC is running.

```
library(viridisLite)
library(viridis)
library(phangorn
library(permute)
library(lattice)
library(vegan)
library(ape)
library(ggplot2)
``` 

The morphological matrix used for this tutorial can be downloaded from [here]({{site.baseurl}}/data/7_phylogenetics/Agnolin_2007a_paleobiodb.nex). This matrix is from Agnolin et al 2007 Revista del Museo Argentino de Ciencias Naturales. It consists of 51 characters and 12 taxa of Brontornis.

We will start with the analysis in RevBayes and then use R to look at the output. 

### How to organise your code 

For this exercise create a new folder with two subdirectories, `data` and `scripts`.
We will create 3 scripts:

1. `main_mcmc.rev` for reading in the data and setting up the tree model and MCMC settings
2. `Mk_model_1.rev` this script is for the first version of the Mk model you choose to run. 
3. `Mk_model_2.rev` this script is for the second version of the Mk model you choose to run. 
4. `Visual_model_comparison.r` this is to see the effects of the different models on parameter estimates

### Tree inference
The set-up of this model is similar to the JC model set up in exercise 4.

In the `main_mcmc.rev` scripts start by reading in your data 

`morpho <- readDiscreteCharacterData("data/Agnolin_2007a_paleobiodb.nex")`

Next, we need to define some helper variables that will come in handy later for setting up our model.

```
taxa <- morpho.names()
num_taxa <- taxa.size()
num_branches <- 2 * num_taxa - 3
```

We also need to define variables for the moves and monitors.

```
moves = VectorMoves()
monitors = VectorMonitors()
```


#### Uniform Tree Prior
We will use the same uniform tree prior for both Mk runs. It is important that when you are comparing two models that only one part is changed - in this case we are only changing the substitution model. Everything else should be left constant. 

```
br_len_lambda <- 10
phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(br_len_lambda)) 
tree_length := phylogeny.treeLength()
```

Add moves to explore tree space and branch lengths 

```
moves.append( mvNNI(phylogeny, weight=num_branches/2.0) )
moves.append( mvSPR(phylogeny, weight=num_branches/10.0) )
moves.append( mvBranchLengthScale(phylogeny, weight=num_branches) )

```

#### Substitution model

Now open one of your `Mk_model_*.rev` scripts and construct your substitution model.
The Mk model used the same function `fnJC()` as the JC model. We need to specify the correct number of character states here. 

`Q <- fnJC(4)`

Next we’ll define a stochastic node representing a sequence alignment and “clamp” that variable to our morphological data.

```
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="Standard")
seq.clamp(morpho)
```

This is the most simple version of the Mk model. In order to relax some of the assumptions of this model we can go on to add extensions that allow the model to more accurately represent morphological evolution.

##### Across site transition variation (+Gamma)
Allowing sites to transition from one state to another at different rates may be more realistic than the assumption that they would all have equal rates. In order to incorporate this into our model we need to first set up the Gamma-distributed rate variation

```
alpha_morpho ~ dnUniform( 0.0, 1E6 )
rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
```
Add moves on the parameters to the Gamma distribution.

`moves.append( mvScale(alpha_morpho,lambda=1, weight=2.0) )`

Then, add this to our model along with our rate matrix as before

```
seq ~ dnPhyloCTMC(tree=phylogeny, siteRates=rates_morpho, Q=Q, type="Standard")
seq.clamp(morpho)
```

##### MkV model 
This model corrects for the collection biases that are inherent in morphological data. To add this parameter, we need to tell the model that all of the characters in our matrix are variable. 

```
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="Standard", coding="variable")
seq.clamp(morpho)
```

##### Partitioning the data based on the number of character states per trait 
By partitioning the data we are allowing different characters to evolve according to different rate matrices. Using a loop, we can create N number of Q matrices for our data. 

```
n_max_states <- (4)
idx = 1
morpho_bystate[1] <- morpho
for (i in 2:n_max_states) {
    morpho_bystate[i] <- morpho                                # make local tmp copy of data
    morpho_bystate[i].setNumStatesPartition(i)                 # only keep character blocks with state space equal to size i
    nc = morpho_bystate[i].nchar()                             # get number of characters per character size with i-sized states

    if (nc > 0) {                                              # for non-empty character blocks
        q[idx] <- fnJC(i)                                      # make i-by-i rate matrix
        m_morph[idx] ~ dnPhyloCTMC( tree=phylogeny,
                                    Q=q[idx],
                                    nSites=nc,
                                   type="Standard")           # create model of evolution for the character block

        m_morph[idx].clamp(morpho_bystate[i])                  # attach the data

        idx = idx + 1                                          # increment counter
        idx
    }
}

```

This code creates an evolutionary model and clamps our model to the data. 


Use `source()` to read in your substitution models from the `main_mcmc.rev script`
### MCMC Settings 

##### Create the model variable

 

`mymodel = model(phylogeny)`
First specify which version of the Mk model you chose, for example 

`model_name="MkV"`

This is just the name of the exercise
`analysis_name="morpho_sub"`

Next, we setup our monitors, like in our previous MCMC analyses.

```
monitors.append( mnModel(filename="output/" + model_name  + "/output/" + analysis_name + "_posterior.log",printgen=10, separator = TAB) )
monitors.append( mnFile(filename="output/" + model_name +  "/output/" + analysis_name + "_posterior.trees",printgen=10, separator = TAB, phylogeny) )
monitors.append( mnScreen(printgen=1000, tree_length) )
monitors.append( mnStochasticVariable(filename="output/"  + model_name +  "/output/" + analysis_name + "_posterior.var",printgen=10) )
```


Finally, we’ll set up the MCMC run using the mcmc function, specifying our model, the vector of monitors and the vector of moves. And then we’ll run the chain for 2000 generations.

```
mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")
mymcmc.burnin(generations=200,tuningInterval=200)
mymcmc.run(generations= 2000,tuningInterval=200)

```

### Evaluating the output

In R try to recreate plots showing the difference in tree length between your two models.

The code to make the distance plots can be downloaded [here]({{site.baseurl}}/data/7_phylogenetics/Visual_model_comparison.r). 






 


