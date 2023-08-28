---
title: ""
layout: "post" 
permalink: "phylogenetics/morpho_inf"
---

## Tree Inference using Mk substitution model

In this exercise we will  estimate trees using the Mk model and some common extensions. The extensions we will explore here include allowing for unequal transition rates (+gamma), correcting for collection biases (MkV) and partitioning the data.
There is then a total of 8 models to choose.

```
Mk
Mk+G
Mk+V
Mkmultistate
Mk+GV
Mk+Vmultistate
Mk+Gmultistate
Mk+GVmultistate
```

Choose two substitution models and see if they change key parameter estimates. 
 
### Tree inference

In the `MCMC.rev` script start by reading in your data. You can use the data provided or you can use your own data or a data set from Graeme T Llyod website [here](http://www.graemetlloyd.com/matr.html).

```
morpho <- readDiscreteCharacterData("data/Egi_etal_2005a_paleobiodb.nex")
```

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

#### Tree Prior

We will use the same tree prior for both Mk runs. It is important that when you are comparing two models that only one part is changed - in this case we are only changing the substitution model. Everything else should be left constant. 

```
tree_length ~ dnExponential(1)
relative_bl ~ dnDirichlet(rep(1, num_branches))
bl := tree_length * relative_bl
topology ~ dnUniformTopology(taxa)
phylogeny := treeAssembly(topology, bl)

```

Add moves to explore tree space and branch lengths 

```
moves.append( mvSlide(tree_length) )
moves.append( mvBetaSimplex(relative_bl, weight = num_branches / 2.0) )
moves.append( mvNNI(topology, weight = num_branches / 2.0) )
moves.append( mvSPR(topology, weight = num_branches / 10.0) )

```

#### Substitution model

Now open one of your `Mk_model_*.rev` scripts and construct your substitution model.
The Mk model used the same function `fnJC()` as the JC model. We need to specify the correct number of character states here. 

```
Q <- fnJC(5)
```

Next we’ll define a stochastic node representing a sequence alignment and “clamp” that variable to our morphological data.

```
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="Standard")
seq.clamp(morpho)
```

This is the most simple version of the Mk model. In order to relax some of the assumptions of this model we can go on to add extensions that allow the model to more accurately represent morphological evolution.

##### Across site transition variation (+Gamma)

Allowing sites to transition from one state to another at different rates may be more realistic than the assumption that they would all have equal rates. In order to incorporate this into our model we need to first set up the Gamma-distributed rate variation.

```
alpha_inv ~ dnExponential(1)
alpha_morpho := 1 / alpha_inv
rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 ) 
```
Add moves on the parameters to the Gamma distribution.

`moves.append( mvScale(alpha_inv,lambda=1, weight=2.0) )`

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
n_max_states <- (5)
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


Use `source()` to read in your substitution models from the `MCMC.rev` script.

### MCMC Settings 

Note that the model variable is added to the end of each model script.

`mymodel = model(phylogeny)`

In the main MCMC script, specify which version of the Mk model you chose, for example,

`model_name="MkV"`

Next, we setup our monitors, like in our previous MCMC analyses.

```
monitors.append( mnModel(filename= "output/APW-example.log",printgen=10, separator = TAB, exclude = ["bl", "relative_bl"]) )
monitors.append( mnFile(filename= "output/APW-example.trees",printgen=10, separator = TAB, phylogeny) )
monitors.append( mnScreen(printgen=1000, tree_length) )
```

Finally, we’ll set up the MCMC run using the MCMC function, specifying our model, the vector of monitors and the vector of moves. And then we’ll run the chain for 2000 generations.

```
mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")
mymcmc.burnin(generations=200,tuningInterval=200)
mymcmc.run(generations= 2000,tuningInterval=200)

```

