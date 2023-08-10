---
title: ""
layout: "post" 
permalink: "phylogenetics/fbd_sim"
---

## Simulating trees and fossils 

Simulations are a valuable tool for exploring and understanding models.
In this exercise we'll simulate trees and fossils under birth-death sampling processes. We'll focus on constant rate models (i.e. rates won't change over time) and explore the complete vs. reconstructed tree.

We'll use the R packages `TreeSim` and `FossilSim`. `TreeSim` can be used to simulate trees under a range of models trees, while `FossilSim` can be used to simulate and plot fossil data under a range of preservation models.

First install the packages.

```
install.packages("TreeSim")
install.packages("FossilSim")
```

Some questions to think about as you go through the examples:

> Can you see the differences between trees based on the model assumptions?
> How do different parameters (e.g. birth, death) impact the trees?


## The pure birth process

Let's start with a pure birth process. In this model we have some speciation rate and assume extinction is zero. 

In the following we specify the `birth` (> 0) and `death` (= 0) parameters (which you can think of as the speciation and extinction rates). We can control the size of tree by condition the process on the number of extant tips.

```
set.seed(123)

birth = 0.1
death = 0
tips = 10

num_sim = 1

trees = TreeSim::sim.bd.taxa(tips, num_sim, birth, death)

t1 = trees[[1]]
```

Note that although we don't state this explicity, the above assumes that sampling at the present (i.e. extant taxon sampling) is complete.

Next let's plot the output. We'll construct an empty fossils vector, so we can use the `FossilSim` plotting functions. 

```
f1 = FossilSim::sim.extant.samples(FossilSim::fossils(), t1)

plot(f1, t1) # complete tree
```

You can add the argument `reconstructed=TRUE` to the plotting function to show the reconstructed tree, although the complete and reconstructed trees are identical in the case of a prue birth process. 

## The birth-death process

Next let's incorporate extinction, `
death` (> 0). And we'll continue to assume that extant species sampling is complete.

```
birth = 0.1
death = 0.05
tips = 10

num_sim = 1

trees = TreeSim::sim.bd.taxa(tips, num_sim, birth, death)

t1 = trees[[1]]
```

This time you should see a difference between the complete and reconstructed trees.

```
plot(f1, t1) # complete tree 

plot(f1, t1, reconstructed = TRUE) # reconstructed tree
```

## The birth-death sampling process

We can relax the assumption that extant speciation was complete using the `FossilSim` function `sim.extant.samples` and introducing the sampling probability parameter rho, for extant species. 

We can start with the birth-death tree we simulated above.

```
rho = 0.5

f2 = FossilSim::sim.extant.samples(f1, t1, rho = rho)

plot(f2, t1)

plot(f2, t1, reconstructed = TRUE)

```

## The fossilised birth-death process

The FBD model also incorporates the fossils sampling process, through the addtional of the fossil sampling rate parameter. 

```
psi = 0.1

f3 = FossilSim::sim.fossils.poisson(psi, fossils = f2, tree = t1)

plot(f3, t1)

plot(f3, t1, reconstructed = TRUE) 
```

> Can you increase / decrease the proportion of sampled ancestors?