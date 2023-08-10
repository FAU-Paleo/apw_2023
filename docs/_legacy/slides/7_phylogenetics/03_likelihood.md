---
title: ""
layout: "post" 
permalink: "phylogenetics/likelihood"
---

## Tree building using maximum likelihood

In this exercise we'll estimate trees using maximum likelihood and look at the impact of different substitution models.

For this we'll use the `phangorn` package and a data set of primates sequences.
The sequences are a short fragment of the mitochondrial genome. 
You can download the data from [here]({{site.baseurl}}/data/7_phylogenetics/primates-data.nex).

First install the necessary R packages if you don't have them aleady. 

```
install.packages("ape")
install.packages("phangorn")
```

Start by reading in sequence the file and converting the data into the right format for `phangorn`.

```
# read in the character matrix in ape 
dat = ape::read.nexus.data(file = "data/primates-data.nex")

# & convert to phanghorn object
dat = phangorn::phyDat(dat, type = "DNA")
```

As [before]({{site.baseurl}}/phylogenetics/likelihood) we need to generate a starting tree.
We can do this using an approach called **neighbour-joining**. This simply calculates the evolutionary distance between the taxa in your matrix (using the `dist.hamming` function) and then groups taxa together based on their relative distances (using the `NJ` function).

*In case of interest, for two sequences the [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) is simply the proportion of sites at which they differ. It's therefore one of the simplest ways to approximate evolutionary distances and generate a starting tree. But note this approach won't account for convergent evolution or homoplasy.*

```
# create a neighbour joining starting tree
dm = phangorn::dist.hamming(dat)
starting_tree = phangorn::NJ(dm)
```

The `phangorn` function `pml` (pml = phylogenetic maximum likelihood) calculates the likelihood of a tree given the data under a model specified by the user. We'll start with the simplest model of sequence evolution, the Jukes-Cantor (JC) model.

```
# calculate the likelihood of the tree given the data
fitJC = phangorn::pml(starting_tree, data=dat, model = "JC")
fitJC
```

> What do you notice about the base frequencies?

Then we'll use the function `optim.pml` to find the tree with the best likelihood score under the Jukes-Cantor model. This function takes the tree and the starting values obtaining using the `pml` function and explores the tree space by switching around the branches. The argument `optNni = TRUE` tells the function to optimize the tree topology (i.e. find the tree with the highest likelihood), as well as the branch lengths and the model parameters. The parameters of the JC model are fixed, so the model parameters won't actually change in this case.

```
# estimate the tree using ML
fitJC_opt = phangorn::optim.pml(fitJC, model = "JC", optNni = TRUE)
fitJC_opt
```

As the ML algorithm proceeds it outputs information about the progress. You can see how the likelihood changes as the tree topology and the branch lengths (edge weights) are updated. You should also see that the likelihood is improving, i.e. getting larger.

The function outputs a list with a bunch of values that might be of interest, including the tree.

```
fitJC_opt$tree
```

Next we'll do the same with the GTR substitution model. 

```
fitGTR = phangorn::pml(starting_tree, data=dat, model = "GTR")
fitGTR_opt = phangorn::optim.pml(fitGTR, model = "GTR", optNni = TRUE)
fitGTR_opt
```

> What differences do you notice in the output?

Let's look at the trees. Remember that by default trees generated based on a substitution model will be unrooted. First, let's root the tree using the Lemur.

```
# identify the outgroup
outgroup = "Lemur"

# root the tree
JC_tree_rooted = ape::root(fitJC_opt$tree, outgroup = outgroup, resolve.root = TRUE)
GTR_tree_rooted = ape::root(fitGTR_opt$tree, outgroup = outgroup, resolve.root = TRUE)

# plot the trees
plot(JC_tree_rooted)
plot(GTR_tree_rooted)
```

> Do you notice any differences in the topology or branch lengths? Which one is correct?

Since we have the maximum likelihood estimate under both models, we can use a straightforward model testing to identify the best fitting model. 

```
AIC(fitJC_opt, fitGTR_opt)
```

> How would you interpret these results? [Hint](https://www.scribbr.com/statistics/akaike-information-criterion/)

## Output a tree from R

To output a tree from R and save it we can use the `ape` function `write.tree`. This outputs the tree in Newick format.

```
# write to screen
ape::write.tree(tree)
```

```
# write to file
ape::write.tree(tree, file = "my-tree.nex")
```

## Other resources for buidling trees using maximum likelihood

`phangorn` is a very useful package but you won't often see published trees generated using this package. There are other more specialized software packages for tree building using maximum likelihood. Some of these feature additional substitution models or have been optimized for the analysis of large data sets. The analysis in this exercise is very fast because our data set has a small number of species and characters (only 150) -- most DNA datasets are *much* larger than this. Two popular alternative likelihood based programs for tree building are  [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) and [IQTREE](http://www.iqtree.org/).

A complete script for this exercise can be downloaded [here]({{site.baseurl}}/data/7_phylogenetics/ML_example.R).