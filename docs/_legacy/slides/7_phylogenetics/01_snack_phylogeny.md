---
title: ""
layout: "post" 
permalink: "phylogenetics/snack_phylogeny"
---

## Snack phylogeny

In this exercise, we'll explore the process of building a phylogeny.

You have 7 snacks and your task is to build a **rooted tree**, assumming the snacks have a shared phylogenetic history, just like species.

How you do this is completely up to you and your group. It might help to start by choosing an outgroup and then considering what characters support each node in your phylogeny.

**Write down your logic** as you go along!

### Create a Newick version of tree

The most widely used format for storing trees is [Newick format](https://en.wikipedia.org/wiki/Newick_format). 

Brackets are used to indicate which taxa group together, representing internal nodes. Descendants are separated by commas.

```
(A, B); # rooted
((A, B), (C, D)); # rooted
(((A, B), (C, D)), E); # rooted
(A, B, (C, D)); # unrooted
```

Note the string must end in a semicolon.
In a rooted binary tree all nodes will have two descendants. 

**Store your tree in Newick format**. You can use Tim Vaughan's online tool [IcyTree](https://icytree.org) to plot the tree and see that everything is correct. You can use 'File > Enter tree directly' to copy & paste your tree.

When you're done email us your tree (rachel.warnock@fau.de).

### Acknowledgements

This exercise is based on the one described [here](https://labroides.wordpress.com/candy-phylogeny/) by Josh Drew. 
