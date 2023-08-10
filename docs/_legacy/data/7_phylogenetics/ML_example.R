
# read the data file  
dat = ape::read.nexus.data(file = "data/primates-data.nex")

# convert to phanghorn object
dat = phangorn::phyDat(dat, type = "DNA")

# create a neighbour joining starting tree
dm = phangorn::dist.hamming(dat)
starting_tree = phangorn::NJ(dm)

# calculate the likelihood of the tree given the data
fitJC = phangorn::pml(starting_tree, data=dat, model = "JC")
fitJC

# let's see if we can find a better tree
fitJC_opt = phangorn::optim.pml(fitJC, model = "JC", optNni = TRUE)

# let's do the same under the GTR model
fitGTR = phangorn::pml(starting_tree, data=dat, model = "GTR")
fitGTR_opt = phangorn::optim.pml(fitGTR, model = "GTR", optNni = TRUE)

# identify the outgroup
outgroup = "Lemur"

# root the tree
JC_tree_rooted = ape::root(fitJC_opt$tree, outgroup = outgroup, resolve.root = TRUE)
GTR_tree_rooted = ape::root(fitGTR_opt$tree, outgroup = outgroup, resolve.root = TRUE)

# plot the trees
plot(JC_tree_rooted)
plot(GTR_tree_rooted)

# you can also compare the ML trees to the starting tree
starting_tree_rooted = ape::root(starting_tree, outgroup = outgroup, resolve.root = TRUE)
plot(starting_tree_rooted)

# compare the models using the AIC
AIC(fitJC_opt, fitGTR_opt)

# write to screen
ape::write.tree(tree)
# write to file
ape::write.tree(tree, file = "my-tree.nex")

