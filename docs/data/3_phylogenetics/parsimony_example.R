# read in a character matrix in ape 
dat = ape::read.nexus.data(file = "data/Frankie-data.nex")

# identify character states
chars = unique(unlist(unique(dat)))

# & convert to phanghorn object
dat = phangorn::phyDat(dat, type = "USER", levels = chars)

# create a neighbour joining starting tree
dm = phangorn::dist.hamming(dat)
tree = phangorn::NJ(dm)

# estimate the parsimony tree
mp = phangorn::optim.parsimony(tree, dat)

# identify the outgroup - we'll use the fungus Laccaria
outgroup = "Laccaria"

# root the tree
mp = ape::root(mp, outgroup = outgroup, resolve.root = TRUE)

plot(mp)
