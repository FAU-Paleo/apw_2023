---
title: ""
layout: "post" 
permalink: "phylogenetics/PPS"
---

# Model adequacy using posterior predictive simulations in RevBayes

In this exercise we will use posterioir predictive simulations to assess the fit of an evolutionary model to our data. Again, you can test two models.
The code highlighted in this tutorial is just to show the more important parts of a script. You do not need to copy and paste them into the terminal. Instead, edit the scripts in your editor and use `source(" ")`to read them into RevBayes.

The PPS workflow is as follows
 
1. `pps_MCMC.rev` runs an MCMC on your empirical dataset and samples parameters from the posterior
2.  `MorphoSim.r` simulates data in r using information sampled in 1. 
3. `pps_MCMC_Simulations.Rev` runs MCMC on simulated data sets
4. `Tree_Summary.rev` generates test statistics
5. `SummaryStats.r` Calculates effect sizes

Use the file `run.sh` for the commands to run this from your terminal. 

__Make sure the filename, model, and number of states is correct for each of the rev scripts used. You need to manually change this in each of the three files__


### Read in the data and define the subsitution model to test 
Start with `pps_MCMC.rev` script. We first ned to read in our data and define which substitution model we will use for the analysis.

```
morpho <- readDiscreteCharacterData("data/Egi_etal_2005a_paleobiodb.nex)
model_name = model # what model are you testing? e.g. "mk"
model_file_name = "scripts/"+model_name+"_Model.Rev"
num_states="5"
```
**Note**: These varaibles are at the start of all rev scripts. Make sure you add the model you want to use for each script.

This monitor samples the posterioir and saves it to a file that can be used to simulate new datasets under the chosen model 

```
monitors.append( mnStochasticVariable("output/" + model_name + "/output/APW-example.var",printgen=10) )
```

### Mk substitution model
Helper variables and branch lengths are set up as in previous exercises. 

Specify the Jukes-Cantor substitution model applied uniformly to all sites. Remeber the mk model is a generalization of the jukes-Cantor model. 

```
Q := fnJC(init(num_states) 
seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="Standard")
```
 
### MkV + Gamma substitution model
Set up gamma-distributed rate variation 

```
alpha_inv ~ dnExponential(1)
alpha_morpho := 1 / alpha_inv
rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 ) 

# Moves on the parameters to the Gamma distribution.
moves.append( mvScale(alpha_inv,lambda=1, weight=2.0) )
``` 
   
To account for asertainment bias. Change the coding to `variable`.

```   
seq ~ dnPhyloCTMC(tree=phylogeny, siteRates=rates_morpho, Q=Q, type="Standard",  coding="variable")
``` 

    
Once you have chosen your model you can source the first script.

``` 
source("scripts/pps_MCMC.Rev)
``` 
  
### Simulating new datasets
We can now simulate new datasets using the `MorphoSim.r`script. You have to specify how many simulated data sets you want. The standard is 1000 but for the class you can use 100 so it will finish faster. 

``` 
Rscript MorphoSim.r Egi_etal_2005a_paleobiodb.nex mk+GV 5
```

### MCMC for simulated datasets
Next we need to run an MCMC for all of our simulated datasets. This function runs an MCMC for all datasets in the directory file path. 

``` 
my_pps_mcmc = posteriorPredictiveAnalysis(mymcmc, directory)
my_pps_mcmc.run(generations=1000)
```

``` 
source("scripts/pps_MCMC_Simulations.Rev")
```

### Generating test statistics
We can then generate some test statisitcs to asses how well models fit our data. This script caluclates the tree length and the Robinson foulds distance for both the empirical dataset and simulated datasets. 

``` 
source("scripts/pps_TreeSummary.Rev")
```

Now run the whole analysis again using the other morphoylogical model

### Producing figures
Use the rscript `SummaryStats.r ` to generate plots of the effect sizes. __Remember, effect sizes closer to zero indicate a good fit__.
This script will take which ever models you analysed (and are now in the output directory) and caluclate effect sizes. 


### Analysing the output 

The rscript should produce a figures directory containing 2 files:

1. `Test_Statistics.pdf` This figure contains boxplots of the effect sizes for the 6 test statistics used here. The closer the plot is to zero the more similar the simulated results are to the empirical 

2. `Test_Statisitcs.csv` Stores the effect sizes for each model.




