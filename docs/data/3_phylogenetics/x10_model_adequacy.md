---
title: ""
layout: "post" 
permalink: "phylogenetics/model_adequacy"
---

## Model adequacy using posterior predictive simulations in RevBayes

In this exercise we will use posterioir predictive simulations to assess the fit of an evolutionary model to our data. Again, use the two models you chose for the past two exercises.
The code highlighted in this tutorial is just to show the more important parts of a script. You do not need to copy and paste them into the terminal. Instead, edit the scripts in your editor and use `source(" ")`to read them into RevBayes.

There are [4 rev scripts]({{site.baseurl}}/data/7_phylogenetics/model_adequacy.zip) for the analysis:
 
1. `MCMC_Simulation.rev` runs an MCMC on your empirical dataset and samples parameters from the posterior
2. `PP_Simulation.rev` takes the output from the first MCMC to simulated new datasets under the chosen model
3. `PP_MCMC.rev` runs an MCMC for all the new datasets
4. `Tree_Summary.rev` generates test statistics


and 1 Rscript `PPS_Analsis.r` to summarise the output 

The data used here is a morphological matrix Brontornis taxa, an extinct genus of giant flightless birds from Argentina. 

### Read in the data and define the subsitution model to test 
Start with `MCMC_Simulation.rev` script. We first ned to read in our data and define which substitution model we will use for the analysis, either Mk or Mkv+G.

```
morpho <- readDiscreteCharacterData("data/Agnolin_2007a_paleobiodb.nex")
analysis_name = "model_adequacy" 
model_name = model # what model are you testing? e.g. "mk"
model_file_name = "scripts/"+model_name+"_Model.Rev"
```
**Note**: These varaibles are at the start of the rev scripts 1-4. Make sure you add the model you want to use for each script.

This monitor samples the posterioir and saves it to a file that can be used to simulate new datasets under the chosen model 

```
monitors.append( mnStochasticVariable(filename="output/"  + model_name +  "/output/" + analysis_name + "_posterior.var",printgen=10) )
```

### Mk substitution model
Helper variables and branch lengths are set up as in previous exercises. 

Specify the Jukes-Cantor substitution model applied uniformly to all sites. Remeber the mk model is a generalization of the jukes-Cantor model. 

```
Q := fnJC(4) # this matrix has 4 different character states. 

seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type="Standard")
```
 
### MkV + Gamma substitution model
Set up gamma-distributed rate variation 

```
alpha_morpho ~ dnUniform( 0.0, 1E6 )
rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )
moves.append( mvScale(alpha_morpho,lambda=1, weight=2.0) )
``` 
   
To account for asertainment bias. Change the coding to `variable`.

```
seq ~ dnPhyloCTMC(tree=phylogeny, siteRates=rates_morpho, Q=Q, type="Standard",  coding="variable")
``` 

    
Once you have chosen your model you can source the first script.

```
source("scripts/1_MCMC_Simulation.Rev)
``` 
  
### Simulating new datasets
We can now simulate new datasets. In this script pay attention to the `thinning` variable. 
This number controls the number of datasets simulated. We set it too 100 here so every 100th tree is used in the .var file. In a standard analysis you would set this a lot lower (e.g. 2). 

```
pps.run(thinning=100)
```

```
source("scripts/2_PP_Simulation")
```

### MCMC for simulated datasets
Next we need to run an MCMC for all of our simulated datasets. This function runs an MCMC for all datasets in the directory file path. 

```
my_pps_mcmc = posteriorPredictiveAnalysis(mymcmc, directory)
my_pps_mcmc.run(generations=1000)
```

```
source("scripts/3_PP_MCMC.Rev")
```

### Generating test statistics
We can then generate some test statisitcs to asses how well models fit our data. This script caluclates the tree length and the Robinson foulds distance for both the empirical dataset and simulated datasets. 

```
source("scripts/4_Tree_Summary.Rev")
```

Now run the whole analysis again using the other morphoylogical model

### Producing figures
Use the rscript `PPS_Analysis.r`
If everything has been run correctly up to now you only need to set your directory and run this file.

### Analysing the output 

The rscript should produce a figures directory containing 3 files:

1. `Boxplots.pdf` This figure contains boxplots of the effect sizes for the 4 test statistics used here. The closer the plot is to zero the more similar the simulated results are to the empirical 

 
2. `Histograms.pdf` This figure shows the distribution of tree length and Robinson Foulds distance across the simulated datasets. The blue line shows the empirical values and the pink  broken line shows the mean of the simulated. 


3. `results.csv` Stores the effect sizes for each model.
