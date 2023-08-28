---
title: ""
layout: "post" 
permalink: "phylogenetics/model_selection"
---

## Model Selection using Stepping Stone 

In this exercise we will use stepping stone to see which of our models is a better fit for our data. This is done using a stepping stone approach. 


### How to organise your code
For this tutorial copy your `mcmc_main.rev`script from exercise 1 and rename it to `stepping_stone.rev`

We can leave most of this file the same as before but we need to remove the MCMC settings from the code.

Under where we defined the model add these lines to set up the stepping stone analysis. 

```
### this parameter cats should be changed. Maybe to 20 something. 
### 127 would be recommended for actually use but could take a long time to run
pow_p = powerPosterior(mymodel, moves, monitors, "output/" + model_name  + "/output_model_selection/" + analysis_name +"_model1.out", cats=127) 
pow_p.burnin(generations=10000,tuningInterval=1000)
pow_p.run(generations=1000)
```

Calculate the marginal likelihood from the stepping stone analysis 

```
ss = steppingStoneSampler(file= "output/" + model_name  + "/output_model_selection/" + analysis_name + "_model1.out", powerColumnName="power", likelihoodColumnName="likelihood")
ss_mk = ss.marginal()

```

Save this number as a csv file so you can compare it to the output from your other model.

```
write(ss.marginal() , filename= "output/" + model_name  + "/output_model_selection/" + analysis_name + "_marginal_likelihood.csv", "\n", append=TRUE)
```


When you have the marginal likelihoods for both models you can caluclate the bayes factor in support for mk_one as

    ss_mk_one - ss_mk_two 
    
Interpret whether this difference is meaningful using this table. 


![img]({{site.baseurl}}/images/Interpreting_Bayes_factors.png)
