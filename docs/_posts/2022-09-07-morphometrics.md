---
title: "Module 8: Geometric morphometrics"
layout: "post" 
permalink: "morphometrics/"
---

Welcome to Module 8: Geometric Morphometrics! Across these two days, we'll focus on why and how to collect different types of shape data and the challenges of analysing high-dimensional trait data. 


| When   | What                                                                 |
|--------|----------------------------------------------------------------------|
| Sept 7 | Intro and landmark-based GMM                                         |
| Sept 8 | Outlines and beyond                                                  |


<br>

#### Purpose
- Learn how to collect 2D landmark and 2D outline data
- Learn how to explore geometric morphometric data 
- Learn how to test hypotheses with GMM data
- Learn how to fisualize shape change

### Slideshows
- [1. Introductions]({{site.baseurl}}/data/8_morphometrics/powerpoints/1_Who_Am_I.pptx)
- [2. A Brief History of Morphometrics]({{site.baseurl}}/data/8_morphometrics/powerpoints/2_A_Brief_History_of_Morphometrics.pdf)
- [3. Outline Data Analysis]({{site.baseurl}}/data/8_morphometrics/powerpoints/3_Outline_Analysis.pptx)
- [4. Advanced Topics in GMM]({{site.baseurl}}/data/8_morphometrics/powerpoints/4_Advanced_GMM.pptx)

### R code

- [0. Setup]({{site.baseurl}}/data/8_morphometrics/exercises/0_Setup.pdf)
- [1. Collecting Outline data]({{site.baseurl}}/data/8_morphometrics/exercises/1_Outline_Data_Collection.pdf)
- [2. Conduct Ellipitcal Fourier Analysis]({{site.baseurl}}/data/8_morphometrics/exercises/2_Ellipitcal_Fourier.pdf)
- [3. Analyze EFA Data]({{site.baseurl}}/data/8_morphometrics/exercises/3_%20Analyzing_Outline_Data.pdf)
- [4. Collecting and Analyzing 2D landmark data]({{site.baseurl}}/data/8_morphometrics/exercises/4_Collecting_and_Analyzing_Landmark_Data.pdf)
- [5. Practice and Synthesis]({{site.baseurl}}/data/8_morphometrics/exercises/5_Practicing_with_3D_Landmark_Data.pdf)

### Data Files

- [Raw Data- Belemnite Outlines]({{site.baseurl}}/data/8_morphometrics/Data/Belemnite_Data.txt)
- [Smoothed Data- Belemnite Outlines]({{site.baseurl}}/data/8_morphometrics/Data/Belemnite_SmoothedOutline.nts)
- [3D Mesh file- Canis lupis]({{site.baseurl}}/data/8_morphometrics/Data/Canis_lupus.ply)
- [2D landmark scheme for mustelids]({{site.baseurl}}/data/8_morphometrics/Data/landmark_scheme.txt)
- [Links for 2D landmark scheme]({{site.baseurl}}/data/8_morphometrics/Data/links.csv)
- [3D landmarks on mammals]({{site.baseurl}}/data/8_morphometrics/Data/mammals.csv)
- [Landmark IDs for 3D data]({{site.baseurl}}/data/8_morphometrics/Data/mammal_3d_fixed_points.csv)

### R Functions

- [MorphoFiles_Function.r]({{site.baseurl}}/data/8_morphometrics/utility_functions/MorphoFiles_Function.r)
- [MorphometricExtraction_Functions.r]({{site.baseurl}}/data/8_morphometrics/Data/MorphometricExtraction_Functions.r)
- [OutlineAnalysis_Functions.r]({{site.baseurl}}/data/8_morphometrics/Data/OutlineAnalysis_Functions.r)
