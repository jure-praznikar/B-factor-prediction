# **B-factor-prediction**

A linear model based on graphlet degree vectors proves to be effective <br>
not only for the prediction of B-factors and the validation of deposited <br> 
protein structures but also for the qualitative estimation of <br>
root-mean-square fluctuations derived from molecular dynamics.


This R scripts:
1) Read the PDB file "protein.pdb" <br>
2) Convert the protein to the graph (cutoff distnce is 7Å) <br>
3) Creates Graphlet Degree Vector matrix<br>
5) Use the linear GDV model and GDV matrix to predict the B-factors <br>

**First: Install following R libraries:**
* library(netdist)
* library(bio3d)
* library(igraph)
* library(pdist)
* library(cry)

**To execute (run script from R command line):** 
> source("main.r")

or run it (main.r) in RStudio.

Expected output is: <br>
*Correlation (predicted vs pdb):  0.811* <br>

**Further details can be found in the publication:** <br>
*Article* <br>
**Title:** Prediction of temperature factors in proteins: effect of data pre-processing and experimental conditions <br>
Crystals 2025, 15(5), 455; https://doi.org/10.3390/cryst15050455 <br>

