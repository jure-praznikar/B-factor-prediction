# **B-factor-prediction**
This R scripts:
1) Read the PDB file "protein.pdb" <br>
2) Convert the protein to the graph (cutoff distnce is 7Å) <br>
3) Count the orbits for each atom to create a feature vector <br> 
   for each atom. <br>
   The length of the feature vector is 15 (max. graphlet size = 4). <br>
5) Use the linear GDV model and GDV matrix to predict the B-factors <br>

**First: Install following R libraries:**
* library(netdist)
* library(bio3d)
* library(igraph)
* library(pdist)
* library(cry)

**To execute (run R script from command line):** 
> source("main.r")

or run it (main.r) in RStudio.

Expected output is: <br>
*Correlation (predicted vs pdb):  0.811* <br>

**Further details can be found in the publication:** <br>
*Article* <br>
**Title:** Prediction of temperature factors in proteins: effect of data pre-processing and experimental conditions <br>
*Jure Pražnikar*
