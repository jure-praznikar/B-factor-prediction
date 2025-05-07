# B-factor-prediction
This R scripts:
1) Read the PDB file "protein.pdb"
2) Convert the protein to the graph (cutoff distnce is 7Å)
3) Count the orbits for each atom to create a feature vector
   for each atom.
   The length of the feature vector is 15 (max. graphlet size = 4)
5) Use the linear GDV model and GDV matrix to predict the B-factors

First: Install and upload the following libraries:
library(netdist)
library(bio3d)
library(igraph)
library(pdist)
library(cry)

To execute (run R script from command line):
> source("main.r")

or run it (main.r) in Rstudio.

Further details can be found in the publication:
