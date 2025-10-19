rm(list=ls())

### LIBRARY
library(netdist) # GDV, count orbits
library(bio3d) # manipulate with pdb files 
library(igraph) # create graphs, ...
library(pdist) # pairwise distance calcualtion
library(cry) # crystallographic symmetry
library(MASS) # ginv: Generalized Inverse of a Matrix

### FUNCTIONS (required)
source("FUNC_BIOextract.r") # extract/create bio chains/data
source("FUNC_symmetry.r") # symmetry matrix
source("FUNC_add_symm_atoms.r") # add symmetry atoms (include ligands)
source("FUNC_GDV_perpartes_and_all_at_once.r") # create GDV=Graphlet Degree Vector
source("FUNC_scaleBf.r") # scale B-factors

### INPUT file - PDB format
pdbfile<-'protein.pdb' # HYDROLASE(O-GLYCOSYL) PDBid: 153L

### run GDV (create graph, and count orbits per atom >> GDV)
source("GDV.r") # output file: protein.dataGDV >> GDV Graphlet Degree Vector

### Predict; use linear GDV model
source("predict.r")
# ------------------------------
# Correlation (predicted vs pdb):  0.809 
# ------------------------------

