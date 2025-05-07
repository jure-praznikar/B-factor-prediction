### Read data
data<-read.table('protein.dataGDV', header = FALSE, sep = ",")
# GDV are in column from 5 to 19
dataGDV<-scale(data[,5:19]) # GDV
# coulmn 4 > scaled B-factors
BFs<-data[,4] # scaled Bfactors
# column 3 > original "PDB" B-factors
# column 2 > number of residues
# column 1 > the name of pdb file

### Predict (linear GDV model) ##################
### Regression coefficients (crystall contacts, log-transform, scale-per-chain, ligand)
beta0 <- 0 
beta <- c( 1.129154e+00,-2.017318e+00,2.022880e-01,-2.954882e+00,-7.716097e-01,
           7.196686e-01,4.185388e-01,-3.315839e-01,-2.859857e-02, 3.738797e-01,
           3.275814e-01,2.474782e-01, 6.646913e-01, 5.175808e-01, 7.514679e-01 )

### BFPR predicted B-factors
BFPR <- beta0 + dataGDV %*% beta # predicted B-values (linear model)
CC<-cor(BFPR,BFs) #BFs > scaled and log transform !!!
cat(data[1,1],' \n')
cat('------------------------------\n')
cat('Correlation (predicted vs pdb): ',sprintf("%.3f", CC),'\n')
cat('------------------------------\n')
cat('\n')
