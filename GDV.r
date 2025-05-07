
##################################################################################
start_time <- Sys.time()    
cat(' -----------------------------------------------\n')
PDBid<-paste(substring(pdbfile, first=1, last=4))
cat('>>> ',pdbfile,'\n')
pdb<-read.pdb(pdbfile)
PDBchain<-unique(pdb$atom$chain)
cat('PDB chains: ',PDBchain,'\n')
##################################################################################
# BIO unit section
BIO<-NA
flagBIO<-FALSE
BIO<-BIOextract(pdb) # function BIOextract
if (is.list(BIO)){
   flagBIO<-TRUE    
   CHAIN<-unique(BIO[[1]]$atom$chain) # extract BIO chains
   CHAIN<-intersect(PDBchain,CHAIN) # is in "both" BIO and PDB
   cat('BIO unit chains: ',CHAIN,'\n')
} else {
    cat('Continue whitout biounit. \n')
}
cat('flagBIO is ',flagBIO,'\n')    
## BIO end
##################################################################################
## SYMM function checkSYM
## check if it is possible to extract symmetry matrices directly from PDB header OR
## by using cell data and cry lib  
flagSYM<-checkSYM(pdbfile) # if true, consider close crtsyall contacts
#flagSYM<-FALSE  # !! do not use symm. contacts !! Force not to use symmetry contacts

# LIGAND flag
flagLIGNUC<-TRUE # TAKE into account ligand and nucleic acid
#flagLIGNUC<-FALSE # DO NOT take into account ligand and nucleic acid  

cat('flagSYM is ',flagSYM,'\n')    
cat('flagLIGNUC is ',flagLIGNUC,'\n')
        
if (flagSYM){
      cat('Add close contact atoms.\n')
      pdb<-addSMYTatoms(pdb,pdbfile,flagLIGNUC) # function addSMYTatoms
      cat('Protein atoms: ',sum(pdb$atom$chain!="9"),'\n') # symmetry have chain id="9"
      cat('Close contact atoms: ',sum(pdb$atom$chain=="9"),'\n') # symmetry have chain id="9"
    } else {    
      cat('No symmetry was found. Continoue whiteout close contacts. \n')     
}    
## SYM end  
#################################################################################
#
end_time1 <- Sys.time()
cat('Time addSMYTatoms:\n')
print(end_time1-start_time)
################################################################################
## protein atoms, and no Hydrogens
# biounit (header info) is lost after trim, that's way is here and not at the begining
indsP <- atom.select(pdb, "protein",chain="9",operator="OR") # !!! LIGAND and NUCLEIC = "9"
indsH <- atom.select(pdb, "h") # no hydrogen atoms
indsTRIM <- combine.select(indsP, indsH, operator="-")
pdb<-trim.pdb(pdb, inds = indsTRIM) 
PDBprotein<-setdiff(unique(pdb$atom$chain),"9") #exclude cryst. "contact" chain=9
cat('Protein chains: ',PDBprotein,'\n')

### select biounit atoms/chains
if (flagBIO){ # if BIOunit=True
   kb<-numeric(0)
   for (b in 1:length(CHAIN)) {
        kb<-c(kb,which(pdb$atom$chain==CHAIN[b])) #only bio chains
        kb<-sort(kb) # in case of A B C A B C >> for A chain 2211 2212 2213 9714 9715 9716
   }   
} else {
        kb<-which(pdb$atom$chain!="9") # all except cryst contact chain="9"
}
cat(length(kb),' atoms was selected. \n')
### END select biounit atoms
#
### GDV, Graphlet Degree Vector
# !! kb is needed to discriminate between the protein and symmetry protein atoms
if ( length(kb) > 10000 ){ # for large proteins (Nr. atoms>10000) do "per partes" otherwise "all at once"    
# PER PARTES ############################
  cat('------------ PER PARTES \n')
  goSM<-GDVperpartes(pdb,kb)
  cat('Done ------- PER PARTES \n')
# end PER PARTES ########################
} else {
# NEW ALL AT ONCE #######################
  cat('------------ ALL AT ONCE \n')
  goSM<-GDVallatonce(pdb,kb)
  cat('Done ------- ALL AT ONCE \n')
  # end ALL AT ONCE #######################
}

# select BIO CHAIN (and not "9"; symmetry protein atoms)
if (flagBIO){
    CHAIN<-unique(BIO[[1]]$atom$chain)
    CHAIN<-intersect(PDBprotein,CHAIN) # is in "both" BIO and PDB-protein
    } else {
    CHAIN<-setdiff(unique(pdb$atom$chain),"9") #exclude cryst. "contact" chain=9 
}
cat('Chains (final) selected for BF scaling: ',CHAIN,'\n')
# end select BIO CHAIN (and not "9")
    
### LOG and SCALE per chain flags
### LOG flag  
flagLOG<-TRUE # default
#flagLOG<-FALSE 
### ScaleCHAIN flag
flagScaleCHAIN<-TRUE # default
#flagScaleCHAIN<-FALSE
###   
cat('flagLog is ',flagLOG,'\n')    
cat('flagScaleCHAIN is ',flagScaleCHAIN,'\n')
    
# four possible options
if (flagScaleCHAIN) { # flagScaleCHAIN is TRUE; flagLOG is TRUE or FALSE
      BFs<-BFscale(pdb,CHAIN,flagLOG) # function BFscale  (chain and/or log)
    } else { # flagScaleCHAIN is FALSE; flagLOG is TRUE or FALSE
        inds <- atom.select(pdb,chain=CHAIN) # CHAIN (bio but not "9")
        if (flagLOG) { # scale+log
            BFs<-scale(log(pdb$atom[inds$atom,"b"])) 
        } else { # scale (no log)
            BFs<-scale(pdb$atom[inds$atom,"b"])
        }        
}
## end BF, log and scale per CHAIN

#### write file
# count residues (number of CA atoms) 
inds <- atom.select(pdb,chain=CHAIN,elety="CA",opreator="AND")
RESIN<-length(inds$atom) # total number of residues
# orignal/deposited Bfactors
inds <- atom.select(pdb,chain=CHAIN)
BForg<-pdb$atom[inds$atom,"b"]
df<-data.frame(pdbfile,RESIN,BForg,BFs,goSM)
#myfile<-paste(PDBid,".dataGDV",sep="")
myfile<-"protein.dataGDV"
write.table(df, file=myfile, sep=",", eol="\n", append=TRUE, row.names=FALSE, col.names=FALSE)
# end write

# time
end_time2 <- Sys.time()
cat('Time Graphlet Degree Vector: \n')
print(end_time2-end_time1)
cat('Time-total: \n')
print(end_time2-start_time)

