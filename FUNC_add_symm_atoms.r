
### symm function >> ADD close contact atoms
addSMYTatoms <- function(pdb,pdbfile,flagLIGNUC) {
# lig+nuc ########################################################################## ratio 0.1 >> 10%
      if (flagLIGNUC) {
      ## set ligand and nucleic acid to chainID "9" > close contacts
         indsNOT1 <- atom.select(pdb, "ligand")
         indsNOT2 <- atom.select(pdb, "nucleic")
         indsNOTprot <- combine.select(indsNOT1, indsNOT2, operator="+")
         indsPROT<-atom.select(pdb,"protein") # to calculate ratio
         #ratio1 > compare number of protein and lig-nuc atoms 
         ratio1<-length(indsNOTprot$atom)/length(indsPROT$atom)
         # ratio2 calculate contact residues
         indsPROTcb<-atom.select(pdb,"protein",elety="CB")
         Pcb<-t(matrix(pdb$xyz[indsPROTcb$xyz],nrow=3))
         Plig<-t(matrix(pdb$xyz[indsNOTprot$xyz],nrow=3))
         MAT <- as.matrix(pdist(Pcb,Plig)) #distance between CB and hetero atoms
         v<-which(MAT<5,arr.ind = TRUE) # if less than 5, is contact
         v<-unique(v[,1])
         ratio2<-length(v)/length(indsPROTcb$atom)
         cat('Ratio1-number  (ligand&|nuc)/prot: ',round(ratio1*100,1),' %\n')
         cat('Ratio2-contact (ligand&|nuc)/prot: ',round(ratio2*100,1),' %\n')
         #if (length(indsNOTprot$atom)>0){
         if (ratio1>=0.1 & ratio2>=0.1){ # ratio 0.1 >> 10%
           cat('>> Ligand and nucleic atoms: ',length(indsNOTprot$atom),'\n')
           pdbNOTprot<-trim.pdb(pdb,inds=indsNOTprot)
           pdbNOTprot$atom$chain<-"9"
         } # later append PDBlist<-list(), see below
      ## end set ligand and nucleic acid to chainID "9" > close contacts
      }    
# lig+nuc ##########################################################################
      # select protein atoms
      idx<-atom.select(pdb,"protein")
      pdb<-trim.pdb(pdb,inds=idx)
      pdbORG<-pdb # store for later use to write final pdb file > write.pdb
      # CLEAN it > to have a unique residue number overall chains
      pdb<-clean.pdb(pdb, consecutive = TRUE,force.renumber = TRUE)
      pdb$atom$chain[1:length(pdb$atom$chain)]<-"9" # set chain name to "9" >> discriminate from non-symm resi.
      # select CB atoms - the first stage of searching close contacts
      idxCB<-c(which(pdb$atom$elety=="CB"),which(pdb$atom$resid=="GLY" & pdb$atom$elety=="CA"))
      # geometric center
      gc<-colMeans(pdb$atom[,c("x","y","z")])
      ## max deviation in x, y, z
      dxyz<-numeric(3)
      dxyz[1]<-max(abs(gc[1]-pdb$atom[,c("x")]))
      dxyz[2]<-max(abs(gc[2]-pdb$atom[,c("y")]))
      dxyz[3]<-max(abs(gc[3]-pdb$atom[,c("z")]))
      # d is max distance between the geometric center and the most distant atom
      d<-sqrt(sum(dxyz^2))
      #
      VdW<-3 # approx. two times Van der Waals radii (3Å)
      cutoff<-5+VdW # distance = 5Å
      # if the distance between the centers of the two molecules is
      # more than 2*d+cutoff, then are not in contact
      dc<-2*d+cutoff
      #Cbeta cutoff - to search nearby residues CB-CB is less than 15Å
      CBcutoff<-15 # 15Å
      #CBcutoff<-7.5 # 7.5Å
      cat('*** CB cutoff (sym contact shell): ',CBcutoff,' ***\n')
      #
      xyz0 <- rbind(matrix(pdb$xyz, nrow=3), 1) # initial coordinates as vector
      xyz30<-t(matrix(pdb$xyz,nrow=3)) # initial coordinates as Nx3 matrix
      # extract symmetry data and Cell size (the data must be in PDB header)
      temp<-readLines(pdbfile,n=10000) # read first 10000 lines of header
      k<-grep('CRYST1',temp)
      # if k>1 example ARK 247  THAT CRYST1 AND SCALE REC
      k<-k[length(k)] # take last one    
      #12345678901234567890123456789012345678901234567890123456789012345678901234567890
      #CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8 
      temp<-substring(temp[k], first=1, last=65)
      A<-unlist(strsplit(trimws(temp), "\\s+"))  
      CELL<-as.numeric(A[2:4])
      ANGLE<-as.numeric(A[5:7])
      ## find SMTRY else use cry-lib, or find MTRX (cryEM)
      temp<-readLines(pdbfile,n=10000)
      k1<-grep('REMARK 290   SMTRY1',temp)
      k2<-grep('REMARK 290   SMTRY2',temp)
      k3<-grep('REMARK 290   SMTRY3',temp)
      k<-sort(c(k1,k2,k3))
      x<-length(k)/3
      #x%%1==0 #to check the fractional part (is fractional part zero)      
      if (x>0 & x%%1==0 & sum(CELL)>3.01 ) {
          #get SMTRYn from PDB header
          M<-getSYMfromPDBheader(temp,k) #function
          # set dummy center / as list 
          C<-list(matrix(0,ncol=1,nrow=3))
          ### end REMARK 290   SMTRY
      } else if (sum(CELL)>3.01) { 
          # CRY lib###
          MC<-getSYMuseCRYlib(A,ANGLE,CELL) #function
          M<-MC[[1]] # rotation and translation
          C<-MC[[2]] # center (needed to construct all sym)
          # end CRY lib ####
          cat('Total symmetries: ',length(M)*length(C),'\n')
          # end CRY ####
      } else if  (abs(sum(CELL)-3)<1.0E-3 & 
                  length(grep('MTRIX1',temp))>=1  &
                  length(grep('MTRIX2',temp))>=1  &
                  length(grep('MTRIX3',temp))>=1   ) { #other than X-ray crystallography, for example cryEM      
          M<-getMTRIXfromPDBheader(temp) #temp<-readLines(pdbfile,n=10000)
          # set cell lengths >> all equal!
          abc<-numeric(0)
          for (m in 1:length(M) ){
              abc<-c(abc,M[[m]][,4]) # column 4
          }
          abc<-abc[abc>1.0E-5] # remove zeros
          if ( sd(abc)<1.0E-5 ) {
               CELL[1]<-abc[1]
               CELL[2]<-abc[1]
               CELL[3]<-abc[1]    
              }
      
          # set dummy center / as list 
          C<-list(matrix(0,ncol=1,nrow=3))             
      }
# END ## find SMTRY else use cry-lib, or find MTRX (example is cryEM)    
      a<-deg2rad(ANGLE[1]) # alpha
      b<-deg2rad(ANGLE[2]) # beta
      g<-deg2rad(ANGLE[3]) # gamma
      temp<-(cos(b)*cos(g)-cos(a))/(sin(b)*sin(g))
      a_st<-acos(temp) #alpha "star", equation (7)
      # convert CELL to orthogonal 
      # https://www.ucl.ac.uk/~rmhajc0/frorth.pdf
      #  equation 14,15, 16, from #https://www.ucl.ac.uk/~rmhajc0/frorth.pdf
      # X <- xi*CELL[1] + yi*CELL[2]*cos(g) + zi*CELL[3]*cos(b)
      # Y <- yi*CELL[2]*sin(g) - zi*CELL[3]*sin(b)*cos(a_st)
      # Z <- zi*CELL[3]*sin(b)*sin(a_st)     
      mat<-matrix(0,ncol=4,nrow=3)
      PDBlist<-list() # store contacts in list
# lig+nuc ##########################################################################
      if (flagLIGNUC) {
      ## append ligand and/or nuleic acid atoms   
         if (ratio1>=0.1 & ratio2>=0.1){ # ratio 0.1 >> 10%
             PDBlist<-append(PDBlist,list(pdbNOTprot)) # append ligand and nucleic
         }
      }    
      ## end append ligand and nuleic
# lig+nuc ##########################################################################
      # -2,-1,0,1,2 all combinations in x, y,z
      myshift<-seq(-5,5) # shift along each axis (only translation)
      if  (abs(sum(CELL)-3)<1.0E-3) {myshift<-0} #other than X-ray crystallography, for example cryEM
      for (i in 1:length(M)) {
           mat<-M[[i]]
           for (xi in myshift) {
             for (yi in myshift) {
               for (zi in myshift) {
                 # to orthogonal  (translation vector) 
                 mat[1,4] <- M[[i]][1,4] + ( xi*CELL[1] + yi*CELL[2]*cos(g) + zi*CELL[3]*cos(b) )
                 mat[2,4] <- M[[i]][2,4] + ( yi*CELL[2]*sin(g) - zi*CELL[3]*sin(b)*cos(a_st) )
                 mat[3,4] <- M[[i]][3,4] + ( zi*CELL[3]*sin(b)*sin(a_st) )
                 m1<-mat[1,4] 
                 m2<-mat[2,4]
                 m3<-mat[3,4]
                 # CENTER shift - already inlcuded in PDB header (but not in cryst_symm >> construction needed)  
                 for (c in 1:length(C)) {
                     Cx<- CELL[1]*C[[c]][1,1] + CELL[[2]]*cos(g)*C[[c]][2,1] + CELL[3]*cos(b)*C[[c]][3,1] 
                     Cy<- CELL[2]*sin(g)*C[[c]][2,1] - CELL[3]*sin(b)*cos(a_st)*C[[c]][3,1]
                     Cz<- CELL[3]*sin(b)*sin(a_st)*C[[c]][3,1]
                     mat[1,4] <- m1 + Cx
                     mat[2,4] <- m2 + Cy
                     mat[3,4] <- m3 + Cz 
                     #init pdb
                     pdb1<-pdb
                     pdb2<-pdb                     
                     #xyz <- matrix(mat %*% xyz0, nrow = 1)
                     pdb2$xyz<-matrix(mat %*% xyz0, nrow = 1)
                     pdb2ALL<-pdb2
                     xyzT<-t(matrix(pdb2$xyz,nrow=3))
                     XYZgc<-colMeans(xyzT)
                     test<-dist(rbind(gc,XYZgc))
                     #if the distance between the centers of the two molecules is more than 2*d+cutoff >> no contact
                     if ( test>0.5 & test<dc ) {
                        MAT <- as.matrix(pdist( xyz30[idxCB,], xyzT[idxCB,] ))
                        #v<-which(MAT<=(cutoff+7), arr.ind = TRUE) # CB atoms cutoff+7Å >> 15Å
                        v<-which(MAT<=CBcutoff, arr.ind = TRUE) # CBcutoff=15Å
                        #v1<-unique(v[,1])
                        v2<-unique(v[,2]) # column2 CB atoms
                        if (length(v2)>=1){
                            #pdb1<-trim.pdb(pdb1,resno=pdb1$atom$resno[idxCB[v1]])
                            pdb2<-trim.pdb(pdb2,resno=pdb2$atom$resno[idxCB[v2]])
                            PDBlist<-append(PDBlist,list(pdb2))
                        }
                     }
                 }    
             }   
            }
           }         
      }
      #### Write/create pdb file
      if (length(PDBlist)>=1) { 
        for ( f in 1:length(PDBlist)) {
            # add/append PDBlist >> symmetry atoms
            suppressWarnings(pdbORG <- cat.pdb(pdbORG,PDBlist[[f]],renumber=TRUE,rechain=FALSE)) 
        }
      } else {pdbORG<-clean.pdb(pdbORG, consecutive = TRUE,force.renumber = TRUE)}
 return(pdbORG)
}   
#########################
