
### check Space Group, cryst_symm(SG)
SGcheck <- function(SG){
    tryCatch(
        expr = {
            flag<-TRUE
            crsym <-cryst_symm(SG)
            return(flag)
        },
        error = function(e){
            flag<-FALSE
            message('Space Group Error!!')
            print(e)
            return(flag)
        },
        warning = function(w){
            flag<-FALSE
            message('Space Group Error!!')
            print(w)
            return(flag)
        }
    )
    #return(flag)
}
#

### function convert deg2grad
deg2rad <- function(deg) {(deg * pi) / (180)}

###
# Converting fractional coordinates into cartesian coordinates for crystallography
# Write down the symmetry operation in fractional coordinates.
# Determine the transformation matrix F that converts fractional to Cartesian coordinates
# for the given unit cell parameters.
# Invert the transformation matrix F.
# Apply the formula Mcartesian = F ⋅ Mfractional ⋅ F^(−1) to obtain the Cartesian symmetry matrix.

Frac_to_cartesian<-function(ANGLE,CELL,Mf) {
     F<-matrix(ncol=3,nrow=3) #F is the fractional-to-Cartesian transformation matrix
     a<-deg2rad(ANGLE[1]) # alpha
     b<-deg2rad(ANGLE[2]) # beta
     g<-deg2rad(ANGLE[3]) # gamma
#
##    | CELL[1] CELL[2]*cos(g) CELL[3]*cos(b)              |
## F= |   0     CELL[2]*sin(g) CELL[3]*n2                  |
##    |   0          0         CELL[3]*sqrt(sin(b)^2-n2^2) |
#     
     n2<-(cos(a)-cos(g)*cos(b))/sin(g)
     F[1,1]<-CELL[1] # cell length "a"
     F[1,2]<-CELL[2]*cos(g) # length "b" * cos(gamma)
     F[1,3]<-CELL[3]*cos(b) # length "c" * cos(beta)
     F[2,1]<-0
     F[2,2]<-CELL[2]*sin(g) # length "b" * sin(gamma)
     F[2,3]<-CELL[3]*n2     # length "c" * n2
     F[3,1]<-0
     F[3,2]<-0
     F[3,3]<-CELL[3]*sqrt(sin(b)^2-n2^2) # # length "c"*sqrt(sin(b)^2-n2^2)
     # to obtain the Cartesian symmetry matrix
     Mc<-F %*% Mf %*% ginv(F)  # Mc = F * Mf * F^(-1)
     Mc[abs(Mc)<1.0E-8]<-0
     return(Mc)     
}
#END ############# Converting fractional coordinates into cartesian coordinates

### check SYM in inputfile
## symmetry matrix can be extracted from PDB header file or
## using CELL data and cryst_symm function
checkSYM<-function(pdbfile){
      flagSYM<-FALSE
      temp<-readLines(pdbfile,n=10000)
      k<-grep('CRYST1',temp)
      # if k>1 example ARK 247  THAT CRYST1 AND SCALE REC
      k<-k[length(k)] # take last one    
      #12345678901234567890123456789012345678901234567890123456789012345678901234567890
      #CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8 
      temp2<-substring(temp[k], first=1, last=65)
      A<-unlist(strsplit(trimws(temp2), "\\s+"))     
      CELL<-as.numeric(A[2:4])
      ANGLE<-as.numeric(A[5:7])
      # Space Group
      if (length(A)==9){SG<-paste(A[8],A[9])}  # Space Group
      if (length(A)==10){SG<-paste(A[8],A[9],A[10])}  # Space Group
      if (length(A)==11){SG<-paste(A[8],A[9],A[10],A[11])}  # Space Group
      #https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html
      #Example of experimental method other than X-ray crystallography or fiber diffraction
      #CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1 
      # sum(CELL) should not be 3.0
      if (SGcheck(SG)){
          crsym <- cryst_symm(SG)
          if (!is.na(sum(CELL)) & !is.na(sum(ANGLE)) & is.character(crsym$SG) & sum(CELL)>3.01 ){
          cat(crsym$SG,' Found symmetry. \n')
          flagSYM<-TRUE
          #print(crsym)
          } else if (abs(sum(CELL)-3)<1.0E-3 & 
            length(grep('MTRIX1',temp))>=1  &
            length(grep('MTRIX2',temp))>=1  &
            length(grep('MTRIX3',temp))>=1   ) { #other than X-ray crystallography, like EM
              cat('Example of experimental method other than X-ray crystallography or fiber diffraction \n')
              cat(A, '\n')
              k1<-grep('MTRIX1',temp)
              k2<-grep('MTRIX2',temp)
              k3<-grep('MTRIX3',temp)
              k<-sort(c(k1,k2,k3))
              x<-length(k)/3
              if ( x>0 & x%%1==0 ) { #x%%1==0 #to check the fractional part (is fractional part zero)
                cat('>> ',x,' matrices (MTRIXn) were found.\n')   
                flagSYM<-TRUE
              }

          } else {
            cat('No MTRIX and/or symmetry was found! \n')
          }
          
      }    
      return(flagSYM)
}

### get SYM matrix from PDB header
getSYMfromPDBheader<-function(temp,k) {
    cat('REMARK 290   SMTRY from the PDB header\n')
    temp<-temp[k]
    #print(temp)
    s1<-seq(1,length(temp),3)
    # create crystallographic  matrices (rotation and translation) - data from PDB header
    ## REMARK 290   SMTRY in PDB are are already cartesian coordinates
    M<-list()
    for (i in 1:(length(temp)/3)){
         mat0<-matrix(0,ncol=4,nrow=3)
         k<-i*3-2
         #1
         A<-unlist(strsplit(temp[k],split = ' '))
         A<-A[A != ""]
         mat0[1,]<-as.numeric(A[5:8])
         #2
         A<-unlist(strsplit(temp[k+1],split = ' '))
         A<-A[A != ""]
         mat0[2,]<-as.numeric(A[5:8])
         #3
         A<-unlist(strsplit(temp[k+2],split = ' '))
         A<-A[A != ""]
         mat0[3,]<-as.numeric(A[5:8])
         # append
         M<-append(M,list(mat0))
    }
    return(M)
}

### create symm matrix, use CRY lib and cell size, angle
getSYMuseCRYlib<-function(A,ANGLE,CELL) {
    # CRY ###
    # Space Group
    cat('SMTRY using cry-library \n')
    if (length(A)==9){SG<-paste(A[8],A[9])}  # Space Group
    if (length(A)==10){SG<-paste(A[8],A[9],A[10])}  # Space Group
    if (length(A)==11){SG<-paste(A[8],A[9],A[10],A[11])}  # Space Group
    crsym <- cryst_symm(SG) #symmetry operations
    # C 2 2 2 Error !!!!!!!!!!!!!!!! error in cryst_symm ?????
    if (SG=="C 2 2 2"){
       crsym$T[[2]][3,1]<-0 
       crsym$T[[4]][3,1]<-0
    }
    # C 2 2 2 Error !!!!!!!!!!!!!!!!    
    SN<-length(unlist(crsym$PG))/9 # all unlisted elemnts over 9 (3x3) = number of matrices
    #print(crsym)
    M<-crsym$PG #!!! this are fractional coordinates (in PDB header are already cartesian)
    # Converting fractional coordinates into cartesian coordinates
    for (i in 1:length(M)) {
         M[[i]]<-Frac_to_cartesian(ANGLE,CELL,M[[i]]) # function Frac_to_cartesian
    }
    # create crystallographic  matrices (rotation and translation) M=3x4    
    for (i in 1:SN) {
        #cbind(M[[1]],unlist(crsym$T[1])*CELL)
        M[[i]]<-cbind(M[[i]],unlist(crsym$T[i])*CELL)
    }
    # center
    C<-crsym$C
    MC<-list(M,C)
    return(MC)    
}

### get MATRX from PDB header !!!! other than X-ray crystallography, like EM !!!!
getMTRIXfromPDBheader<-function(temp) {
    #temp<-readLines(pdbfile,n=10000)
    k1<-grep('MTRIX1',temp)
    k2<-grep('MTRIX2',temp)
    k3<-grep('MTRIX3',temp)
    k<-sort(c(k1,k2,k3))
    #x<-length(k)/3
    cat('MTRIX tranformation (close contacts)\n')
    temp<-temp[k]
    #print(temp) # print on screen
    s1<-seq(1,length(temp),3)
    # create matrices (rotation and translation) - data from PDB header
    M<-list()
    abc<-numeric(0)
    for (i in 1:(length(temp)/3)){
         mat0<-matrix(0,ncol=4,nrow=3)
         k<-i*3-2
         #1
         A<-unlist(strsplit(temp[k],split = ' '))
         A<-A[A != ""]
         mat0[1,]<-as.numeric(A[3:6])
         #2
         A<-unlist(strsplit(temp[k+1],split = ' '))
         A<-A[A != ""]
         mat0[2,]<-as.numeric(A[3:6])
         #3
         A<-unlist(strsplit(temp[k+2],split = ' '))
         A<-A[A != ""]
         mat0[3,]<-as.numeric(A[3:6])
         # append
         M<-append(M,list(mat0))
         #abc<-c(abc,mat0[(mat0[,4]>0.1),4]) # append not zero translations
    }
    return(M)
}
