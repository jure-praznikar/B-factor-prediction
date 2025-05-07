### scale (per-chain) B-factors 
BFscale <- function(pdb,CHAIN,flagLOG) {
   indsTRIM <- atom.select(pdb,chain=CHAIN)     
   pdb<-trim.pdb(pdb, inds = indsTRIM)
   inds <- atom.select(pdb) # re-index
   IDX<-matrix(nrow=length(pdb$atom$chain),ncol=length(CHAIN))
   IDX[,]<-FALSE
   if (length(CHAIN)>1) {
     for (ch in 1:length(CHAIN)) {
       if (ch<=9){chnum<-paste('0',ch,sep='')}else{chnum<-ch} # to have name _01 _02 _10 _33
       var_name <- paste('inds_CH_',chnum,sep='') # create var name
       indsCHtemp<-atom.select(pdb,chain=CHAIN[ch])
       #cat(CHAIN[ch],'\n')
       assign(var_name,indsCHtemp) # set given var name to atom select
       # idx's to merge chain scaled BF toghether it could be A C and D F
       tempidx<-which(pdb$atom$chain==CHAIN[ch])
       IDX[tempidx,ch]<-TRUE
      }
   }

   flagDIFFbf<-FALSE
   # "significant" difference of mean Bf , then flagDIFFbf<-TRUE
   # TWO chains   
   if (length(CHAIN)==2) {
       ratioL<-max(length(inds_CH_01$atom),length(inds_CH_02$atom)) /
               min(length(inds_CH_01$atom),length(inds_CH_02$atom))
       Bx<-pdb$atom[inds_CH_01$atom,"b"]
       By<-pdb$atom[inds_CH_02$atom,"b"]
       ratioB<-max(mean(Bx),mean(By)) / min(mean(Bx),mean(By))
       dB<-abs(mean(Bx)-mean(By))
       cat('Chain length ratio: ',round(ratioL,3),' \n')
       cat('Bf mean ratio: ',round(ratioB,3),' \n')
       cat('Bf mean diff: ',round(dB,3),' \n')
       if (ratioL<1.5 & ratioB>1.25 & dB>5){# large difference betweeen chains, do scale per chain
          flagDIFFbf<-TRUE
          cat("!!!!!! Enable scaling per chain !!!!!!!! (2) \n")
       } 
   }   
   ### THREE chains
   if (length(CHAIN)==3) {
      r<-0
      ratio3L<-numeric(0)
      ratio3B<-numeric(0)
      dB3<-numeric(0)  
      for (i in 2:length(CHAIN)) { # Columns > from 2 to N
           r<-r+1
           for (j in 1:r) {
               cat(CHAIN[j],CHAIN[i]," >> ")
               tempCHindsX<-get(ls(pattern='inds_CH_')[j]) 
               tempCHindsY<-get(ls(pattern='inds_CH_')[i])
               Bx<-pdb$atom[tempCHindsX$atom,"b"]
               By<-pdb$atom[tempCHindsY$atom,"b"]
               ratioL<-max(length(Bx),length(By))/min(length(Bx),length(By))
               ratio3L<-c(ratio3L,ratioL)
               ratioB<-max(mean(Bx),mean(By)) / min(mean(Bx),mean(By))
               ratio3B<-c(ratio3B,ratioB)
               dB<-abs(mean(Bx)-mean(By))
               dB3<-c(dB3,dB)
               cat('Chain length ratio: ',round(ratioL,3),' \n')
               cat('Bf mean ratio: ',round(ratioB,3),' \n')
               cat('Bf mean diff: ',round(dB,3),' \n')
           }    
      }

          cat('ratio3L ',ratio3L,sum(ratio3L<1.5),'\n')
          cat('ratio3B ',ratio3B,sum(ratio3B>1.25),'\n')
          cat('dB3 ',dB3,sum(5>dB3),'\n')
      if (sum(ratio3L<1.5)>=3  & #large difference betweeen chains, do scale per chain
          sum(ratio3B>1.25)>=1 &
          sum(5>dB3)>=1       ){
          flagDIFFbf<-TRUE
          cat("!!!!!! Enable scaling per chain !!!!!!!! (3)\n")
      } 
   }

   ### perform scaling per chain (if two or three chains exist!! and flagDIFFbf<-TRUE )
   if (length(CHAIN)>1 & length(CHAIN)<=3 & flagDIFFbf==TRUE) {
      BF<-numeric(length=length(pdb$atom$chain))
      BF[]<-NA
      cat('Scale per chain. \n')
      for (ch in 1:length(CHAIN)) {
        tempCHinds<-get(ls(pattern='inds_CH_')[ch]) # get var name: this is "atom.select" 
        BFch<-pdb$atom[tempCHinds$atom,"b"]
        cat(ls(pattern='inds_CH_')[ch],CHAIN[ch],mean(BFch),sd(BFch),
            length(BFch),round(length(BFch)/7),'\n') # round(length(BFch)/7) estimated num. of CA atoms
        if (flagLOG) { 
            BFch<-scale(log(BFch))# log and scale per chain !!!!!!!!!!
        } else {  
            BFch<-scale(BFch)# scale per chain (but no log-transform) !!!!!!!!!!
        }    
        BF[IDX[,ch]]<-BFch  
      }
      # check correlation between raw Bfs and scaled per-chain, it should be "high" but
      # not = 1 !!!
      if ( sum(is.na(BF))>0 ) {
          cat('!!!!!!!!!!!!!!!!!!!! NOT sum(is.na(BF))==0 \n') 
          Sys.sleep(1.0E6)
      }
      tempcc<-cor(pdb$atom[inds$atom,"b"],BF,method='spearman')
      cat('>>>>>>>>> Correlation test after scaling (2 or 3 chains): ',tempcc,'\n')
      
   } else {
      cat('Scale as single chain. \n')
      BF<-pdb$atom[inds$atom,"b"]
      cat('Mean Bf=',mean(BF),' stdev(Bf)=',sd(BF),'\n')
      if (flagLOG) {
          BF<-scale(log(BF))# log-transform and scale !!!!!!!!!!
          cat('<<CHAIN: scale and log>>\n')
      } else {   
          BF<-scale(BF)# only scale !!!!!!!!!!
          cat('<<CHAIN: scale (not log)>>\n')
      }   
      tempcc<-cor(pdb$atom[inds$atom,"b"],BF,method='spearman')
      cat('>>>>>>>>> Correlation test after scaling (single chain): ',tempcc,'\n')
   }  
   ### 
  return(BF) # scaled B-factors
}
