###### ALL AT ONCE
GDVallatonce <- function(pdb,kb) {  
  xyz<-t(matrix(pdb$xyz,nrow=3)) #vector to Nx3 matrix
  cat('///////////////// ',dim(xyz)[1],'\n')
  dd<-10 # cutoff around protein for (graph) orbit count   
  Xmin<-min(xyz[kb,1])-dd # kb are biounit atoms
  Ymin<-min(xyz[kb,2])-dd
  Zmin<-min(xyz[kb,3])-dd
  Xmax<-max(xyz[kb,1])+dd
  Ymax<-max(xyz[kb,2])+dd
  Zmax<-max(xyz[kb,3])+dd
  IDX<-logical(length(pdb$atom$x)) 
  IDX[ (xyz[,1] > Xmin & xyz[,1] < Xmax) &
       (xyz[,2] > Ymin & xyz[,2] < Ymax) &
       (xyz[,3] > Zmin & xyz[,3] < Zmax) ]<-TRUE  
  IDX[kb]<-FALSE # da ni dvojno
  my.atoms<-rbind(xyz[kb,],xyz[IDX,])
  DDsave<-dist(my.atoms)
  DD<-DDsave*0
  threshold<-7.0 #treshold=7.0Å to create graph
  DD[DDsave<threshold]<-1
  adj<-as.matrix(DD)
  g<-graph_from_adjacency_matrix(adj,mode="undirected",weighted=NULL)
  g<-delete_vertices(g,which(degree(g)==0)) #delete isolated vertices
  rm(adj,DD)
  go<-count_orbits_per_node(g,max_graphlet_size=4) # orbit
  #smooth
  DD<-DDsave*0
  DD[DDsave<=2.1]<-1 # covalent bond cutoff (2.1A)for smoothing
  adj<-as.matrix(DD)
  goSM<-go*0 # for smoothed
  for (c in 1:length(kb)){
       n<-which(adj[c,]==1)
       if (length(n)>=2) {goSM[c,]<-(go[c,] + colMeans(go[n,]))*0.5}
            else if (length(n)==1) {goSM[c,]<-(go[c,] + go[n,])*0.5}
            else if (length(n)==0) {goSM[c,]<-go[c,]}
  }
  # end smooth
  goSM<-goSM[1:length(kb),] # smoothed GDV
  #
  return(goSM)
}
################################################


# PER PARTES ####################
GDVperpartes <- function(pdb,kb) {   
  dp<-10000 # the size of part
  L<-length(kb)
  N<-floor(L/dp) # split N+1 parts
  cat('Number of parts: ',N+1,'\n')
  ap<-1:dp
  al<-1:(L-N*dp) # the last part (L %% 1000)
  
  goSM<-numeric(0)
  for (r in 1:(N+1)) { 
      cat(r,'\n')
      if (r<=N){
          #cat(r,a1000+(r-1)*1000,'\n')
          myATOM<-ap+(r-1)*dp
          } else {
          #cat(r,al+(r-1)*1000,'\n')
          myATOM<-al+(r-1)*dp
      }  
  
      xyz<-t(matrix(pdb$xyz,nrow=3)) #vector to Nx3 matrix
      dd<-10 # cutoff around protein for (graph) orbit count      
      Xmin<-min(xyz[kb[myATOM],1])-dd
      Ymin<-min(xyz[kb[myATOM],2])-dd
      Zmin<-min(xyz[kb[myATOM],3])-dd
      Xmax<-max(xyz[kb[myATOM],1])+dd
      Ymax<-max(xyz[kb[myATOM],2])+dd
      Zmax<-max(xyz[kb[myATOM],3])+dd
      IDX<-logical(length(pdb$atom$x)) 
      IDX[ (xyz[,1] > Xmin & xyz[,1] < Xmax) &
           (xyz[,2] > Ymin & xyz[,2] < Ymax) &
           (xyz[,3] > Zmin & xyz[,3] < Zmax) ]<-TRUE
      IDX[kb[myATOM]]<-FALSE # da ni dvojno
      my.atoms<-rbind(xyz[kb[myATOM],],xyz[IDX,])
      DDsave<-dist(my.atoms)
      DD<-DDsave*0
      threshold<-7.0 #treshold=7.0Å to create graph
      DD[DDsave<threshold]<-1
      adj<-as.matrix(DD)
      g<-graph_from_adjacency_matrix(adj,mode="undirected",weighted=NULL)
      g<-delete_vertices(g,which(degree(g)==0)) #delete isolated vertices
      rm(adj,DD)
      goRES<-count_orbits_per_node(g,max_graphlet_size=4) # orbit
      #smooth
      DD<-DDsave*0
      DD[DDsave<=2.1]<-1 # covalent bond cutoff (2.1A)for smoothing
      adj<-as.matrix(DD)
      goRES_SM<-goRES*0 # for smoothed
      for (c in 1:length(myATOM)){
          n<-which(adj[c,]==1)
          if (length(n)>=2) {goRES_SM[c,]<-(goRES[c,] + colMeans(goRES[n,]))*0.5}
          else if (length(n)==1) {goRES_SM[c,]<-(goRES[c,] + goRES[n,])*0.5}
          else if (length(n)==0) {goRES_SM[c,]<-goRES[c,]}
      }
      #end smooth
      goRES_SM<-goRES_SM[1:length(myATOM),] # only selected atoms "per partes"
      goSM<-rbind(goSM,goRES_SM)
      cat(myATOM[1],myATOM[length(myATOM)],'\n')
  }
  return(goSM)
}
###########################


