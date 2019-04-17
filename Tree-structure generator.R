er###############################################################################################
#
#
# NPI activity “Stochastic Modelling of Atmospheric Re-entry Highly Energetic Break-up Events”
# Fragmentation model
#
#This file contains the function that generates the tree-structure of the parameter of the Dirichlet distribution of the masses of the fragments. The tree structure is a representation of the proposed alternative stick breaking process.
# To draw the resulting tree-structure uncomment the lines 73-79, 108-115.
#The input variables are defined in the main file.
#OUTPUT 
#masses:
#g.data= tree structure
#dir.vec.data= 
#alpha.data= splitting proportions vector
#vert.deg= auxiliary variable
#deg.ind= auxiliary variable
#adjm1= adjacency matrix
#adjm1.data= weighted adjacency matrix
#leaves= indices of the leaves



###############################################################################################
vert.deg<-vector() #                            #
deg.ind<-vector() #


Create.data <- function (deg1,deg2,n.leaves,height_tree,alpha.data,alpha) {
  
  if (n.leaves<=height_tree*(deg1-1)+1) 
    number_nodes<-n.leaves+height_tree #caterpillar tree-structure
  if(n.leaves>height_tree*(deg1-1)+1) #lobster tree-structure
  {
    ngroups<-(n.leaves-height_tree*(deg1-1)-1)/(deg2-1)
    res<-(n.leaves-height_tree*(deg1-1)-1)%%(deg2-1)
    if (res==1) res<-res+1
    number_nodes<-1+height_tree+(height_tree*(deg1-1))+
      floor(ngroups)*deg2+res
  }
  number_nodes<-number_nodes-1
  adjm1.data<-matrix(data=0,nrow=number_nodes,ncol=number_nodes) #adjacency matrix
  
  #Building the tree
  
  for (i in 1:(height_tree-1))
  {
    adjm1.data[i,i+1]<-1
  }
  
  for (i in (height_tree+1):min(number_nodes,deg1*
                                height_tree))
    adjm1.data[(i%%height_tree+1),i]<-1     
  
  if (number_nodes>deg1*height_tree)
  {
    ind<-sample((height_tree+1):(deg1*height_tree),
                ceiling((number_nodes-deg1*height_tree)/deg2)) 
  #choose the parents of the remaining fragments
    for (i in 1:(number_nodes-deg1*height_tree))
    {
      adjm1.data[ind[ceiling(i/deg2)],(deg1*height_tree+i)]<-1
    }
  }
  adjm1.data<-rbind(cbind(adjm1.data,rep(0,number_nodes)),
                    rep(0,number_nodes+1))
  
  #the main body at the end
  adjm1.data[height_tree,number_nodes+1]<-1
  number_nodes<-number_nodes+1
  
  g.data<-graph.adjacency(adjm1.data, mode="directed",weighted=TRUE)
  
###This code could be used to plot the tree#####################################################################
#   co<-layout.fruchterman.reingold.grid (g.data,flip.y=TRUE,circular=TRUE,root=1)
#   plot(g.data, layout=co,
#        vertex.color=c(rep("red",height_tree),rep("blue",(number_nodes-
#                                                            height_tree-1)),rep("red",number_nodes)),vertex.size=20, 
#        vertex.label.cex=1, vertex.label.color="white",edge.arrow.size=0.5)
#   #############################################################################################################
  
  
  if (n.leaves!=length(which(rowSums(adjm1.data)==0))) print("ERROR") #check 
  
  leaves<-which(rowSums(adjm1.data)==0) #vector of the fragments
  adjm1<-adjm1.data#store not weighted matrix 
  
  #memorize the structure of the tree
  for(j in 2:max(deg1,deg2))   {
    vert.deg<-c(vert.deg,which(rowSums(adjm1.data)==j)) 
    #all the nodes which are not leaves
    if(length(which(rowSums(adjm1.data)==j))>0)   
      deg.ind<-c(deg.ind,rep(j,length(which(rowSums(adjm1.data)==j)))) 
    #number of children of the nodes in vert.deg (all the nodes except the leaves)
  }
  
  #make the adjacency matrix a weigthed adiacency matrix
  for(j in 2:max(deg1,deg2))
  {
    for (i in vert.deg[which(deg.ind==j)])
      adjm1.data[i,which(adjm1.data[i,]==1)]<-c(alpha.data[1:(j-1)],1
                                                -sum(alpha.data[1:(j-1)]))
  }
  
  
  g.data<-graph.adjacency(adjm1.data, mode="directed",weighted=TRUE)
  
  
  ##Plot weighted tree###################################################################
  #   co<-layout.fruchterman.reingold.grid (g.data,flip.y=TRUE,circular=TRUE,root=1)
  #   plot(g.data, layout=co,
  #        vertex.color=c(rep("red",height_tree),rep("blue",(number_nodes-
  #                                                            height_tree-1)),rep("red",number_nodes)), vertex.size=10, 
  #        vertex.label.cex=0.8, vertex.label.color="white",edge.arrow.size=0.5,
  #        edge.label=round(E(g.data)$weight,digits=2))
  #######################################################################################
  
  dir.vec.data<-vector(mode="numeric",length=n.leaves)
  #dirichlet distribution alpha vector
  
  for (i in 1:n.leaves)
    dir.vec.data[i]<-prod(E(g.data)[get.shortest.paths(g.data,1,leaves[i],
                                                       output=c("both"))$epath[[1]]]$weight) 
  #generate alpha vector from the adjacency matrix
  dir.vec.data<-alpha*sort(dir.vec.data)
  
  
  repeat
  {
    masses<-rdirichlet(1,dir.vec.data) #check 
    print("Please wait")
    if (length(which(masses==0)) == 0)
      break
  }
  masses<-sort(masses)
  results<-list(masses,g.data,dir.vec.data,alpha.data,vert.deg,deg.ind,adjm1.data,leaves,adjm1)
  return(results)
}