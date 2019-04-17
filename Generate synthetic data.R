###############################################################################################
#
#
# NPI activity “Stochastic Modelling of Atmospheric Re-entry Highly Energetic Break-up Events”
# Fragmentation model
#
#Goal: generate not complete synthetic observations
#The output is a set of sorted masses Dirichlet distributed              
#The parameter of the Dirichlet distribution depends on the structure of tree  
#and on the splitting parameter.                                               
#                                
#################################################################################################


while(length(which(alpha.data==0))>0){
  alpha.data<-rdirichlet(1, rep(alpha,max(deg1,deg2))) #generation of the splitting proportions vector 
}


masses<-function(n.leaves){
  return(Create.data(deg1,deg2,n.leaves,height_tree,alpha.data,alpha))}
n<-rep(n.leaves,k)
masses_list<-lapply(n,masses) #complete set of observations
truealphavalues<-masses_list[[1]][[3]]
truealphavalues2<-masses_list[[2]][[3]]

fobserved<-function(nu){
  return(rbern(nu,prob_observed))
}
observed<-apply(cbind(n),1,fobserved)

n_observed<-function(nk)
{  return(sum(observed[,nk]))
}
s<-apply(cbind(1:k),1,n_observed)



ma_observed<-function(nk)
{
  return(masses_list[[nk]][[1]][which(observed[,nk]==1)])
}

masses_observed<-apply(cbind(1:k),1,ma_observed) #final output of this function: not complete observations of masses



