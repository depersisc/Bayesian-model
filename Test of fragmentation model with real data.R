###############################################################################################
#
# NPI activity “Stochastic Modelling of Atmospheric Re-entry Highly Energetic Break-up Events”
# Fragmentation model
#
#This code contains a test of the Fragmentation model when we used a set of real data coming from the
#ATV 1 re-entry observation campaign, assuming a total number of fragments equal to 20 and that the predicted masses sum to 1000.
################################################################################################
rm(list=ls(all=TRUE))
library(igraph)
library(xtable)#to generate latex tables
library(MASS) 
library(MCMCpack)
library(coda)
library(Rlab)

source("Tree-structure generator.R")
################################################################################
###Input tree-settings data#####################################################
height_tree<-3                                                     #
deg1<-20                                                                       #
deg2<-5       
alpha<-1000
n.leaves<-20 #number of fragments      
mu<-0
while(length(which(mu==0))>0){
mu<-rdirichlet(1, rep(alpha,max(deg1,deg2)))
}
k=1 #number of observations 
##################Read data from a *.csv file#####################################################
s<-vector(length=k)
masses_observed<-read.csv("realdata.csv")$Mass #known data
totalmass<-1000 #test value
for (i in 1:k)
  s[k]<-length(masses_observed) #in this case only one
masses_observed<-masses_observed/1000
masses_observed<-sort(masses_observed)
n<-n.leaves 
z<-vector("list", k)
#Gibbs sampler
################################################################################
#Initial values
################################################################################
startingvalues<-Create.data(deg1,deg2,n,height_tree,mu,alpha)
vec.alpha<-startingvalues[[3]] #alpha vector
lbeta_alpha<-sum(lgamma(vec.alpha))-lgamma(sum(vec.alpha))
adjm1.weighted<-startingvalues[[7]] #weighted adjacency matrix (tree)
adjm1<-startingvalues[[9]] # adjacency matrix (tree)
vert.deg<-startingvalues[[5]]
deg.ind<-startingvalues[[6]]
leaves<-startingvalues[[8]]
creatervector<-function(nu){
  r<-vector(mode="integer",length=n)
  r[sample(1:n, s[nu], replace = FALSE, prob = NULL)]<-1
  return(r)
}
r<-lapply(cbind(1:k),creatervector) #indicator variable (a list of vectors, each of them of size s_j, with 1 and 0. 1 for observed masses, 0 missing masses)

accept_mu<-0
accept_r<-0
accept.tree<-0
niterations<-1000
memorize.vec.alpha<-vec.alpha
memorize.mu<-mu
memorize.z<-z[[1]]

for (iterations in 1:niterations){
  ################################################################################
  #1)sample each z_j independently, z is a list, z[[1]]=z_1
  #needed variables: vec_alpha, r
  #each row for each observation, vec_alpha's length is n
  for (i in 1:k){
    indices<-which(r[[i]]==0)
    z[[i]]<-0
    while(length(which(z[[i]]==0))>0){
      new_masses<-rdirichlet(1,vec.alpha[indices])
      new_masses_normalized<-new_masses*(1-sum(masses_observed[[i]]))/sum(new_masses)
      z[[i]]<-sort(new_masses_normalized)
    }
  } 
  #############################################################################################################################
  #2)Sample each r_j. Metropolis
  #Propose a new candidate keeping constant the number of 1
  r.new<-lapply(cbind(1:k),creatervector)
  #Calculate the acceptance probability
  #Accept probability
  prob_acc<-runif(1,0,1)
  acceptance_prob_num<-0
  acceptance_prob_den<-0
  
  for (i in 1:k){
    acceptance_prob_num<-sum(acceptance_prob_num,log(masses_observed[[i]])*(vec.alpha[which(r.new[[i]]==1)]-1), log(z[[i]])*(vec.alpha[which(r.new[[i]]==0)]-1))
    acceptance_prob_den<-sum(acceptance_prob_den, log(masses_observed[[i]])*(vec.alpha[which(r[[i]]==1)]-1), log(z[[i]])*(vec.alpha[which(r[[i]]==0)]-1))
  }
  acceptance_prob<-exp(acceptance_prob_num-acceptance_prob_den)
  if(prob_acc<min(1, acceptance_prob)) 
  {
    r<-r.new
    accept_r<-accept_r+1
  }
  #the new candidate is accepted
  ########################################################################################
  #3)Sample mu splitting proportions vector
  mu.new<-0
  while(length(which(mu.new==0))>0){
    #mu.new<-rdirichlet(1, rep(alpha,max(deg1,deg2))) #generate a new splitting parameters vector # check if there are zeros
    #mu.new<-rdirichlet(1, alpha*rep(seq(1,max(deg1,deg2)/2,length=4),5))
    mu.new<-rdirichlet(1, alpha*rep(c(1,4),10))
  }
  
  adjm1.weighted.new<-adjm1.weighted #new adiacency matrix to update
  
  #update the weighted adjacency matrix
  for(j in 2:max(deg1,deg2)){
    for (i in vert.deg[which(deg.ind==j)])
      adjm1.weighted.new[i,which(adjm1.weighted.new[i,]!=0)]<-c(mu.new[1:(j-1)],1-sum(mu.new[1:(j-1)]))
  }
  
  #update the tree
  g.new<-graph.adjacency(adjm1.weighted.new, mode="directed",weighted=TRUE)
  
  vec.alpha.new<-vector(mode="numeric",length=n.leaves)#dirichlet distribution alpha vector
  for (i in 1:n)
    vec.alpha.new[i]<-prod(E(g.new)[get.shortest.paths(g.new,1,leaves[i],output=c("both"))$epath[[1]]]$weight)   #generate alpha vector from the adjacency matrix
  vec.alpha.new<-alpha*sort(vec.alpha.new)
  #update the log-likelihood (and then the acceptance probability), considering the already accepted new z_i and new r
  
  llikelihood<-0
  for (i in 1:k){
    partllikelihood<-sum(llikelihood, (vec.alpha[which(r[[i]]==1)]-1)*log(masses_observed[[i]]), log(z[[i]])*(vec.alpha[which(r[[i]]==0)]-1) )
  }
  llikelihood<-partllikelihood-lbeta_alpha
  
  #update the log-likelihood (and then the acceptance probability), considering the proposed mu.new
  lbeta_alpha.new<-sum(lgamma(vec.alpha.new))-lgamma(sum(vec.alpha.new))
  llikelihood.new<-0
  for (i in 1:k){
    partllikelihood.new<-sum(llikelihood.new, (vec.alpha.new[which(r[[i]]==1)]-1)*log(masses_observed[[i]]), log(z[[i]])*(vec.alpha.new[which(r[[i]]==0)]-1) )
  }
  llikelihood.new<-partllikelihood.new-lbeta_alpha.new
  
  #Accept probability
  prob_acc<-runif(1,0,1)
  if(prob_acc<min(1,exp(llikelihood.new-llikelihood)*ddirichlet(mu,rep(1,max(deg1,deg2)))/ddirichlet(mu.new,rep(1,max(deg1,deg2))))) 
  {
    vec.alpha<-vec.alpha.new 
    lbeta_alpha<-lbeta_alpha.new
    mu<-mu.new
    adjm1.weighted<-adjm1.weighted.new #change the weights of the matrix but not the structure
    llikelihood<-llikelihood.new
    accept_mu<-accept_mu+1
  }
  ################################################################################################################
  ########################################################################################
  #4)Sample the tree structure
  #change tree move a leaf
  selection<-sample(leaves,1)
  sel.newfather<-which(rowSums(adjm1)!=0)
  father<-which(adjm1[,selection]!=0)
  
  adjm1.st<-adjm1
  adjm1.new<-adjm1
  
  adjm1.new[father,selection]<-0
  feasible_nodes.old<-0
  for (i in sel.newfather){
    if (length(which(adjm1.new[i,]!=0))<deg1)
      feasible_nodes.old<-feasible_nodes.old+1
  }
  
  
  exit<-"FALSE"
  while(exit=="FALSE"){
    newfather<-sample(sel.newfather,1)
    if(length(which(adjm1.new[newfather,]!=0))<deg1){
      adjm1.new[newfather,selection]<-1
      exit<-"TRUE"
    }
  }
  
  
  vert.deg.new<-vector()
  deg.ind.new<-vector()
  for(j in 2:max(deg1,deg2)){
    vert.deg.new<-c(vert.deg.new,which(rowSums(adjm1.new)==j)) #the adjacency matrix here needs to be not weighted
    if(length(which(rowSums(adjm1.new)==j))>0)   
      deg.ind.new<-c(deg.ind.new,rep(j,length(which(rowSums(adjm1.new)==j))))
  }
  
  #re-weight the matrix considering the new structure
  for(j in 2:max(deg1,deg2)){
    for (i in vert.deg.new[which(deg.ind.new==j)])
      adjm1.new[i,which(adjm1.new[i,]!=0)]<-c(mu[1:(j-1)],1-sum(mu[1:(j-1)]))
  }
  
  g.new<-graph.adjacency(adjm1.new, mode="directed",weighted=TRUE)
  
  vec.alpha.new<-vector(mode="numeric",length=n.leaves)#dirichlet distribution alpha vector, the variable has already been used
  
  for (i in 1:n)
    vec.alpha.new[i]<-prod(E(g.new)[get.shortest.paths(g.new,1,leaves[i],output=c("both"))$epath[[1]]]$weight) #generate alpha vector from the adjacency matrix
  vec.alpha.new<-alpha*sort(vec.alpha.new)
  
  #update the log-likelihood (and then the acceptance probability), considering the proposed tree structure
  lbeta_alpha.new<-sum(lgamma(vec.alpha.new))-lgamma(sum(vec.alpha.new))
  llikelihood.new<-0
  for (i in 1:k){
    llikelihood.new<-sum(llikelihood.new, (vec.alpha.new[which(r[[i]]==1)]-1)*log(masses_observed[[i]]), log(z[[i]])*(vec.alpha.new[which(r[[i]]==0)]-1) )
  }
  llikelihood.new<-llikelihood.new-lbeta_alpha.new
  
  feasible_nodes.new<-0
  for (i in which(rowSums(adjm1)!=0)){
    if (length(which(adjm1[i,]!=0))<deg1)
      feasible_nodes.new<-feasible_nodes.new+1
  }
  
  #Accept probability
  prob_acc<-runif(1,0,1)
  if(prob_acc<min(1,exp(llikelihood.new-llikelihood)*feasible_nodes.old/feasible_nodes.new)) 
  {
    vec.alpha<-vec.alpha.new
    lbeta_alpha<-lbeta_alpha.new
    adjm1.weighted<-adjm1.new
    vert.deg<-vert.deg.new
    deg.ind<-deg.ind.new
    adjm1<-adjm1.st #not weighted matrix
    accept.tree<-accept.tree+1
  }
  
  ##############################################################################################################
  
  memorize.vec.alpha<-rbind(memorize.vec.alpha,vec.alpha)
  memorize.mu<-rbind(memorize.mu,mu)
  memorize.z<-rbind(memorize.z,z[[1]])
  
} 

#End Gibbs sampler

#Save the outcomes

save(memorize.mu, file="musyntheticdata.RData") #values of the splitting proportions vector generated by the Gibbs sampler
save(memorize.vec.alpha, file= "iterationsalpha.RData") #values of the alpha parameter generated by the Gibbs sampler
################################Masses prediction###########################################################
exp_value<-vector(length=niterations-100)
percentili90<-matrix(nrow=niterations-100, ncol=2) #confidence interval
percentili90_average<-matrix(nrow=n.leaves,ncol=2)
expected_mass<-vector(length=n.leaves)
for (j in 1:n)
{
  for (i in 1:(niterations-100)) #evaluate the mean over the outcomes of the Gibbs sampler starting from the 100th iterations 
  {  exp_value[i]<-memorize.vec.alpha[i+100,j]/sum(memorize.vec.alpha[i+100,]) #expected values of the Dirichlet distribution for all the outcomes of the Gibbs sampler iterations
  percentili90[i,]<-qbeta(c(0.05,0.955),memorize.vec.alpha[i+100,j],(sum(memorize.vec.alpha[i+100,])-memorize.vec.alpha[i+100,j])) #using the quantile function of the beta distribution considering for each element a beta distribution with parameter (alpha_i, sum(alpha)-alpha_i) 
  }
  expected_mass[j]<-sum(exp_value)*1000/(niterations-100)
  percentili90_average[j,1]<-sum(percentili90[,1])*1000/(niterations-100) 
  percentili90_average[j,2]<-sum(percentili90[,2])*1000/(niterations-100)
}


mean<-expected_mass
lowervalue<-percentili90_average[,1] #confidence interval
uppervalue<-percentili90_average[,2] 

masses<-c(1:n)
expecteddata<-data.frame(masses,mean,lowervalue,uppervalue)
save(expecteddata, file="predictedmasses.RData") #masses with confidence interval
print(expecteddata)

























