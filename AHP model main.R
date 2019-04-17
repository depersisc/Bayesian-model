###############################################################################################
# NPI activity “Stochastic Modelling of Atmospheric Re-entry Highly Energetic Break-up Events”
# Belief-network model for failure prediction
#
#Goal: Generate a plot of the probabilities for the failure to occur in various events, given pairwise comparisons from experts and observations. Example of application.
################################################################################################
rm(list=ls())
library(psych)
########################################################################################################
#Input variables
#Collect pairwise comparisons

#Set A 
#A1 Burst of a pressure vessel
#A2 Chemical reaction propellant+air
#A3 Chemical reaction between hypergolic propellants
#A4 Burst of a battery cell
#a_ij=Ai vs Aj

a12<-5
a13<-7
a14<-7
a23<-5
a24<-1
a34<-1

#SetB
#B1 Sudden release of propellant
#B2 Slow release of propellant
b12<-3

#SetC
#C1 Valve leakage
#C2 Tank destruction
#C3 Pipe rupture

c12<-1
c13<-9
c23<-9

#SetD
#D1 Chemical reactions
#D2 Overpressure
#D3 Short circuit
#D4 Corrosion
#D5 Overcharge
#D6 Overdischarge

d12<-1
d13<-1
d14<-1
d15<-1
d16<-1
d23<-1
d24<-1
d25<-1
d26<-1
d34<-1
d35<-1
d36<-1
d45<-1
d46<-1
d56<-1

number_reentries<-3 #number of observations. Default:3
number_events<-13 #number of basic events
#Observations vector for nominal case. k_vector[i]=1 if the ith event occurs
k_vector<-vector( mode="numeric",length=number_events)
k_vector[5:8]<-1 #to evaluate the likelihood of nominal case

##########################################################################################
#Collect tendency-to-occur quotients for new re-entries
#
q<-matrix(nrow=number_events,ncol=(number_reentries+1))

q[,1]<-c(1,rep(1,number_events-1))
q[,2]<-c(1,1,3,3,3,rep(1,7),1)
q[,3]<-c(1,1,1/5,1/5,1/5,rep(1,8))
q[,4]<-c(1,1,1,1,1,1,1,1,7,1,7,1,1)
###End of user input################################################
A<-matrix(nrow=4,ncol=4)  

A[1,]<-c(1,a12,a13,a14)
A[2,]<-c(1/a12,1,a23,a24)
A[3,]<-c(1/a13,1/a23,1,a34)
A[4,]<-c(1/a14,1/a24,1/a34,1)

#consistency check
CI_A<-(eigen(A)$values[1]-dim(A)[1])/(dim(A)[1]-1)

B<-matrix(nrow=2,ncol=2)  
B[1,]<-c(1,b12)
B[2,]<-c(1/b12,1)

#consistency check
CI_B<-(max(eigen(B)$values)-dim(B)[1])/(dim(B)[1]-1)

C<-matrix(nrow=3,ncol=3)  

C[1,]<-c(1,c12,c13)
C[2,]<-c(1/c12,1,c23)
C[3,]<-c(1/c13,1/c23,1)

#consistency check
CI_C<-(max(eigen(C)$values)-dim(C)[1])/(dim(C)[1]-1)

D<-matrix(nrow=8,ncol=8)  
D[1,]<-c(1,d12,d13,d14,d15,d16,rep(1,2))
D[2,]<-c(1/d12,1,d23,d24,d25,d26, rep(1,2))
D[3,]<-c(1/d13,1/d23,1,d34,d35,d36,rep(1,2))
D[4,]<-c(1/d14,1/d24,1/d34,1,d45,d46,rep(1,2))
D[5,]<-c(1/d15,1/d25,1/d35,1/d45,1,d56,rep(1,2))
D[6,]<-c(1/d16,1/d26,1/d36,1/d46,1/d56,1,rep(1,2))
D[7,]<-rep(1,8)
D[8,]<-rep(1,8)

#consistency check
CI_D<-(max(eigen(D)$values)-dim(D)[1])/(dim(D)[1]-1)
########################################################################################################
#Tendency-to-occur weights######################################################
weights_a<-geometric.mean(A)/sum(geometric.mean(A)) #weights sum to 1
weights_b<-geometric.mean(B)/sum(geometric.mean(B))
weights_c<-geometric.mean(C)/sum(geometric.mean(C))
weights_d<-geometric.mean(D)/sum(geometric.mean(D))

weights_base<-c(weights_a[1],weights_a[3],weights_c,weights_d)
CR<-c(0,0,0.58,0.9,1.12,1.24,1.32,1.41) #taken from literature
CR_A<-CI_A/CR[dim(A)[1]]
CR_C<-CI_C/CR[dim(C)[1]]
CR_D<-CI_D/CR[dim(D)[1]]

##########Prior quantification under nominal conditions#####################################################################################
#Assume the expert is more confident about knowing the probability of the event C2 (tank destruction)
#Define a range for theta_j tank destruction set C 
theta_min_tank<-0.01
theta_max_tank<-0.04
#Fit a beta distribution to this range by matching the mean +-2 standard deviation
mean_tankdestruction<-mean(c(theta_min_tank,theta_max_tank))
dev_stand_tankdestruction<-(theta_max_tank-mean_tankdestruction)/2
#Find shape1 and shape2 for the beta distribution, given mean and variance of the sample. 

fun<-function(mu=mu,sigmasquare=sigmasquare){
  f<-numeric(2)
  f[1]<-(mu-mu^2-sigmasquare)*mu/(sigmasquare)
  f[2]<-f[1]/mu-f[1]
  f
}
shape_tankdestruction<-fun(mu=mean_tankdestruction, sigmasquare=dev_stand_tankdestruction^2)

ranges<-matrix(nrow=number_events,ncol=2)
ranges<-t(t(weights_base))%*%(c(theta_min_tank,theta_max_tank)/weights_c[2]) #evaluation of the ranges for all the base events 
means<-rowMeans(ranges)
dev_stand<-(ranges[,1]-means)/2
shapes_nominalprior<-matrix(nrow=number_events,ncol=2)
shapes<-matrix(nrow=number_events,ncol=2)
for (i in 1:number_events)
  shapes_nominalprior[i,]<-fun(mu=means[i],sigmasquare=dev_stand[i]^2) #parameters of beta prior distribution for all the base events under nominal conditions

#End of Prior elicitation step under nominal conditions 
########################################################################################################

#Compute posterior in the nominal case
#Parameters of the beta posterior in the nominal case
shapes[,1]<-shapes_nominalprior[,1]+k_vector
shapes[,2]<-shapes_nominalprior[,2]-k_vector+1 

#################################################################
#Probability of explosion 
nsimul<-100000 #number of samples of the probabilities for each event to occur and number of simulation of the importance sampling algorithm

Prob_expl<-matrix(nrow=nsimul, ncol=(number_reentries+1))
PE_part<-vector(length=nsimul)

priorsimul<-array(0, dim=c(nsimul,number_events, (number_reentries+1)))
for (ind in 1:number_events){
  priorsimul[,ind,1]<-rbeta(nsimul,shapes[ind,1],shapes[ind,2]) #sampling probabilities from nominal case
}

########################################################################################################
#Use tendency-to-occur quotients and current posterior to make prediction for P(explosion) for new re-entry
#Update the likelihood with AHP scores

###########Likelihood updates#################### Importance sampling implementation##################
priorsim<-matrix(nrow=nsimul,ncol=number_events)
for (ind in 1:number_events){
  priorsim[,ind]<-rbeta(nsimul,shapes_nominalprior[ind,1],shapes_nominalprior[ind,2]) #sample from nominal prior 
}

weight_imp<-array(0, dim=c(nsimul, number_reentries)) 
weight_imp_norm<-array(0, dim=c(nsimul, number_reentries))

#this loop needs to be generalized
for (i in 1:nsimul){
  weight_imp[i,1]<-prod(priorsim[i,5:8])*(1-prod(1-priorsim[i,3:5]^(1/q[3:5,2]))) #1st re-entry 
  weight_imp[i,2]<-weight_imp[i,1]*(1-prod(1-priorsim[i,3:5]^(1/q[3:5,3]))) #2nd re-entry
  weight_imp[i,3]<-weight_imp[i,2]*prod(1-priorsim[i,3:5]^(1/q[3:5,4])) #3rd re-entry
}

for (i in 1:number_reentries)
  weight_imp_norm[,i]<-weight_imp[,i]/sum(weight_imp[,i])

for(k in 2:(number_reentries+1)){
  for (f in 1:number_events)
  {priorsimul[,f,k]<-sample(priorsim[,f],nsimul,replace=TRUE,prob=weight_imp_norm[,(k-1)])}
}


for (k in 1:(number_reentries+1)){
  for (i in 1:nsimul){
    PE_part[i]<-(1-priorsimul[i,1,k]^(1/q[1,k]))
    for(j in 2:number_events){
      PE_part[i]<-(1-priorsimul[i,j,k]^(1/q[j,k]))*PE_part[i]}
    Prob_expl[i,k]<-1-PE_part[i]
  }
  PE_part<-rep(0,nsimul)
  
}

######################################Evaluation of the marginal likelihoods#############################
llikelihood<-function(comb)
{
  return(log(1-prod(1-comb^(1/q[3:5,2])))+log(1-prod(1-comb^(1/q[3:5,3])))+sum(log(1-comb^(1/q[3:5,4]))))
}

n_thetac<-50  #number of thetac values
st_thetac<-0.001
end_thetac<-0.03
thetac_values<-seq(st_thetac,end_thetac,length=n_thetac)


deltathetac<-thetac_values[2]-thetac_values[1]
combinations<-cbind(thetac_values,sort(rep(thetac_values,n_thetac)),sort(rep(thetac_values,n_thetac^2))) 
#loglikelihood of B2
loglikelihood<-apply(combinations,1,llikelihood)
norm_const<-sum(exp(loglikelihood))*(deltathetac^3)
log_likelihood_norm<-loglikelihood-log(norm_const)

#marginal distributions
d<-split(log_likelihood_norm,ceiling(seq_along(log_likelihood_norm)/n_thetac^2))
f<-split(log_likelihood_norm,ceiling(seq_along(log_likelihood_norm)/n_thetac))
md<-do.call(rbind,d)

matsplitter<-function(M, r, c) {
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
} 

rd<-as.vector(matsplitter(md,n_thetac,n_thetac))
r<-split(rd,ceiling(seq_along(rd)/n_thetac^2))


mf<-do.call(rbind,f)
e<-split(mf,col(mf))

#marginal distribution C1

marg1<-function(nr){
  return(sum(exp(e[[nr]])))
}
marginal_thetac1<-apply(cbind(1:n_thetac),1,marg1)
marginal_thetac1<-marginal_thetac1*(deltathetac)^2 

#marginal distribution C2

marg2<-function(nr){
  return(sum(exp(r[[nr]])))
}
marginal_thetac2<-apply(cbind(1:n_thetac),1,marg2)
marginal_thetac2<-marginal_thetac2*(deltathetac)^2

#marginal distribution C3

marg3<-function(nr){
  return(log(sum(exp(d[[nr]]))))
}
log_marginal_thetac3<-apply(cbind(1:length(thetac_values)),1,marg3)
marginal_thetac3<-exp(log_marginal_thetac3)
marginal_thetac3<-marginal_thetac3*(deltathetac)^2
########################End of evaluation of marginal likelihoods###############################

#Generate plot for thesis#####################################################################
pdf(file="Failureprobabilities.pdf",height=50,width=50)
par(mfrow=c(1,1),cex="5")
hist(Prob_expl[,1],freq=FALSE,col=rgb(0,0,1,0.3),xlim=c(0.1,1),ylim=c(0,50),main="P(TE=1)", xlab="",ylab="",lwd="3",cex.axis="2",breaks=50)#old posterior
hist(Prob_expl[,2],freq=FALSE,col=rgb(1,0,1,1),breaks=50, add=T)#new posterior
hist(Prob_expl[,3],freq=FALSE,col=rgb(0,1,1,0.5),breaks=50, add=T)#new posterior
hist(Prob_expl[,4],freq=FALSE,col=rgb(0,0,1,1),breaks=50, add=T)#new posterior
legend("topleft", c("Nominal", "1st observation", "2nd observation", "3rd observation"), col=c(rgb(0,0,1,0.3), rgb(1,0,1,1), rgb(0,1,1,0.5),rgb(0,0,1,1)), lwd=7)
dev.off()


