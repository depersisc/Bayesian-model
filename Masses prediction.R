###############################################################################################
#Author: Cristina De Persis
# NPI activity “Stochastic Modelling of Atmospheric Re-entry Highly Energetic Break-up Events”
# 
#masses prediction (see section 6.8.2 of the thesis)
################################################################################################
exp_value<-vector(length=niterations-100)
percentili90<-matrix(nrow=niterations-100, ncol=2) #confidence interval
percentili90_average<-matrix(nrow=n.leaves,ncol=2)
expected_mass<-vector(length=n.leaves)
for (j in 1:n.leaves)
{
  for (i in 1:(niterations-100)) #evaluate the mean over the outcomes of the Gibbs sampler starting from the 100th iterations 
  {  exp_value[i]<-memorize.vec.alpha[i+100,j]/sum(memorize.vec.alpha[i+100,]) #expected values of the Dirichlet distribution for all the outcomes of the Gibbs sampler iterations
  percentili90[i,]<-qbeta(c(0.05,0.955),memorize.vec.alpha[i+100,j],(sum(memorize.vec.alpha[i+100,])-memorize.vec.alpha[i+100,j]))
  }
  expected_mass[j]<-sum(exp_value)/(niterations-100)
  percentili90_average[j,1]<-sum(percentili90[,1])/(niterations-100)
  percentili90_average[j,2]<-sum(percentili90[,2])/(niterations-100)
}

mean<-expected_mass
lowervalue<-percentili90_average[,1] #confidence interval
uppervalue<-percentili90_average[,2] 

masses<-c(1:n.leaves)
expecteddata<-data.frame(masses,mean,lowervalue,uppervalue)
print(expecteddata)