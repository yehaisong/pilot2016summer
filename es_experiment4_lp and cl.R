library(metaSEM)
library("lavaan")
##Multiple-endpoint studies

#control group
lower <- '10.24
0, 1.9321'

lower
##Convert a lower triangle data into a covariance matrix
Cov1<-getCov(lower, diag=TRUE, name=c("LP","CL"))

##treatment group
lower <- '12.6736
0, 2.25'
##Convert a lower triangle data into a covraiance matrix
Cov2<-getCov(lower, diag=TRUE,name=c("LP","CL"))

#Convert covariance matrices into a list
Cov<-list(Cov1,Cov2)

#Means for the two groups
Mean<-list(c(5.50,4.05),c(5.55,3.43))

#Sample sizes for the gorups
N <-c(20,16)

#Assuming homogeneity of covariance matrices by 
#using the same lables: "sd1", "sd2", and "r"
model4<-'eta1=~c("sd1","sd1")*LP
eta2=~c("sd2","sd2")*CL
eta1~~c("r","r")*eta2
LP~c("m1_1","m1_2")*1
CL~c("m2_1","m2_2")*1
LP~~0*LP
CL~~0*CL
# Multiple endpoint effect size 1
MELP:=(m1_2 - m1_1)/sd1
#Multiple endpoint effect size 2
MECL:=(m2_2 - m2_1)/sd2'
fit4<-sem(model4, sample.cov = Cov, sample.mean = Mean, 
          sample.nobs = N,std.lv = TRUE, 
          sample.cov.rescale = FALSE)
##Obtain the free parameters in the model
(x<-fit4@Fit@x)
##Obtain the sampling covariance matrix of the parameter estimates
(VCOV<-vcov(fit4))
##Compute the multivariate effect sizes
(ME<-fit4@Model@def.function(x))
##Compute the jacobian for 'defined parameters'
JAC<-lavaan:::lavJacobianD(func=fit4@Model@def.function, x=x)
##compute the sampling covariance matrix using delta method
ME.VCOV<-JAC%*% VCOV %*% t(JAC)
##Add the variable names for ease of reference
dimnames(ME.VCOV) <- list(names(ME),names(ME))
ME.VCOV