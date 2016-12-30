library(metaSEM)
library(meta)
RTCL <- read.csv("pilot data_r and t and cl.csv")

#print(RTCL)
#calculate RECALL performance variance group va
RTCL$GVA_R_V<-RTCL$GVA_LP_sd^2
#calculate recall performance variance group vv
RTCL$GVV_R_V<-RTCL$GVV_LP_sd^2
#calculate recall performance sample covariance
RTCL$R_V<-(RTCL$GVA_R_sd^2*(RTCL$GVA_n-1)+RTCL$GVV_R_sd^2*(RTCL$GVV_n-1))/(RTCL$GVA_n+RTCL$GVV_n-2)
#calculate transfer performance variance group va
RTCL$GVA_T_V<-RTCL$GVA_LP_sd^2
#calculate transfer performance variance group vv
RTCL$GVV_T_V<-RTCL$GVV_LP_sd^2
#calculate transfer performance sample covariance
RTCL$T_V<-(RTCL$GVA_T_sd^2*(RTCL$GVA_n-1)+RTCL$GVV_T_sd^2*(RTCL$GVV_n-1))/(RTCL$GVA_n+RTCL$GVV_n-2)
#calculate cL variance group va
RTCL$GVA_CL_V<-RTCL$GVA_CL_sd^2
#calculate CL varnce group vv
RTCL$GVV_CL_V<-RTCL$GVV_CL_sd^2
#calculate CL sample covariance
RTCL$CL_V<-(RTCL$GVA_CL_sd^2*(RTCL$GVA_n-1)+RTCL$GVV_CL_sd^2*(RTCL$GVV_n-1))/(RTCL$GVA_n+RTCL$GVV_n-2)

#calculate effect sizes
ME<-matrix(c((RTCL$GVA_R_m-RTCL$GVV_R_m)/sqrt(RTCL$R_V)
             ,(RTCL$GVA_T_m-RTCL$GVV_T_m)/sqrt(RTCL$T_V)
             ,(RTCL$GVA_CL_m-RTCL$GVV_CL_m)/sqrt(RTCL$CL_V))
           ,ncol=3)
#insert effect sizes to RTCL
RTCL$ES_R<-ME[,1]
RTCL$ES_T<-ME[,2]
RTCL$ES_CL<-ME[,3]

#assume r=0
r<-0
ME.VCOV<-matrix(c(((1/RTCL$GVA_n)+(1/RTCL$GVV_n)+(ME[,1]^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#RR
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_T_m-RTCL$GVV_T_m)/sqrt(RTCL$T_V)*(RTCL$GVA_R_m-RTCL$GVV_R_m)/sqrt(RTCL$R_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#RT
                 ,((1/RTCL$GVA_n)+(1/RTCL$GVV_n)+(ME[,2]^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#TT
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_CL_m-RTCL$GVV_CL_m)/sqrt(RTCL$CL_V)*(RTCL$GVA_R_m-RTCL$GVV_R_m)/sqrt(RTCL$R_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#RCL
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_CL_m-RTCL$GVV_CL_m)/sqrt(RTCL$CL_V)*(RTCL$GVA_T_m-RTCL$GVV_T_m)/sqrt(RTCL$T_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#TCL
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n+((ME[,3])^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#CLCL
                 ),ncol=6)

#insert covariate to RTCL
RTCL$var_R<-ME.VCOV[,1]
RTCL$cov_R_T<-ME.VCOV[,2]
RTCL$var_T<-ME.VCOV[,3]
RTCL$cov_R_CL<-ME.VCOV[,4]
RTCL$cov_T_CL<-ME.VCOV[,5]
RTCL$var_CL<-ME.VCOV[,6]

#ME<-ME[,c(1,2)]
#ME.VCOV<-ME.VCOV[,c(1,2,3)]

result<-meta(y=cbind(ES_R,ES_T,ES_CL),v=cbind(var_R,cov_R_T,cov_R_CL,var_T,cov_T_CL,var_CL),data=RTCL, model.name = "Random effects model")
result1<-meta(y=cbind(ES_R,ES_T,ES_CL),v=cbind(var_R,cov_R_T,cov_R_CL,var_T,cov_T_CL,var_CL),data=RTCL, model.name = "Random effects model",RE.constraints = Diag(paste("0.1*Tau2",1:3,1:3,sep="_")))

#Extract the variance component of the random effects
(T2<-vec2symMat(coef(result1,select="random")))
#Convert the covariance matrix to a correlation matrix
cov2cor(T2)

#lower<-'5.906044
#1.810071 0.777015'
#Cov1<-getCov(lower, diag=TRUE)
#cov2cor(Cov1)

print(ME)
print(ME.VCOV)

summary(result)
summary(result1)
print(RTCL)

plot(result1)

#Add moderator variables
MV<-read.csv("Workbook1.csv")
RTCL$MV_TL<-MV$textlength 
RTCL$MV_PR<-MV$peerreview
RTCL$MV_PP<-MV$presentationpace


result2<-meta(y=cbind(ES_R,ES_T,ES_CL),v=cbind(var_R,cov_R_T,cov_R_CL,var_T,cov_T_CL,var_CL),x=cbind(MV_PP), data=RTCL, model.name = "Random effects model",RE.constraints = Diag(paste("0.1*Tau2",1:3,1:3,sep="_")))
summary(result2)

#meta package
#Add short-term retention score
temp<-read.csv("pilot data1.csv")
RTCL$GVA_SR_m<-temp$GVA_SR_m
RTCL$GVA_SR_sd<-temp$GVA_SR_sd
RTCL$GVV_SR_m<-temp$GVV_SR_m
RTCL$GVV_SR_sd<-temp$GVV_SR_sd
result3<-metacont(GVA_n,GVA_SR_m,GVA_SR_sd,GVV_n,GVV_SR_m,GVV_SR_sd,data=RTCL,sm="SMD", method.smd="Cohen")
print(result3)
result3.mr<-metareg(result3,MV_PR+MV_PP)
print(result3.mr)
bubble(result3.mr, ylim=c(-4,10),xlim=c(-1,2),regline=TRUE)
 