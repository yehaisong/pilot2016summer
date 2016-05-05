RTCL <- read.csv("C:/Users/hsye/Desktop/pilot050220156/pilot data_r and t and cl.csv")
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

#assume r=0
r<-0
ME.VCOV<-matrix(c(((1/RTCL$GVA_n)+(1/RTCL$GVV_n)+(ME[,1]^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#RR
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_T_m-RTCL$GVV_T_m)/sqrt(RTCL$T_V)*(RTCL$GVA_R_m-RTCL$GVV_R_m)/sqrt(RTCL$R_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#RT
                 ,((1/RTCL$GVA_n)+(1/RTCL$GVV_n)+(ME[,2]^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#TT
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_CL_m-RTCL$GVV_CL_m)/sqrt(RTCL$CL_V)*(RTCL$GVA_R_m-RTCL$GVV_R_m)/sqrt(RTCL$R_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#RCL
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n)*r+(((RTCL$GVA_CL_m-RTCL$GVV_CL_m)/sqrt(RTCL$CL_V)*(RTCL$GVA_T_m-RTCL$GVV_T_m)/sqrt(RTCL$T_V))/2*(RTCL$GVA_n+RTCL$GVV_n))*(r^2)#TCL
                 ,(1/RTCL$GVA_n+1/RTCL$GVV_n+((ME[,3])^2/(2*(RTCL$GVA_n+RTCL$GVV_n))))#CLCL
                 ),ncol=6)
#ME<-ME[,c(1,2)]
#ME.VCOV<-ME.VCOV[,c(1,2,3)]
result<-meta(ME,ME.VCOV, model.name = "Random effects model")
lower<-'5.906044
1.810071 0.777015'
Cov1<-getCov(lower, diag=TRUE)
cov2cor(Cov1)

print(ME)
print(ME.VCOV)

summary(result)
print(RTCL)




 