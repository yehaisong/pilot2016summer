pilot.data<-read.csv("~/pilot050220156/pilot data_lp and cl.csv")
#print(pilot.data)
#calculate learning performance variance group va
pilot.data$GVA_LP_V<-pilot.data$GVA_LP_sd^2
#calculate learning performance variance group vv
pilot.data$GVV_LP_V<-pilot.data$GVV_LP_sd^2
#calculate learning performance sample covariance
pilot.data$LP_V<-(pilot.data$GVA_LP_sd^2*(pilot.data$GVA_n-1)+pilot.data$GVV_LP_sd^2*(pilot.data$GVV_n-1))/(pilot.data$GVA_n+pilot.data$GVV_n-2)
#calculate learning performance variance group va
pilot.data$GVA_CL_V<-pilot.data$GVA_CL_sd^2
#calculate learning performance varnce group vv
pilot.data$GVV_CL_V<-pilot.data$GVV_CL_sd^2
#calculate cognitive load sample covariance
pilot.data$CL_V<-(pilot.data$GVA_CL_sd^2*(pilot.data$GVA_n-1)+pilot.data$GVV_CL_sd^2*(pilot.data$GVV_n-1))/(pilot.data$GVA_n+pilot.data$GVV_n-2)

#calculate effect sizes
ME<-matrix(c((pilot.data$GVA_LP_m-pilot.data$GVV_LP_m)/sqrt(pilot.data$LP_V),(pilot.data$GVA_CL_m-pilot.data$GVV_CL_m)/sqrt(pilot.data$CL_V)),ncol=2)
#assume r=0
r<-0
MEVCOV<-matrix(c(((1/pilot.data$GVA_n)+(1/pilot.data$GVV_n)+(ME[,1]^2/(2*(pilot.data$GVA_n+pilot.data$GVV_n))))
                 ,(1/pilot.data$GVA_n+1/pilot.data$GVV_n)*r+(((pilot.data$GVA_CL_m-pilot.data$GVV_CL_m)/sqrt(pilot.data$CL_V)*(pilot.data$GVA_LP_m-pilot.data$GVV_LP_m)/sqrt(pilot.data$LP_V))/2*(pilot.data$GVA_n+pilot.data$GVV_n))*(r^2)
                 ,(1/pilot.data$GVA_n+1/pilot.data$GVV_n+((ME[,2])^2/(2*(pilot.data$GVA_n+pilot.data$GVV_n))))
                 ),ncol=3)

result<-meta(ME,MEVCOV,model.name="Random effects model")
lower<-'5.906044
1.810071 0.777015'
Cov1<-getCov(lower, diag=TRUE)
cov2cor(Cov1)

print(ME)
print(MEVCOV)

summary(result)
print(pilot.data)




 