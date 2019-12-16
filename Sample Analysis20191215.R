library(plyr)
library(survey)
library(boot)

options(stringsAsFactors = FALSE)

setwd("~/Box Sync/AuditRC/Revision/SIM submission/SIM Revision/simulated data")

#### Read in simulated data
data<-read.csv("simulated_data.csv")
head(data)
nrow(data)  ### number in cohort N=25000
table(data$Sr)  ### number in biomarker subset n=500
table(data$V)  ### number with repeat biomarkers nr=100

### Assume data not in biomarker subset is missing
data$X_b<-ifelse(data$Sr==1,data$X_b,NA)
data$Y_b<-ifelse(data$Sr==1,data$Y_b,NA)
data$X_b2<-ifelse(data$V==1,data$X_b2,NA)
data$Y_b2<-ifelse(data$V==1,data$Y_b2,NA)

##################################################################
###	 True answer: Since simulated data can consider true answer
##################################################################
Tfit<-lm(Y~X+Z,data=data)
Tans<-coef(Tfit)
summary(Tfit)
##################################################################
###	 Naive answer: using self-report only
##################################################################
#####
Nfit<-lm(Y_s~X_s+Z,data=data)
Nans<-coef(Nfit)
summary(Nfit)

###
##################################################################
###	 complete case RC using biomarker data only (correct SE using bootstrap below)
### note a 2nd method of moments esimtator for this is also at end of file. Results are similar.
##################################################################
#####

### complete case: RCBiomarker only 
	bioxRC<-lm(X_b2~X_b+Z,data=data)  ###  calibrated biomarker model
	data$xhat_bio<-predict(bioxRC,newdata=data) 
	RCfit<-lm(Y_b ~ xhat_bio  + Z, data=data)
	RCbmans<-coef(RCfit)
	summary(RCfit)

##################################################################
##### proposed method (correct bootstrap SE calculated below)
##################################################################

# regress biomarker Xbm on selfreport Xs and precise covariate Z
xRC<-lm(X_b~ X_s+Z, data=data)
#summary(xRC)
data$xhat<-predict(xRC,newdata=data)

### get regression calibration yhat
data$ydiff<-data$Y_s - data$Y_b  ### Y self-report - Y biomarker 
cstarRC<-lm(ydiff~ X_s + Z,data=data)
data$cstar<-predict(cstarRC,newdata=data)
data$yhat<- data$Y_s- data$cstar

summary(lm(yhat ~ xhat + Z, data=data))
pans<-coef(summary(lm(yhat ~ xhat + Z, data=data)))

##################################################################
##### Raking: more efficient estimator that uses proposed 
##### and biomarker data
##################################################################

### survey calibration
	fitpop<-lm(yhat~xhat+Z,data=data)  ### get influence functions of proposed method
	myinffun<-dfbeta(fitpop)
	data$infint<-myinffun[,1]
	data$infxhat<-myinffun[,2]
	data$infz<-myinffun[,3]

	myclus<-twophase(id=list(~1,~1),strata=list(NULL,NULL),subset=~I(Sr==1),data=data,method="simple")
	ifcal<-calibrate(myclus,formula=~infint+infxhat+infz,phase=2,calfun="raking") 
 	m2<-svyglm(Y_b~xhat_bio + Z ,design=ifcal)
	pans_cal<-coef(m2)
	summary(m2)
	
####### Bootstrap esitmator for proposed and biomarker RC. 

### Bootstrap
NBOOT<-500
set.seed(1223)
bans<-NULL
bansCal<-NULL
bRCbmans<-NULL
bioindex<-(1:nrow(data))[data$Sr==1]
nonbioindex<-(1:nrow(data))[data$Sr==0]

for(b in 1:NBOOT){
	bindex<-c(sample(bioindex,length(bioindex),replace=T), sample(nonbioindex,length(nonbioindex),replace=T))
	bdata<-data[bindex,]

	###RC X
	xRC<-lm(X_b~ X_s+Z, data=bdata)
	bdata$xhat<-predict(xRC,newdata=bdata)

	###RC Y
	bdata$ydiff<-bdata$Y_s - bdata$Y_b
	cstarRC<-lm(ydiff~ X_s + Z,data=bdata)
	bdata$cstar<-predict(cstarRC,newdata=bdata)
	bdata$yhat<- bdata$Y_s- bdata$cstar

	###proposed method
	bans<-rbind(bans,coef(lm(yhat ~ xhat + Z, data=bdata)))
	
	####calibrate proposed method with raking
	xRC<-lm(X_b2~X_b+Z,data=bdata)  ###  calibrated biomarker model
	bdata$xhat_bio<-predict(xRC,newdata=bdata) 

### survey calibration
	fitpop<-lm(yhat~xhat+Z,data=bdata)  ### get influence functions of proposed method
	myinffun<-dfbeta(fitpop)
	bdata$infint<-myinffun[,1]
	bdata$infxhat<-myinffun[,2]
	bdata$infz<-myinffun[,3]

	myclus<-twophase(id=list(~1,~1),strata=list(NULL,NULL),subset=~I(Sr==1),data=bdata,method="simple")
	ifcal<-calibrate(myclus,formula=~infint+infxhat+infz,phase=2,calfun="raking") 
 	m2<-svyglm(Y_b~xhat_bio  + Z,design=ifcal)
	pans_cal<-coef(m2)
	bansCal<-rbind(bansCal,pans_cal)
	
	
### complete case: RCBiomarker only
	RCfit<-lm(Y_b ~ xhat_bio  + Z, data=bdata)
	bRCbmans<-rbind(bRCbmans,coef(RCfit))
	
print(b)
	
}

### proposed
best<-apply(bans,2,mean)
bsd<-apply(bans,2,sd)
L<-best +qnorm(.025)*bsd
r<-best + qnorm(.975)*bsd
cbind(best,L,r)


###### calibrated proposed
best<-apply(bansCal,2,mean)
bsd<-apply(bansCal,2,sd)
L<-best +qnorm(.025)*bsd
r<-best + qnorm(.975)*bsd
cbind(best,L,r)


####
### biomarker rc (complete case)
best<-apply(bRCbmans,2,mean)
bsd<-apply(bRCbmans,2,sd)
L<-best +qnorm(.025)*bsd
r<-best + qnorm(.975)*bsd
cbind(best,L,r)


###### function for biomarker data only regression calibration

estimation2_b<-function(data, numiter){
  data1<-data[data$Sr==1, ]
  M<-dim(data1)[1] 
  
  ## regression calibration using repeated biomarker ##
   
  Xstar<-c(data1$X_b, data1$X_b2[data1$V==1])
  Ystar<-c(data1$Y_b, data1$Y_b2[data1$V==1])
  id<-c(1:M, (1:M)[data1$V==1])
  
  X_ave<-(data1$X_b+data1$X_b2*data1$V)/(1+data1$V)
  Y_ave<-(data1$Y_b+data1$Y_b2*data1$V)/(1+data1$V)
  
  mu_x<-mean(X_ave,na.rm=T)
  mu_y<-mean(Y_ave,na.rm=T)
  
  sigma_xstar<-sum((Xstar-mu_x)^2)/(length(Xstar)-1)
  cov_xy<-sum((Xstar-mu_x)*(Ystar-mu_y))/(length(Xstar)-1)
  sigma_t<-sum(data1$V*(data1$X_b-data1$X_b2)^2,na.rm=T)/sum(data1$V)/2
  cov_t<-sum(data1$V*(data1$X_b-data1$X_b2)*(data1$Y_b-data1$Y_b2),na.rm=T)/sum(data1$V)/2
  
  Zstar<-c(data1$Z, data1$Z[data1$V==1])
  mu_z<-mean(data$Z)
  sigma_xz<-sum((Xstar-mu_x)*(Zstar-mu_z))/(length(Xstar)-1)
  cov_zy<-sum((Ystar-mu_y)*(Zstar-mu_z))/(length(Xstar)-1)
  sigma_zhat<-sum((data1$Z-mu_z)^2)/(M-1)
  
  ## moment based estimator of beta
  b1<-solve(matrix(c(sigma_xstar-sigma_t, sigma_xz, sigma_xz, sigma_zhat), 2, 2))%*%c(cov_xy - cov_t, cov_zy)
  
  I12<-c(sigma_xstar - sigma_t, sigma_xstar - sigma_t, sigma_xz)
  I22<-matrix(c(sigma_xstar, sigma_xstar - sigma_t, sigma_xz, sigma_xstar - sigma_t, sigma_xstar, sigma_xz, sigma_xz, sigma_xz, sigma_zhat), 3,3)
  
  I12_2<-c(sigma_xstar - sigma_t, sigma_xz)
  I22_2<-matrix(c(sigma_xstar, sigma_xz, sigma_xz, sigma_zhat), 2,2)
  
  Xhat<-Cstar1<-Cstar2<-numeric(M)
  for (k in 1: M){
    
    if (data1$V[k]==1){
      Xhat[k] = mu_x + I12%*%solve(I22)%*%(c(data1$X_b[k], data1$X_b2[k], data1$Z[k]) - c(mu_x, mu_x, mu_z))
      
    } else if(data1$V[k]==0){
      Xhat[k] = mu_x + I12_2%*%solve(I22_2)%*%(c(data1$X_b[k], data1$Z[k]) - c(mu_x, mu_z))
    }
  }   
  
  xhat<-c(Xhat, Xhat[data1$V==1])
  m_regcali_bm<-lm(Ystar~xhat+Zstar)
  
  ## bootstrap standard error for regcali estimate ##
  
  regcali<-function(data1, index){
    data1<-data1[index, ]
    M<-dim(data1)[1]
    
    Xstar<-c(data1$X_b, data1$X_b2[data1$V==1])
    Ystar<-c(data1$Y_b, data1$Y_b2[data1$V==1])
    id<-c(1:M, (1:M)[data1$V==1])
    
    X_ave<-(data1$X_b+data1$X_b2*data1$V)/(1+data1$V)
    Y_ave<-(data1$Y_b+data1$Y_b2*data1$V)/(1+data1$V)
    
    mu_x<-mean(X_ave,na.rm=T)
    mu_y<-mean(Y_ave,na.rm=T)
    
    sigma_xstar<-sum((Xstar-mu_x)^2)/(length(Xstar)-1)
    cov_xy<-sum((Xstar-mu_x)*(Ystar-mu_y))/(length(Xstar)-1)
    sigma_t<-sum(data1$V*(data1$X_b-data1$X_b2)^2,na.rm=T)/sum(data1$V)/2
    cov_t<-sum(data1$V*(data1$X_b-data1$X_b2)*(data1$Y_b-data1$Y_b2),na.rm=T)/sum(data1$V)/2
    
    Zstar<-c(data1$Z, data1$Z[data1$V==1])
    mu_z<-mean(data$Z)
    sigma_xz<-sum((Xstar-mu_x)*(Zstar-mu_z))/(length(Xstar)-1)
    cov_zy<-sum((Ystar-mu_y)*(Zstar-mu_z))/(length(Xstar)-1)
    sigma_zhat<-sum((data1$Z-mu_z)^2)/(M-1)
    
    ## moment based estimator of beta
    b1<-solve(matrix(c(sigma_xstar-sigma_t, sigma_xz, sigma_xz, sigma_zhat), 2, 2))%*%c(cov_xy - cov_t, cov_zy)
    
    I12<-c(sigma_xstar - sigma_t, sigma_xstar - sigma_t, sigma_xz)
    I22<-matrix(c(sigma_xstar, sigma_xstar - sigma_t, sigma_xz, sigma_xstar - sigma_t, sigma_xstar, sigma_xz, sigma_xz, sigma_xz, sigma_zhat), 3,3)
    
    I12_2<-c(sigma_xstar - sigma_t, sigma_xz)
    I22_2<-matrix(c(sigma_xstar, sigma_xz, sigma_xz, sigma_zhat), 2,2)
    
    Xhat<-Cstar1<-Cstar2<-numeric(M)
    for (k in 1: M){
      
      if (data1$V[k]==1){
        Xhat[k] = mu_x + I12%*%solve(I22)%*%(c(data1$X_b[k], data1$X_b2[k], data1$Z[k]) - c(mu_x, mu_x, mu_z))
        
      } else if(data1$V[k]==0){
        Xhat[k] = mu_x + I12_2%*%solve(I22_2)%*%(c(data1$X_b[k], data1$Z[k]) - c(mu_x, mu_z))
      }
    }   
    
    xhat<-c(Xhat, Xhat[data1$V==1])
    m_regcali_bm<-lm(Ystar~xhat+Zstar)
    return(m_regcali_bm$coef[2:3])
  }
  
  output<-boot(data1, regcali, strata=data1$V, R=numiter)
  m_regcali_bm_SE<-apply(output$t, 2, sd)
  conf_x<-boot.ci(output, index=1, type = c("perc"))$perc[4:5]
  conf_z<-boot.ci(output, index=2, type = c("perc"))$perc[4:5]
  
  return(list(m_regcali_bm$coef[2:3], m_regcali_bm_SE, conf_x, conf_z))
}

  output2<-estimation2_b(data, NBOOT)


