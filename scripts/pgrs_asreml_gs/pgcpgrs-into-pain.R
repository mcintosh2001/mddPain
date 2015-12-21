setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")
source("scripts/include/run_analysis.R")

traits<-c("CPG2cases","CPGquant")
scores <-c("bdldscore","sczldscore","mddldscore")

sink(paste0(locations$results,"/pgrs-gs/pgcpgrs_into_pain",today,".txt"),append=FALSE)
cat("PGC polygenic prediction of pain in GS\n")
sink()

for(i in 1:(length(traits))){
  for(j in 1:length(scores)){
  
    
  asr.formula <- paste("fixed=",traits[i]," ~sex + age + I(age^2) +C1+C2+C3+C4 +", scores[j])
  
  asr.model<- asreml(fixed=as.formula(asr.formula)
                     ,random= ~ped(id, var=T, init=1) 
                     ,ginverse=list(id=ainv)
                     ,data=dataframe
                     ,maxiter=15000,trace=FALSE,na.method.X="omit",workspace=8e+08)
  
  sink(paste0(locations$results,"/pgrs-gs/pgcpgrs_into_pain",today,".txt"),append=TRUE)
  
  cat("\n")
  print(asr.model$fixed.formula)
  print(summary(asr.model, all=T)$coef.fixed)
  
  beta_ldpred <- summary(asr.model, all=T)$coef.fixed[1,1]
  print(paste0("beta for ", scores[j], " = ", beta_ldpred))
  
  print(wald.asreml(asr.model, ssType="conditional", denDF="numeric"))

#calculate predicted phenotype value 
dataframe$temp <- dataframe[,scores[j]]*beta_ldpred

var_explained <- var(dataframe$temp, na.rm=TRUE)
var_total <- var(dataframe[,traits[i]], na.rm=TRUE)
r2<- var_explained/var_total
print(paste0("variance explained by ",scores[j], " in ",traits[i]," = ",r2))
r2 <- var_explained <- var_total <- NULL 
  cat("\n")
  sink()



}
}
