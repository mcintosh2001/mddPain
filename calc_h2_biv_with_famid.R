traits<-c("dep_status","CPGquant","CPG2cases")

sink("FamID_Sibs_rg_rp.txt",append=TRUE)
for(i in 1:(length(traits) - 1)){
  for(j in (i+1):length(traits)){
    
    formula <- paste("fixed=cbind(",traits[i],",",traits[j],") ~trait + trait:sex + trait:age + trait:I(age^2)")
    
    print(formula)
    
    asr.model<- asreml(fixed=as.formula(formula)
                       ,random= ~us(trait,init=c(1,0.1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(famid)
                       ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                       ,ginverse=list(id=ainv)
                       ,data=dataframe
                       ,family=asreml.gaussian(link="identity",dispersion=NA)
                       ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
    ##Assuming no covariance for the null model
    
    asr.model.null<- asreml(fixed=as.formula(formula)
                            ,random= ~diag(trait,init=c(1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(famid)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
  
    print(summary(asr.model)$varcomp)
    
    print(asr.model$converge)
    print(asr.model.null$converge)
    
    #Additive covariance
    covA<-summary(asr.model)$varcomp[2,2]
    
    #Additive variance of trait[i] 
    Vg1A<-summary(asr.model)$varcomp[1,2]
    
    #Additive variance of trait[j]
    Vg2A<-summary(asr.model)$varcomp[3,2]
    
    #Genetic correlation between trait[i] and trait[j]
    gCORR<-covA/(sqrt(Vg1A*Vg2A))
#     
#     print(paste("Genetic correlation of ", traits[i]," and ",traits[j],"is = ",gCORR))
#     
#     #calculating significance between model and null model  
#     diffloglik <- 2*(asr.model$loglik-asr.model.null$loglik)
#     signif <- 1-pchisq(diffloglik,1)
#     
#     print(paste0("Chi2 = ",diffloglik,", p = ",signif))
#     
#     print(pin(asr.model, ~ V2/sqrt(V1*V3)))
#     
#     ##Phenotypic correlation
#     
#     ##Phenotypic covariance
#     #####Additive covariance + sib covariance + spouse covariance + residual covariance
#     covP<-summary(asr.model)$varcomp[2,2]+summary(asr.model)$varcomp[5,2]+summary(asr.model)$varcomp[8,2]+summary(asr.model)$varcomp[12,2]
#     
#     ##Phenotypic variance of trait[i]
#     #####Additive variance + sib variance + spouse variance + residual variance for trait[i]
#     Vp1<-summary(asr.model)$varcomp[1,2]+summary(asr.model)$varcomp[4,2]+summary(asr.model)$varcomp[7,2]+summary(asr.model)$varcomp[11,2]
#     
#     
#     ##Phenotypic variance of trait[j] 
#     #####Additive variance + sib variance + spouse variance + residual variance for trait[j]
#     Vp2<-summary(asr.model)$varcomp[3,2]+summary(asr.model)$varcomp[6,2]+summary(asr.model)$varcomp[9,2]+summary(asr.model)$varcomp[13,2]
#     
#     #Phenotypic correlation between trait[i] and trait[j]
#     pCORR<-covP/(sqrt(Vp1*Vp2))
#     
#     
#     print(paste("Phenotypic correlation of ", traits[i]," and ",traits[j],"is = ",pCORR))
#     
#     print(pin(asr.model, ~ (V2+V5+V8+V12)/sqrt((V1+V4+V7+V11)*(V3+V6+V9+V13))))
#     
#     ##Calculating significance of phenotypic correlation
#     ## Standard error of phenotypic covariance component
#     SE_P<-CovP<-summary(asr.model)$varcomp[2,3]+summary(asr.model)$varcomp[5,3]+summary(asr.model)$varcomp[8,3]+summary(asr.model)$varcomp[12,3]
#     
#     ## Z Ratio of phenotypic covariance component
#     z <- covP/SE_P
#     
#     ## Corresponding p value for Z ratio
#     p <- 2*pnorm(-abs(z)) 
#     
#     print(paste(" p = ",p))
    

    ##hij=sqrt(hi)*sqrt(hj)*gCORR
    
  }
}
sink()
