traits<-c("dep_status","CPGquant","CPG2cases","likert_total","eysenck_N","gfactor")


today <- format(Sys.Date(),format="%Y%B%d")

sink(paste0(locations$results,"/",today,"sib_rg_rp.txt"),append=FALSE)
for(i in 1:(length(traits) - 1)){
  for(j in (i+1):length(traits)){
    
    formula <- paste("fixed=cbind(",traits[i],",",traits[j],") ~trait + trait:sex + trait:age + trait:I(age^2)")
    
    print(formula)
    
    asr.model<- asreml(fixed=as.formula(formula)
                       ,random= ~us(trait,init=c(1,0.1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(sib)
                       ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                       ,ginverse=list(id=ainv)
                       ,data=dataframe
                       ,family=asreml.gaussian(link="identity",dispersion=NA)
                       ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
    ##Assuming no genetic covariance for the null model
    asr.model.nullg<- asreml(fixed=as.formula(formula)
                            ,random= ~diag(trait,init=c(1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(sib)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
    ##Assuming no environmental covariance for the null model
    asr.model.nulle<- asreml(fixed=as.formula(formula)
                            ,random= ~us(trait,init=c(1,0.1,1)):ped(id)+diag(trait,init=c(1,1)):id(sib)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
  
    print(summary(asr.model)$varcomp)
    
    print(paste0("Did full model converge: ",asr.model$converge))
    print(paste0("Did null-genetic model converge: ",asr.model.nullg$converge))
    print(paste0("Did null-environment model converge: ",asr.model.nulle$converge))
    
    #Additive covariance
    covA<-summary(asr.model)$varcomp[2,2]
    #Environmental covariance
    eCov<-summary(asr.model)$varcomp[5,2]
          
    #Additive variance of trait[i] 
    Vg1A<-summary(asr.model)$varcomp[1,2]
    #Environmental variance of trait[i]
    Ve1<-summary(asr.model)$varcomp[4,2]
          
    #Additive variance of trait[j]
    Vg2A<-summary(asr.model)$varcomp[3,2]
    #Environmental variance of trait[i]
    Ve2<-summary(asr.model)$varcomp[6,2]
    
    #Genetic correlation between trait[i] and trait[j]
    gCORR<-covA/(sqrt(Vg1A*Vg2A))
    eCORR<-eCov/(sqrt(Ve1*Ve2))
    
    #calculating significance between model and null-genetic model  
    diffloglikg <- 2*(asr.model$loglik-asr.model.nullg$loglik)
    signifg <- 1-pchisq(diffloglikg,1)
    print(paste0("Genetic Chi2 = ",diffloglikg,", p = ",signifg))
    
    #calculating significance between model and null-environmental model  
    diffloglike <- 2*(asr.model$loglik-asr.model.nulle$loglik)
    signife <- 1-pchisq(diffloglike,1)
    print(paste0("Environment Chi2 = ",diffloglike,", p = ",signife))
    
    print("rg and its standard error")
    print(pin(asr.model, ~ V2/sqrt(V1*V3)))
    
    print("re and its standard error")
    print(pin(asr.model, ~ V5/sqrt(V4*V6)))
    
    print("r and its standard error")
    print(pin(asr.model, ~ (V2+V5+V9)/(sqrt((V1+V4+V8)*(V3+V6+V10)))     ))
    
          
    ##Phenotypic correlation
    
    ##Phenotypic covariance
    #####Additive covariance + environmental covariance + residual covariance
    covP<-summary(asr.model)$varcomp[2,2]+summary(asr.model)$varcomp[5,2]+summary(asr.model)$varcomp[9,2]
    
    ##Phenotypic variance of trait[i]
    #####Additive variance + family variance + residual variance for trait[i]
    Vp1<-summary(asr.model)$varcomp[1,2]+summary(asr.model)$varcomp[4,2]+summary(asr.model)$varcomp[8,2]
    
    
    ##Phenotypic variance of trait[j] 
    #####Additive variance + family variance + residual variance for trait[j]
    Vp2<-summary(asr.model)$varcomp[3,2]+summary(asr.model)$varcomp[6,2]+summary(asr.model)$varcomp[10,2]
    
    #Phenotypic correlation between trait[i] and trait[j]
    pCORR<-covP/(sqrt(Vp1*Vp2))
    
    print("Correlation checks")
    print(paste("Genetic correlation of ", traits[i]," and ",traits[j],"is = ",gCORR))
    print(paste("Environmental correlation of ", traits[i]," and ",traits[j],"is = ",eCORR))
    print(paste("Phenotypic correlation of ", traits[i]," and ",traits[j],"is = ",pCORR))
    
    # print(pin(asr.model, ~ (V2+V5+V8+V12)/sqrt((V1+V4+V7+V11)*(V3+V6+V9+V13))))
    
    ##Calculating significance of phenotypic correlation
    ## Standard error of phenotypic covariance component
    SE_P <- CovP <- summary(asr.model)$varcomp[2,3]+summary(asr.model)$varcomp[5,3]+summary(asr.model)$varcomp[9,3]
    
    ## Z Ratio of phenotypic covariance component
    z <- covP/SE_P
    
    ## Corresponding p value for Z ratio
    p <- 2*pnorm(-abs(z)) 
    
    print(paste(" p = ",p))
    

    #hij=sqrt(hi)*sqrt(hj)*gCORR
    
  }
}
sink()
