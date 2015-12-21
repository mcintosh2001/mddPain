setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")
source("scripts/include/run_analysis.R")

traits<-c("dep_status","likert_total","eysenck_N","gfactor")

sink(paste0(locations$results,"/ace-bivar-gs/spouse_rg_rpCPG2cases",today,".txt"),append=FALSE)
cat("Bivariate correlations with CPG2cases\n")
sink()

for(i in 1:(length(traits))){
  
 
    asr.formula <- paste("fixed=cbind(CPG2cases,",traits[i],") ~trait + trait:sex + trait:age + trait:I(age^2)")
    
    
    asr.model<- asreml(fixed=as.formula(asr.formula)
                       ,random= ~us(trait,init=c(1,0.1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(spouse)
                       ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                       ,ginverse=list(id=ainv)
                       ,data=dataframe
                       ,family=asreml.gaussian(link="identity",dispersion=NA)
                       ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
    ##Assuming no genetic covariance for the null model
    asr.model.nullg<- asreml(fixed=as.formula(asr.formula)
                            ,random= ~diag(trait,init=c(1,1)):ped(id)+us(trait,init=c(1,0.1,1)):id(spouse)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
    ##Assuming no environmental covariance for the null model
    asr.model.nulle<- asreml(fixed=as.formula(asr.formula)
                            ,random= ~us(trait,init=c(1,0.1,1)):ped(id)+diag(trait,init=c(1,1)):id(spouse)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    sink(paste0(locations$results,"/ace-bivar-gs/spouse_rg_rpCPG2cases",today,".txt"),append=TRUE)
    
    print(asr.formula)
    #print(summary(asr.model)$varcomp)
    cat("\n")
    
    print(paste0("Did full model converge: ",asr.model$converge))
    print(paste0("Did null-genetic model converge: ",asr.model.nullg$converge))
    print(paste0("Did null-environment model converge: ",asr.model.nulle$converge))
    cat("\n")
    
    
    #calculating significance between model and null-genetic model  
    diffloglikg <- 2*(asr.model$loglik-asr.model.nullg$loglik)
    signifg <- pchisq(diffloglikg,1,lower.tail=FALSE)
    print(paste0("Genetic Chi2 = ",diffloglikg,", p = ",signifg))
    cat("\n")
    rg <-pin(asr.model, ~ V2/sqrt(V1*V3))[1,1]
    rg_se <- pin(asr.model, ~ V2/sqrt(V1*V3))[1,2]
    print(paste0("rg = ",rg," SE = ",rg_se))
    cat("\n")
    
    #calculating significance between model and null-environmental model  
    diffloglike <- 2*(asr.model$loglik-asr.model.nulle$loglik)
    signife <- pchisq(diffloglike,1,lower.tail=FALSE)
    print(paste0("Environment Chi2 = ",diffloglike,", p = ",signife))
    cat("\n")
    
    re <-pin(asr.model, ~ V5/sqrt(V4*V6))[1,1]
    re_se<-pin(asr.model, ~ V5/sqrt(V4*V6))[1,2]
    print(paste0("re = ",re," SE = ",re_se))
    cat("\n")
    
    rp <- pin(asr.model, ~ (V2+V5+V9)/(sqrt((V1+V4+V8)*(V3+V6+V10))))[1,1]
    rp_se <- pin(asr.model, ~ (V2+V5+V9)/(sqrt((V1+V4+V8)*(V3+V6+V10))))[1,2]
    print(paste0("rp= ",rp," SE = ", rp_se))
    
    cat("\n")
    
    #Z Ratio of phenotypic covariance component
    z <- abs(rp/rp_se)
    
  # Corresponding p value for Z ratio
   p <- 2*pnorm(z, lower.tail=FALSE) 
    
   print(paste("significance of phenotypic correlation = ",p))
  cat("\n")
  cat("\n")
  
  sink()
    
  }

