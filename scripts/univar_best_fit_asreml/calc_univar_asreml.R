setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")
source("scripts/include/run_analysis.R")

# traits<-c("dep_status","CPGquant","CPG2cases")

sink(paste0(locations$results,"/ace-univar-gs/asreml_MDDbestfit",today,".txt"),append=FALSE)
        
    asr.model.sib<- asreml(fixed=dep_status~ sex + age + I(age^2)
                           ,random= ~ped(id) + id(sib) 
                           ,ginverse=list(id=ainv)
                           ,data=dataframe
                           ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
#     
    asr.model.spouse<- asreml(fixed=dep_status~sex + age + I(age^2)
                              ,random= ~ped(id) + id(spouse) 
                              ,ginverse=list(id=ainv)
                              ,data=dataframe
                              ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
#     
    asr.model.sib.spouse<- asreml(fixed=dep_status~sex + age + I(age^2)
                                  ,random= ~ped(id) + id(sib) +id(spouse) 
                                  ,ginverse=list(id=ainv)
                                  ,data=dataframe
                                  ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
#     
#     
    asr.model.sib.spouse.young<- asreml(fixed=dep_status~sex + age + I(age^2)
                                        ,random= ~ped(id) + id(sib) +id(spouse) +id(young_household)
                                        ,ginverse=list(id=ainv)
                                        ,data=dataframe
                                        ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
#     
    asr.model.sib.spouse.old<- asreml(fixed=dep_status~sex + age + I(age^2)
                                      ,random= ~ped(id) + id(sib) +id(spouse) +id(old_household)
                                      ,ginverse=list(id=ainv)
                                      ,data=dataframe
                                      ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
#     
    
    ##Assuming no environmental covariance for the null model
    asr.model.nulle<- asreml(fixed=dep_status ~ sex + age + I(age^2)
                            ,random= ~ped(id)
                            ,ginverse=list(id=ainv)
                            ,data=dataframe
                            ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)
    
     print(summary(asr.model.nulle)) 
    print(paste0("Did no.e model converge: ",asr.model.nulle$converge))
    print(paste0("Did sib only model converge: ",asr.model.sib$converge))
    print(paste0("Did spouse model converge: ",asr.model.spouse$converge))
    print(paste0("Did sib and spouse model converge: ",asr.model.sib.spouse$converge))
    print(paste0("Did sib spouse and old model converge: ",asr.model.sib.spouse.old$converge))
    print(paste0("Did sib spouse and young model converge: ",asr.model.sib.spouse.young$converge))
    
    
    

    
    #calculating significance between sib model and null.e model  
    diffloglikg <- 2*(asr.model.sib$loglik-asr.model.nulle$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of sib effect compared to null model = ",diffloglikg,", p = ",signif))

    #calculating significance between spouse model and null.e model  
    diffloglikg <- 2*(asr.model.spouse$loglik-asr.model.nulle$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of spouse effect compared to null model = ",diffloglikg,", p = ",signif))
   
    #calculating significance between sib+spouse model and sib model  
    diffloglikg <- 2*(asr.model.sib.spouse$loglik-asr.model.sib$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of sib+spouse effect compared to sib model = ",diffloglikg,", p = ",signif))
    
    #calculating significance between sib+spouse model and spouse model  
    diffloglikg <- 2*(asr.model.sib.spouse$loglik-asr.model.spouse$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of sib+spouse effect compared to spouse model = ",diffloglikg,", p = ",signif))
    
    


sink()
