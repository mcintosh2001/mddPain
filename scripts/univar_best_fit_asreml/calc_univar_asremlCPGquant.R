setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")
source("scripts/run_analysis.R")

# traits<-c("CPGquant","CPGquant","CPG2cases")

sink(paste0(locations$results,"/ace-univar-gs/asreml_CPGquantbestfit",today,".txt"),append=FALSE)
        
    asr.model.sib<- asreml(fixed=CPGquant~ sex + age + I(age^2)
                           ,random= ~ped(id) + id(sib) 
                           ,ginverse=list(id=ainv)
                           ,data=dataframe
                           ,maxiter=15000
                           ,trace=FALSE
                           ,na.method.X="omit"
                           ,na.method.Y="omit"
                           ,workspace=8e+08)
#     
    asr.model.spouse<- asreml(fixed=CPGquant~sex + age + I(age^2)
                              ,random= ~ped(id) + id(spouse) 
                              ,ginverse=list(id=ainv)
                              ,data=dataframe
                              ,maxiter=15000
                              ,trace=FALSE
                              ,na.method.X="omit"
                              ,na.method.Y="omit"
                              ,workspace=8e+08)
#     
    asr.model.sib.spouse<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                  ,random= ~ped(id) + id(sib) +id(spouse) 
                                  ,ginverse=list(id=ainv)
                                  ,data=dataframe
                                  ,maxiter=15000
                                  ,trace=FALSE
                                  ,na.method.X="omit"
                                  ,na.method.Y="omit"
                                  ,workspace=8e+08)

asr.model.spouse.old.young<- asreml(fixed=CPGquant~sex + age + I(age^2)
                              ,random= ~ped(id) + id(sib) +id(spouse) +id(young_household)+id(old_household)
                              ,ginverse=list(id=ainv)
                              ,data=dataframe
                              ,maxiter=15000
                              ,trace=FALSE
                              ,na.method.X="omit"
                              ,na.method.Y="omit"
                              ,workspace=8e+08)
#     
#     
    asr.model.sib.spouse.young<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                        ,random= ~ped(id) + id(sib) +id(spouse) +id(young_household)
                                        ,ginverse=list(id=ainv)
                                        ,data=dataframe
                                        ,maxiter=15000
                                        ,trace=FALSE
                                        ,na.method.X="omit"
                                        ,na.method.Y="omit"
                                        ,workspace=8e+08)

asr.model.spouse.young<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                    ,random= ~ped(id) +id(spouse) +id(young_household)
                                    ,ginverse=list(id=ainv)
                                    ,data=dataframe
                                    ,maxiter=15000,trace=FALSE,na.method.X="include",workspace=8e+08)

asr.model.spouse.old<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                    ,random= ~ped(id) +id(spouse) +id(old_household)
                              ,ginverse=list(id=ainv)
                              ,data=dataframe
                              ,maxiter=15000
                              ,trace=FALSE
                              ,na.method.X="omit"
                              ,na.method.Y="omit"
                              ,workspace=8e+08)


#     
    asr.model.sib.spouse.old<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                      ,random= ~ped(id) + id(sib) +id(spouse) +id(old_household)
                                      ,ginverse=list(id=ainv)
                                      ,data=dataframe
                                      ,maxiter=15000
                                      ,trace=FALSE
                                      ,na.method.X="omit"
                                      ,na.method.Y="omit"
                                      ,workspace=8e+08)
#     
#     
asr.model.old<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                  ,random= ~ped(id) +id(old_household)
                       ,ginverse=list(id=ainv)
                       ,data=dataframe
                       ,maxiter=15000
                       ,trace=FALSE
                       ,na.method.X="omit"
                       ,na.method.Y="omit"
                       ,workspace=8e+08)

#     
asr.model.young<- asreml(fixed=CPGquant~sex + age + I(age^2)
                                  ,random= ~ped(id)  +id(young_household)
                         ,ginverse=list(id=ainv)
                         ,data=dataframe
                         ,maxiter=15000
                         ,trace=FALSE
                         ,na.method.X="omit"
                         ,na.method.Y="omit"
                         ,workspace=8e+08)

asr.model.young.old<- asreml(fixed=CPGquant~sex + age + I(age^2)
                         ,random= ~ped(id)  +id(young_household) +id(old_household)
                         ,ginverse=list(id=ainv)
                         ,data=dataframe
                         ,maxiter=15000
                         ,trace=FALSE
                         ,na.method.X="omit"
                         ,na.method.Y="omit"
                         ,workspace=8e+08)
    
    ##Assuming no environmental covariance for the null model
    asr.model.nulle<- asreml(fixed=CPGquant ~ sex + age + I(age^2)
                            ,random= ~ped(id)
                             ,ginverse=list(id=ainv)
                             ,data=dataframe
                             ,maxiter=15000
                             ,trace=FALSE
                             ,na.method.X="omit"
                             ,na.method.Y="omit"
                             ,workspace=8e+08)
    
    #print(summary(asr.model.nulle)) 
    print(paste0("Did no.e model converge: ",asr.model.nulle$converge))
    print(paste0("Did sib only model converge: ",asr.model.sib$converge))
    print(paste0("Did spouse model converge: ",asr.model.spouse$converge))
    print(paste0("Did old model converge: ",asr.model.old$converge))
    print(paste0("Did young model converge: ",asr.model.young$converge))

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

#calculating significance between old model and null.e model  
diffloglikg <- 2*(asr.model.old$loglik-asr.model.nulle$loglik)
signif <- 1-pchisq(diffloglikg,1)
print(paste0("Significance of old effect compared to null model = ",diffloglikg,", p = ",signif))

#calculating significance between young model and null.e model  
diffloglikg <- 2*(asr.model.young$loglik-asr.model.nulle$loglik)
signif <- 1-pchisq(diffloglikg,1)
print(paste0("Significance of young effect compared to null model = ",diffloglikg,", p = ",signif))
   
    #calculating significance between sib+spouse model and sib model  
    diffloglikg <- 2*(asr.model.sib.spouse$loglik-asr.model.sib$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of sib+spouse effect compared to sib model = ",diffloglikg,", p = ",signif))
    
    #calculating significance between sib+spouse model and spouse model  
    diffloglikg <- 2*(asr.model.sib.spouse$loglik-asr.model.spouse$loglik)
    signif <- 1-pchisq(diffloglikg,1)
    print(paste0("Significance of sib+spouse effect compared to spouse model = ",diffloglikg,", p = ",signif))
    
#calculating significance between young+spouse model and spouse model  
diffloglikg <- 2*(asr.model.spouse.young$loglik-asr.model.spouse$loglik)
signif <- 1-pchisq(diffloglikg,1)
print(paste0("Significance of young+spouse effect compared to spouse model = ",diffloglikg,", p = ",signif))

#calculating significance between old+spouse model and spouse model  
diffloglikg <- 2*(asr.model.spouse.old$loglik-asr.model.spouse$loglik)
signif <- 1-pchisq(diffloglikg,1)
print(paste0("Significance of old+spouse effect compared to spouse model = ",diffloglikg,", p = ",signif))

#calculating significance between spouse+old+young model and spouse model  
diffloglik <- 2*(asr.model.spouse.old.young$loglik-asr.model.spouse$loglik)
signif <- 1-pchisq(diffloglik,2)
print(paste0("Significance of spouse+old+young effect compared to spouse model = ",diffloglikg,", p = ",signif))

#calculating significance between old+young model and nulle model  
diffloglik <- 2*(asr.model.young.old$loglik-asr.model.nulle$loglik)
signif <- 1-pchisq(diffloglik,2)
print(paste0("Significance of old + young effect compared to null.e model = ",diffloglikg,", p = ",signif))


sink()
