##Perform univariate and bivariate (with depression) heritability on traits

univtraits<-c("CPGquant","CPG2cases")

traits<-c("dep_status","CPGquant","CPG2cases")
             
## Univariate heritability of traits

for (i in 1:length(univtraits)){
sink(paste0(univtraits[i],"heritability.txt"))
formula <- paste("fixed=",univtraits[i]," ~1 + sex + age + age2")

print(formula)
asr.model<- asreml(fixed= as.formula(formula), random =~ ped(id), ginverse=list(id=ainv), data=pain_mdd, maxiter=500,trace=FALSE,na.method.X="omit")

## Heritability estimate
h2<- (summary(asr.model)$varcomp[1,2])/(summary(asr.model)$varcomp[1,2]+summary(asr.model)$varcomp[2,2])

print(summary(asr.model)$varcomp)
print(paste("Heritability of", univtraits[i]," =", h2))

## Get SE for heritability estimate
print(pin(asr.model, hsw ~ V1/(V1+V2)))
sink()
}



## Bivariate heritability between traits (genetic and phenotypic correlations)
## Genetic correlation

sink("phengenCorr.txt",append=TRUE)
for(i in 1:(length(traits) - 1)){
  for(j in (i+1):length(traits)){
    
    pain_mdd$test1<-pain_mdd[,traits[i]]
    pain_mdd$test2<-pain_mdd[,traits[j]]
    
    asr.model<- asreml(cbind(test1,test2)~trait+trait:sex+trait:age+trait:age2
                       ,random= ~us(trait,init=c(0.1,0.3,0.1)):ped(id)
                       ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                       ,ginverse=list(id=ainv)
                       ,data=pain_mdd
                       ,family=asreml.gaussian(link="identity",dispersion=NA)
                       ,maxiter=15000,trace=FALSE,na.method.X="include")
    
    # No covariance
    asr.model.null<- asreml(cbind(test1,test2)~trait+trait:sex+trait:age+trait:age2
                            ,random= ~diag(trait,init=c(1,1)):ped(id)
                            ,rcov= ~units:us(trait,init=c(0.1,0.1,0.1))
                            ,ginverse=list(id=ainv)
                            ,data=pain_mdd
                            ,family=asreml.gaussian(link="identity",dispersion=NA)
                            ,maxiter=15000,trace=FALSE,na.method.X="include")
    
    	summary(asr.model)$varcomp
    
    	#Additive covariance
    	covA<-summary(asr.model)$varcomp[2,2]
   
    	#Additive variance of trait[i] 
    	Vg1A<-summary(asr.model)$varcomp[1,2]
    
    	#Additive variance of trait[j]
    	Vg2A<-summary(asr.model)$varcomp[3,2]
    
    	#Genetic correlation between trait[i] and trait[j]
    	gCOV<-covA/(sqrt(Vg1A*Vg2A))
    
    	print(paste("Genetic correlation of ", traits[i]," and ",traits[j],"is = ",gCOV))
    
    	#calculating significance between model and null model	
    	diffloglik <- 2*(asr.model$loglik-asr.model.null$loglik)
    	signif <- 1-pchisq(diffloglik,1)
    
    	print(paste0("Chi2 = ",diffloglik,", p = ",signif))
    
    	print(pin(asr.model, ~ V2/sqrt(V1*V3)))


##Phenotypic correlation
	summary(asr.model)$varcomp
    
	#Additive covariance	
	covA<-summary(asr.model)$varcomp[2,2]
	
	#Residual covariance
	covR<-summary(asr.model)$varcomp[6,2]
	
	#Phenotypic covariance
	covP<-covA+covR
	
	#Additive variance of trait[i]
	Vp1A<-summary(asr.model)$varcomp[1,2]

	#Residual variance of trait[i]
	Vp1R<-summary(asr.model)$varcomp[5,2]

	#Phenotypic variance for trait[i]
	Vp1<-Vp1A+Vp1R

	#Additive variance of trait[j]
	Vp2A<-summary(asr.model)$varcomp[3,2]

	#Residual variance of trait[j]
	Vp2R<-summary(asr.model)$varcomp[7,2]

	#Phenotypic variance for trait[j]
	Vp2<-Vp2A+Vp2R
    
	#Phenotypic correlation between trait[i] and trait[j]
	pCOV<-covP/(sqrt(Vp1*Vp2))
    
	
	print(paste("Phenotypic correlation of ", traits[i]," and ",traits[j],"is = ",pCOV))

  	print(pin(asr.model, ~ (V2+V6)/sqrt((V1+V5)*(V3+V7))))

##Calculating significance of phenotypic correlation

	## Standard error of additive covariance component
	SE_A<-summary(asr.model)$varcomp[2,3]

	## Standard error of residual covariance component
	SE_R<-summary(asr.model)$varcomp[6,3]

	## Standard error of phenotypic covariance component
	SE_P<-SE_A+SE_R

	## Z Ratio of phenotypic covariance component
	z <- covP/SE_P

	## Corresponding p value for Z ratio
	p <- 2*pnorm(-abs(z)) 
        
	print(paste(" p = ",p))


## Remove model generated traits    
	pain_mdd$test1 <- NULL
	pain_mdd$test2 <- NULL
          
  }
}
sink()
