# Heritability of CPG2cases with sib, spouse, old_househols, young_household effects

p.var<-var(mcmc_dataframe$CPG2cases,na.rm=TRUE)
prior1.1all <- list(G = list(G1 = list(V = 1, nu = 1), G2 = list(V = 1, nu = 1), G3 = list(V = 1, nu = 1), G4 = list(V = 1, nu = 1), G5 = list(V = 1, nu = 1)), R = list(V = 2, fix=1))
model1.1all <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib + spouse + old_household + young_household, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="threshold", 
                        prior = prior1.1all, nitt=512000, thin=100, burnin=12000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_sib_spouse_household",today,".txt"),append=FALSE)

summary(model1.1all)

print(paste0("DIC for sib spouse household = ", model1.1all$DIC))
print(autocorr(model1.1all$VCV))

sink()

pdf(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_sib_spouse_household",today,".pdf"))
plot(model1.1all$VCV)
plot(model1.1all$Sol)
dev.off()

# Heritability of CPG2cases with sib, spouse,  effects

p.var<-var(mcmc_dataframe$CPG2cases,na.rm=TRUE)
prior1.1all <- list(G = list(G1 = list(V = 1, nu = 1), G2 = list(V = 1, nu = 1), G3 = list(V = 1, nu = 1)), R = list(V = 2, fix=1))
model1.1all <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib + spouse , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="threshold", 
                        prior = prior1.1all, nitt=512000, thin=100, burnin=12000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_sib_spouse",today,".txt"),append=FALSE)

summary(model1.1all)

print(paste0("DIC for sib spouse  = ",model1.1all$DIC))
print(autocorr(model1.1all$VCV))

        sink()
        
pdf(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_sib_spouse",today,".pdf"))
        plot(model1.1all$VCV)
        plot(model1.1all$Sol)
dev.off()
        
# Heritability of CPG2cases with spouse,  effects
        
        p.var <- var(mcmc_dataframe$CPG2cases,na.rm=TRUE)
        
prior1.1all <- list(G = list(G1 = list(V = 40, nu = 1), G2 = list(V = 40, nu = 1)), R = list(V = 2, fix=1))
model1.1all <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + spouse , pedigree = mcmc_pedigree, 
                                data = mcmc_dataframe, family="threshold", 
                                prior = prior1.1all, nitt=212000, thin=100, burnin=12000)              
        
        sink(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_spouse",today,".txt"),append=FALSE)
        
        summary(model1.1all)
        
        print(paste0("DIC for spouse  = ",model1.1all$DIC))
        print(autocorr(model1.1all$VCV))
        
                sink()
                
                pdf(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_spouse",today,".pdf"))
                plot(model1.1all$VCV)
                plot(model1.1all$Sol)
                dev.off()
        
        # Heritability of CPG2cases with sib,  effects
        
        p.var<-var(mcmc_dataframe$CPG2cases,na.rm=TRUE)
        prior1.1all <- list(G = list(G1 = list(V = 1, nu = 1), G2 = list(V = 1, nu = 1)), R = list(V = 2, fix=1))
                            model1.1all <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib , pedigree = mcmc_pedigree, 
                                                    data = mcmc_dataframe, family="threshold", 
                                                    prior = prior1.1all, nitt=512000, thin=100, burnin=12000)              
                            
                            sink(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_spouse",today,".txt"),append=FALSE)
                            
                            summary(model1.1all)
                            
                            print(paste0("DIC for sib  = ",model1.1all$DIC))
                            print(autocorr(model1.1all$VCV))
                            
                            sink()
                            
                            pdf(paste0(locations$results,"/ace-univar-gs/mcmc_THRESH_cpg2cases_sib",today,".pdf"))
                            plot(model1.1all$VCV)
                            plot(model1.1all$Sol)
                            dev.off()