setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD")
source("scripts/include/run_analysis.R")

# Narrow sense eritability of CPGquant 

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))

model1.1 <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal , pedigree = mcmc_pedigree, 
                           data = mcmc_dataframe, family="ordinal", 
                           prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant",today,".txt"),append=FALSE)

posterior.heritability1.1<- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"]+ 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))


print(paste0("model DIC for CPGquant (narrow h2) with animal   =",model1.1$DIC))
sink()



# Heritability of CPG2quant with sib, spouse, old_household, young_household random effects

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002), G5 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1all <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + sib + spouse + old_household + young_household, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant_sib_spouse_household",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1all$VCV[, "animal"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"] +1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1spouse <- model1.1all$VCV[, "spouse"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPGquantspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

posterior.environability1.1old_household <- model1.1all$VCV[, "old_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPGquantold_household= ", posterior.mode(posterior.environability1.1old_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1old_household, 0.95)))

posterior.environability1.1young_household <- model1.1all$VCV[, "young_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPGquantyoung_household= ", posterior.mode(posterior.environability1.1young_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1young_household, 0.95)))

posterior.environability1.1sib <- model1.1all$VCV[, "sib"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPGquantsib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))


print(paste0("model DIC for CPGquant  with animal spouse sib old_household young_household=",model1.1all$DIC))
sink()


# Heritability of CPG2quant with sib, spouse random effects

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sibspouse <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + sib + spouse , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant_sib_spouse",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1sibspouse$VCV[, "animal"]/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"]+ model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+ 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1spouse <- model1.1sibspouse$VCV[, "spouse"]/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"] + model1.1sibspouse$VCV[, "spouse"]+ model1.1sibspouse$VCV[, "sib"]+1)
print(paste0("e2CPGquantspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))


posterior.environability1.1sib <- model1.1sibspouse$VCV[, "sib"]/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"] + model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+1)
print(paste0("e2CPGquantsib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))


print(paste0("model DIC for CPGquant  with animal spouse sib =",model1.1sibspouse$DIC))
sink()

# Heritability of CPGquant with spouse random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1spouse <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + spouse , pedigree = mcmc_pedigree, 
                              data = mcmc_dataframe, family="ordinal", 
                              prior = prior1.1all, nitt=515000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant_spouse",today,".txt"),append=FALSE)

posterior.heritability1.1spouse <- model1.1spouse$VCV[, "animal"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"]+ model1.1spouse$VCV[, "spouse"]+ 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1spouse, 0.95)))

posterior.environability1.1spouse <- model1.1spouse$VCV[, "spouse"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"] + model1.1spouse$VCV[, "spouse"]+1)
print(paste0("c2CPGquantspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

print(paste0("model DIC for CPGquant  with animal spouse  =",model1.1spouse$DIC))

summary(model1.1spouse)

sink()

pdf(paste0(locations$results,"/ace-univar-gs/",today,"plotOfSoVCVl.pdf"))
plot(model1.1spouse$Sol)
plot(model1.1spouse$VCV)
dev.off()

# Heritability of CPG2quant with sib random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + sib , pedigree = mcmc_pedigree, 
                           data = mcmc_dataframe, family="ordinal", 
                           prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant_sib",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"]+ 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2CPGquantsib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))


print(paste0("model DIC for CPGquant with animal sib  =",model1.1sib$DIC))
sink()





# Heritability of CPG2quant with animal and household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1house <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + young_household +old_household , pedigree = mcmc_pedigree, 
                          data = mcmc_dataframe, family="ordinal", 
                          prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpgquant_house",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1house$VCV[, "animal"]/(model1.1house$VCV[,"animal"] + model1.1house$VCV[, "units"]+ model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"]+ 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1house <- (model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"])/(model1.1house$VCV[,"animal"] + model1.1house$VCV[, "units"] + model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"]+1)
print(paste0("c2CPGquanthouse= ", posterior.mode(posterior.environability1.1house)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1house, 0.95)))


print(paste0("model DIC for CPGquant with animal household  =",model1.1house$DIC))
sink()

# Heritability of CPGquant with animal and spouse young_household old_household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1spousehousehold <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + spouse + young_household+ old_household, pedigree = mcmc_pedigree, 
                                    data = mcmc_dataframe, family="ordinal", 
                                    prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_CPGquant_spouse_household",today,".txt"),append=FALSE)

print(paste0("model DIC for CPGquant with spouse+ household =",model1.1spousehousehold$DIC))
sink()

# Heritability of CPGquant with animal and spouse young_household old_household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1spousehousehold <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + spouse + young_household+ old_household, pedigree = mcmc_pedigree, 
                                    data = mcmc_dataframe, family="ordinal", 
                                    prior = prior1.1all, nitt=285000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_CPGquant_spouse_household",today,".txt"),append=FALSE)

posterior.heritability1.1spousehousehold <- model1.1spousehousehold$VCV[, "animal"]/(model1.1spousehousehold$VCV[,"animal"] + model1.1spousehousehold$VCV[, "units"]+ model1.1spousehousehold$VCV[, "young_household"] + model1.1spousehousehold$VCV[, "old_household"]+ model1.1spousehousehold$VCV[, "spouse"] + 1)
print(paste0("h2CPGquant= ", posterior.mode(posterior.heritability1.1spousehousehold)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1spousehousehold, 0.95)))

posterior.environability1.1household <- (model1.1spousehousehold$VCV[, "young_household"]+ model1.1spousehousehold$VCV[, "old_household"])/(model1.1spousehousehold$VCV[,"animal"] + model1.1spousehousehold$VCV[, "units"] + model1.1spousehousehold$VCV[, "young_household"]+ model1.1spousehousehold$VCV[, "old_household"] + model1.1spousehousehold$VCV[,"spouse"] +1)
print(paste0("c2CPGquanthousehold= ", posterior.mode(posterior.environability1.1household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1household, 0.95)))

posterior.environability1.1spouse <- (model1.1spousehousehold$VCV[, "spouse"])/(model1.1spousehousehold$VCV[,"animal"] + model1.1spousehousehold$VCV[, "units"] + model1.1spousehousehold$VCV[, "young_household"]+ model1.1spousehousehold$VCV[, "old_household"]+ model1.1spousehousehold$VCV[, "spouse"]+1)
print(paste0("c2CPGquantspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

print(paste0("model DIC for CPGquant with spouse+ household =",model1.1spousehousehold$DIC))
sink()


