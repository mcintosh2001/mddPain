setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD")
source("scripts/run_analysis.R")


# Heritability of CPG2cases with sib, spouse, old_househols, young_household effects

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002), G5 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1all <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib + spouse + old_household + young_household, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases_sib_spouse_household",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1all$VCV[, "animal"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"] +1)
print(paste0("h2CPG2cases= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1spouse <- model1.1all$VCV[, "spouse"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPG2casesspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

posterior.environability1.1old_household <- model1.1all$VCV[, "old_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPG2casesold_household= ", posterior.mode(posterior.environability1.1old_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1old_household, 0.95)))

posterior.environability1.1young_household <- model1.1all$VCV[, "young_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPG2casesyoung_household= ", posterior.mode(posterior.environability1.1young_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1young_household, 0.95)))

posterior.environability1.1sib <- model1.1all$VCV[, "sib"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2CPG2casessib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

print(paste0("model DIC for CPG2cases with animal spouse sib old_household young_houisehold =",model1.1all$DIC))
sink()




# Heritability of CPG2cases with spouse random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
prior1.1all2 <- list(G = list(G1 = list(V = 1, nu = 1000, alpha.mu=0, alpha.V=1), G2 = list(V = 1, nu = 1000, alpha.mu=0, alpha.V=1)), R = list(V = 1, fix=1))
model1.1spouse <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + spouse , pedigree = mcmc_pedigree, 
                           data = mcmc_dataframe, family="ordinal", 
                           prior = prior1.1all, nitt=4000000, thin=150, burnin=1500000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases_spouse2",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1spouse$VCV[, "animal"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"]+ model1.1spouse$VCV[, "spouse"]+ 1)
print(paste0("h2CPG2cases= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

print("1.5M burnin with 4M itterations")

posterior.environability1.1spouse <- model1.1spouse$VCV[, "spouse"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"] + model1.1spouse$VCV[, "spouse"]+1)
print(paste0("e2CPG2casesspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

summary(model1.1spouse)
print(paste0("model DIC for CPG2cases  with animal spouse  =",model1.1spouse$DIC))
sink()

pdf(paste0(locations$results,"/ace-univar-gs/",today,"final_cpg2casesPlotOfSolVCV.pdf"))
plot(model1.1spouse$Sol)
plot(model1.1spouse$VCV)
dev.off()


# Heritability of CPG2cases with animal and sib random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases_sib",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"]+ 1)
print(paste0("h2CPG2cases= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("c2CPG2casessib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))


print(paste0("model DIC for CPG2cases with animal sib  =",model1.1sib$DIC))
sink()


# Heritability of CPG2cases with animal and household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1house <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + young_household +old_household , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases_house",today,".txt"),append=FALSE)

print(paste0("model DIC for CPG2cases with animal household  =",model1.1house$DIC))
sink()


# Heritability of CPG2cases with animal and sib spouse random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sibspouse <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib + spouse , pedigree = mcmc_pedigree, 
                          data = mcmc_dataframe, family="ordinal", 
                          prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases_sibspouse",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1sibspouse$VCV[, "animal"]/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"]+ model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+ 1)
print(paste0("h2CPG2cases= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1sibspouse <- (model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"])/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"] + model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+1)
print(paste0("c2CPG2casessibspouse= ", posterior.mode(posterior.environability1.1sibspouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sibspouse, 0.95)))


print(paste0("model DIC for CPG2cases with animal sib spouse  =",model1.1sibspouse$DIC))
sink()

# Heritability of CPG2cases with animal random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1narrow <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal , pedigree = mcmc_pedigree, 
                              data = mcmc_dataframe, family="ordinal", 
                              prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_cpg2cases",today,".txt"),append=FALSE)

posterior.heritability1.1 <- model1.1narrow$VCV[, "animal"]/(model1.1narrow$VCV[,"animal"] + model1.1narrow$VCV[, "units"] + 1)
print(paste0("h2CPG2cases= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))


print(paste0("model DIC for CPG2cases with animal only  =",model1.1narrow$DIC))
sink()



