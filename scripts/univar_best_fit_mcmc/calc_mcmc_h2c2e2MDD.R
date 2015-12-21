setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD")
source("scripts/run_analysis.R")
# Crude narrow sense heritability of MDD

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1 <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd",today,".txt"),append=FALSE)
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))

print(paste0("model DIC for MDD  with animal only =",model1.1$DIC))

sink()

# Heritability of MDD with sib effect

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1sib, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_sib",today,".txt"),append=FALSE)

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

print(paste0("model DIC for MDD with animal and sib =",model1.1sib$DIC))
sink()

# Heritability of MDD with spouse effect

prior1.1spouse <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1spouse <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + spouse, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1spouse, nitt=515000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_spouse",today,".txt"),append=FALSE)

posterior.heritability1.1spouse <- model1.1spouse$VCV[, "animal"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"]+ model1.1spouse$VCV[, "spouse"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1spouse, 0.95)))

posterior.environability1.1spouse <- model1.1spouse$VCV[, "spouse"]/(model1.1spouse$VCV[,"animal"] + model1.1spouse$VCV[, "units"] + model1.1spouse$VCV[, "spouse"]+1)
print(paste0("c2MDD= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

print(paste0("model DIC for MDD with animal and spouse =",model1.1spouse$DIC))

summary(model1.1spouse)
sink()

pdf(paste0(locations$results,"/ace-univar-gs/",today,"FinalMDDPlots.pdf"))
plot(model1.1spouse$Sol)
plot(model1.1spouse$VCV)
dev.off()

# Heritability of MDD with sib, spouse, old_househols, young_household effects

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002), G5 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1all <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + sib + spouse + old_household + young_household, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_sib_spouse_household",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1all$VCV[, "animal"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1spouse <- model1.1all$VCV[, "spouse"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2MDDspouse= ", posterior.mode(posterior.environability1.1spouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1spouse, 0.95)))

posterior.environability1.1old_household <- model1.1all$VCV[, "old_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2MDDold_household= ", posterior.mode(posterior.environability1.1old_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1old_household, 0.95)))

posterior.environability1.1young_household <- model1.1all$VCV[, "young_household"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2MDDyoung_household= ", posterior.mode(posterior.environability1.1young_household)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1young_household, 0.95)))

posterior.environability1.1sib <- model1.1all$VCV[, "sib"]/(model1.1all$VCV[,"animal"] + model1.1all$VCV[, "units"] + model1.1all$VCV[, "sib"]+ model1.1all$VCV[, "spouse"]+ model1.1all$VCV[, "old_household"]+ model1.1all$VCV[, "young_household"]+1)
print(paste0("e2MDDsib= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))


print(paste0("model DIC for MDD with animal sib spouse old_household young_household =",model1.1all$DIC))
sink()


# Heritability of dep_status with animal and household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1house <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + young_household +old_household , pedigree = mcmc_pedigree, 
                          data = mcmc_dataframe, family="ordinal", 
                          prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_house",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1house$VCV[, "animal"]/(model1.1house$VCV[,"animal"] + model1.1house$VCV[, "units"]+ model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"]+ 1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1house <- (model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"])/(model1.1house$VCV[,"animal"] + model1.1house$VCV[, "units"] + model1.1house$VCV[, "young_household"]+ model1.1house$VCV[, "old_household"]+1)
print(paste0("c2MDDhouse= ", posterior.mode(posterior.environability1.1house)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1house, 0.95)))


print(paste0("model DIC for MDD with animal household  =",model1.1house$DIC))
sink()

# Heritability of dep_status with animal and sib and spouse random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sibspouse <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + sib + spouse , pedigree = mcmc_pedigree, 
                          data = mcmc_dataframe, family="ordinal", 
                          prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_sibspouse",today,".txt"),append=FALSE)

posterior.heritability1.1all <- model1.1sibspouse$VCV[, "animal"]/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"]+ model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+ 1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1all)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1all, 0.95)))

posterior.environability1.1sibspouse <- (model1.1sibspouse$VCV[, "sib"]+model1.1sibspouse$VCV[, "spouse"])/(model1.1sibspouse$VCV[,"animal"] + model1.1sibspouse$VCV[, "units"] + model1.1sibspouse$VCV[, "sib"]+ model1.1sibspouse$VCV[, "spouse"]+1)
print(paste0("c2MDDsibspouse= ", posterior.mode(posterior.environability1.1sibspouse)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sibspouse, 0.95)))


print(paste0("model DIC for MDD with sib + spouse  =",model1.1sibspouse$DIC))
sink()

# Heritability of dep_status with animal and spouse young_household old_household random effect

prior1.1all <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1spousehousehold <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + spouse + young_household+ old_household, pedigree = mcmc_pedigree, 
                              data = mcmc_dataframe, family="ordinal", 
                              prior = prior1.1all, nitt=85000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/ace-univar-gs/mcmc_h2_mdd_spouse_household",today,".txt"),append=FALSE)

print(paste0("model DIC for MDD with spouse+ household =",model1.1spousehousehold$DIC))
sink()


