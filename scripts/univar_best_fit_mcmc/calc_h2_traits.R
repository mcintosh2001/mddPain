today <- format(Sys.Date(),format="%Y%B%d")

# Heritability of MDD

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1 <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_mdd.txt"),append=FALSE)
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))
sink()

# Heritability of MDD with sib effect

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_mdd_sib.txt"),append=FALSE)

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

# Heritability of CPGcases 

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1 <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_CPG2cases.txt"),append=FALSE)

print("Heritability of CPG2cases with no common environmental effect")
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"] +1)
print(paste0("h2MDD CPG2Cases= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))


sink()

# Heritability of CPGcases with sib effect

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_CPG2cases_sib.txt"),append=FALSE)
print("Heritability of CPG2cases with sib environmental effect")

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2MDD CPG2Cases= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD CPG2Cases= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

# Heritability of CPGquant

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1 <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal , pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_CPGquant.txt"),append=FALSE)
print("Heritability of CPGquant with no common environment")
posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"] +1)
print(paste0("h2 CPGquant= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))

sink()

# Heritability of CPGquant with sib effect

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/",today,"mcmc_h2_CPGquant_sib.txt"),append=FALSE)
print("Heritability of CPGquant with sib common environment")

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2 CPGquant= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD CPGquant= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

