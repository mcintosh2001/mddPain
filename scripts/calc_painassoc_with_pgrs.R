today <- format(Sys.Date(),format="%Y%B%d")

# Heritability of MDD
sink(paste0(locations$results,"/",today,"mcmc_h2_mdd.txt"),append=FALSE)

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1 <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.1 <- model1.1$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1$VCV[, "units"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1, 0.95)))
sink()

# Heritability of MDD with sib effect
sink(paste0(locations$results,"/",today,"mcmc_h2_mdd_sib.txt"),append=FALSE)

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(dep_status ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2MDD= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

# Heritability of CPGcases with sib effect
sink(paste0(locations$results,"/",today,"mcmc_h2_CPG2cases_sib.txt"),append=FALSE)

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPG2cases ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2MDD CPG2Cases= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD CPG2Cases= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

# Heritability of CPGquant with sib effect
sink(paste0(locations$results,"/",today,"mcmc_h2_CPGquant_sib.txt"),append=FALSE)

prior1.1sib <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002)), R = list(V = 1, fix=1))
model1.1sib <- MCMCglmm(CPGquant ~ sex + age +I(age^2), random = ~animal + sib, pedigree = mcmc_pedigree, 
                        data = mcmc_dataframe, family="ordinal", 
                        prior = prior1.1sib, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.1sib <- model1.1sib$VCV[, "animal"]/(model1.1$VCV[,"animal"] + model1.1sib$VCV[, "units"]+ model1.1sib$VCV[, "sib"] +1)
print(paste0("h2 CPGquant= ", posterior.mode(posterior.heritability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.heritability1.1sib, 0.95)))

posterior.environability1.1sib <- model1.1sib$VCV[, "sib"]/(model1.1sib$VCV[,"animal"] + model1.1sib$VCV[, "units"] + model1.1sib$VCV[, "sib"]+1)
print(paste0("e2MDD CPGquant= ", posterior.mode(posterior.environability1.1sib)))
print(paste0("95%CI=", HPDinterval(posterior.environability1.1sib, 0.95)))

sink()

# Effect of ldpred_mdd
sink(paste0(locations$results,"/",today,"mcmc_mddldpred.txt"),append=FALSE)

prior1.2 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,fix = 1))
model1.2 <- MCMCglmm(dep_status ~ sex + age +I(age^2) + mddldscore, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.2, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.2 <- model1.2$VCV[, "animal"]/(model1.2$VCV[,"animal"] + model1.2$VCV[, "units"] +1)
posterior.mode(posterior.heritability1.2)
HPDinterval(posterior.heritability1.2, 0.95)

print(paste0("effect of mddldpred= ", posterior.mode(model1.2$Sol[,"mddldscore"])))
print(paste0("95%CI= ",HPDinterval(model1.2$Sol[,"mddldscore"])))
sink()

# Effect of ldpred_scz
sink(paste0(locations$results,"/",today,"mcmc_sczldpred.txt") ,append=FALSE)
prior1.3 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix = 1))
model1.3 <- MCMCglmm(dep_status ~ sex + age +I(age^2) + sczldscore, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.3, nitt=65000, thin=50, burnin=15000)              

posterior.heritability1.3 <- model1.3$VCV[, "animal"]/(model1.3$VCV[,"animal"] + model1.3$VCV[, "units"] +1)
posterior.mode(posterior.heritability1.3)
HPDinterval(posterior.heritability1.3, 0.95)

print(paste0("effect of sczldpred = ", posterior.mode(model1.3$Sol[,"sczldscore"])))
print(paste0("95%CI= ",HPDinterval(model1.3$Sol[,"sczldscore"])))
sink()