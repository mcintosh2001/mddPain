today <- format(Sys.Date(),format="%Y%B%d")

source("/sdata/images/projects/GENSCOT/1/andrew/painMDD/scripts/include/run_analysis.R")

prior1.3 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix = 1))
model1.3 <- MCMCglmm(dep_status ~ sex + age +I(age^2) +C1 +C2 +C3 +C4 + MDD_pT_0.01, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.3, nitt=250000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.01.txt") ,append=FALSE)
model1.3$Fixed["formula"]
summary(model1.3)
print(paste0("effect of MDD_pT_0.01 on dep_status = ", posterior.mode(model1.3$Sol[,"MDD_pT_0.01"])))
print(paste0("95%CI= ",HPDinterval(model1.3$Sol[,"MDD_pT_0.01"])))
sink()

pdf(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.01.pdf"))
plot(model1.3$Sol)
dev.off()

prior1.3 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix = 1))
model1.3 <- MCMCglmm(dep_status ~ sex + age +I(age^2) +C1 +C2 +C3 +C4 + MDD_pT_0.05, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.3, nitt=250000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.05.txt") ,append=FALSE)
model1.3$Fixed["formula"]
summary(model1.3)
print(paste0("effect of MDD_pT_0.05 on dep_status = ", posterior.mode(model1.3$Sol[,"MDD_pT_0.05"])))
print(paste0("95%CI= ",HPDinterval(model1.3$Sol[,"MDD_pT_0.05"])))
sink()

pdf(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.05.pdf"))
plot(model1.3$Sol)
dev.off()


prior1.3 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix = 1))
model1.3 <- MCMCglmm(dep_status ~ sex + age +I(age^2) +C1 +C2 +C3 +C4 + MDD_pT_0.1, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.3, nitt=250000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.1.txt") ,append=FALSE)
model1.3$Fixed["formula"]
summary(model1.3)
print(paste0("effect of MDD_pT_0.1 on dep_status = ", posterior.mode(model1.3$Sol[,"MDD_pT_0.1"])))
print(paste0("95%CI= ",HPDinterval(model1.3$Sol[,"MDD_pT_0.1"])))
sink()

pdf(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.1.pdf"))
plot(model1.3$Sol)
dev.off()

prior1.3 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, fix = 1))
model1.3 <- MCMCglmm(dep_status ~ sex + age +I(age^2) +C1 +C2 +C3 +C4 + MDD_pT_0.5, random = ~animal, pedigree = mcmc_pedigree, 
                     data = mcmc_dataframe, family="ordinal", 
                     prior = prior1.3, nitt=250000, thin=50, burnin=15000)              

sink(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.5.txt") ,append=FALSE)
model1.3$Fixed["formula"]
summary(model1.3)
print(paste0("effect of MDD_pT_0.5 on dep_status = ", posterior.mode(model1.3$Sol[,"MDD_pT_0.5"])))
print(paste0("95%CI= ",HPDinterval(model1.3$Sol[,"MDD_pT_0.5"])))
sink()

pdf(paste0(locations$results,"/pgrs-gs/mcmc/",today,"mcmc_dep_status_MDD_pT_0.5.pdf"))
plot(model1.3$Sol)
dev.off()

