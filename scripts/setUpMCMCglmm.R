traits<-c("dep_status","CPGquant","CPG2cases","likert_total","eysenck_N","gfactor")
library("MCMCglmm")
library("MasterBayes")
# Make mcmcmglmm pedigree object with required headers and founders mo/fa=NA
mcmc_pedigree <- pedigree
names(mcmc_pedigree) <- c("animal","FATHER","MOTHER")
mcmc_pedigree$FATHER[mcmc_pedigree$FATHER==0]<-NA
mcmc_pedigree$MOTHER[mcmc_pedigree$MOTHER==0]<-NA
mcmc_pedigree<-orderPed(mcmc_pedigree)

# Change id var to animal
mcmc_dataframe<-dataframe
names(mcmc_dataframe)[1]<-"animal"
mcmc_dataframe$units<-NULL

library(MCMCglmm)
