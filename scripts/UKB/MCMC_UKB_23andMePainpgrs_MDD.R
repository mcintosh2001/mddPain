today <- format(Sys.Date(),format="%Y%B%d")
require(MCMCglmm)

ukb_pain <-readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/health/pain.rds")

ukb_mdd <-readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/mood/mdd/putative_mdd.Rds")

ukb_painpgrs <- read.table("/sdata/images/projects/UKBIOBANK/users/mcintosh/PRSice/Pain/PRSice_SCORES_AT_ALL_THRESHOLDS.txt", header=TRUE)

ukb_demography <- readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/assessment/baseline.Rds")

ukb_pcs <- read.table("/sdata/images/projects/UKBIOBANK/incoming/uk4844_fromGail/UKB4844_UKB4723_AMcintosh_QCd_Samples_30_June_2015_GD.csv", sep=",", header=TRUE)

allData <- merge(ukb_pain, ukb_painpgrs, by.x="f.eid", by.y="IID", all.y=TRUE, all.x=T)
allData <- merge (allData,ukb_mdd, by.x="f.eid",by.y="f.eid", all.x.=T, all.y=T)
allData <- merge (allData,ukb_demography, by.x="f.eid",by.y="f.eid", all.x.=T, all.y=T)
allData <- merge (allData, ukb_pcs, by.x="f.eid", by.y="eid")
allData$chronic_group<- ordered(allData$chronic_group, levels=c("No chronic pain","Single site","Two to three sites","Widespead or multi-site pain"))

allData$chronic_group12<-as.numeric(allData$chronic_group)

allData<-subset(allData,is.na(pT_0.01)==F)

PT <- c("pT_0.01","pT_0.05","pT_0.1","pT_0.5","pT_1")

for(i in 1:length(PT)){


allData$score<-scale(allData[,PT[i]])


prior1.1all <- list(R = list(V=1, fix=1))
model1.1all <- MCMCglmm(mdd_narrow ~ sex + age + I(age^2) + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, 
                        data = allData, family="ordinal", 
                        prior = prior1.1all, nitt=212000, thin=50, burnin=12000)             

sink(paste0("/sdata/images/projects/GENSCOT/1/andrew/painMDD/results/pgrs-ukb/mcmc_pgrsPAIN_MDD_narrow",PT[i],"_",today,".txt"),append=FALSE)

print(paste("Chronic pain scores at ",PT[i]))
print(paste0("Summary of results for thtreshold ",PT[i],"on MDD narrow phenotype"))
print(summary(model1.1all))
sink()

pdf(paste0("/sdata/images/projects/GENSCOT/1/andrew/painMDD/results/pgrs-ukb/plot_mcmc_pgrsPAIN_MDDnarrow",PT[i],"_",today,".pdf"))
plot(model1.1all$Sol)
plot(model1.1all$VCV)
dev.off()

}
