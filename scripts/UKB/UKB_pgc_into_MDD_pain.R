ukb_pain <-readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/health/pain.rds")

ukb_mdd <-readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/mood/mdd/putative_mdd.Rds")

ukb_mddpgrs <- read.table("/sdata/images/projects/UKBIOBANK/incoming/PGRS/UKB4844_MDD_SPH_29082015_SCORES_AT_ALL_THRESHOLDS.txt", header=TRUE)
ukb_mddpgrs$pT_0.5<-scale(ukb_mddpgrs$pT_0.5)

ukb_demography <- readRDS("/sdata/images/projects/UKBIOBANK/data/phenotypes/assessment/baseline.Rds")

ukb_pcs <- read.table("/sdata/images/projects/UKBIOBANK/incoming/uk4844_fromGail/UKB4844_UKB4723_AMcintosh_QCd_Samples_30_June_2015_GD.csv", sep=",", header=TRUE)

allData <- merge(ukb_pain, ukb_mddpgrs, by.x="f.eid", by.y="IID", all.y=TRUE, all.x=T)
allData <- merge (allData,ukb_mdd, by.x="f.eid",by.y="f.eid", all.x.=T, all.y=T)
allData <- merge (allData,ukb_demography, by.x="f.eid",by.y="f.eid", all.x.=T, all.y=T)
allData <- merge (allData, ukb_pcs, by.x="f.eid", by.y="eid")
allData$chronic_group<- ordered(allData$chronic_group, levels=c("No chronic pain","Single site","Two to three sites","Widespead or multi-site pain"))

allData$chronic_group12<-as.numeric(allData$chronic_group)

allData<-subset(allData,is.na(pT_0.01)==F)

PT <- c("pT_0.01","pT_0.05","pT_0.1","pT_0.5","pT_1")

for(i in 1:length(PT)){
  
sink(paste0("/sdata/images/projects/GENSCOT/1/andrew/painMDD/results/pgrs-ukb/pgrsALL_",PT[i],"_",today,".txt"),append=FALSE)

print(paste0("Summary of results for thtreshold ",PT[i]))

allData$score<-scale(allData[,PT[i]])

print("Chronic pain yes/no")
test_pain <- glm( I(allData$chronic_group!="No chronic pain") ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, family=binomial(link="logit"), data=allData)
test_pain_null <- glm( I(allData$chronic_group!="No chronic pain") ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre, family=binomial(link="logit"), data=allData)
print(summary(test_pain))
var_explained <- 1-(logLik(test_pain)/logLik(test_pain_null))

print(paste("variance explained by ",PT[i]," in pain is ",var_explained[1]))

print("Chronic pain quantitative")
test_painq <- glm( chronic_group12 ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, family=gaussian(link="identity"), data=allData)
test_painq_null <- glm( chronic_group12 ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre, family=gaussian(link="identity"), data=allData)
print(summary(test_painq))
var_explained <- 1-(logLik(test_painq)/logLik(test_painq_null))



print(paste("variance explained by ",PT[i]," in painq is ",var_explained[1]))

print("MDD broad yes/no")
test_mdd_broad <- glm( allData$mdd_broad ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre+ score, family=binomial(link="logit"), data=allData)
test_mdd_broad_null <- glm( allData$mdd_broad ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre, family=binomial(link="logit"), data=allData)
print(summary(test_mdd_broad))
var_explained <- 1-(logLik(test_mdd_broad)/logLik(test_mdd_broad_null))

print(paste("variance explained by ",PT[i]," in mdd_broad is ",var_explained[1]))

print("MDD narrow yes/no")
test_mdd_narrow <- glm( allData$mdd_narrow ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, family=binomial(link="logit"), data=allData)
test_mdd_narrow_null <- glm( allData$mdd_narrow ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre , family=binomial(link="logit"), data=allData)
print(summary(test_mdd_narrow))
var_explained <- 1-(logLik(test_mdd_narrow)/logLik(test_mdd_narrow_null))

print(paste("variance explained by ",PT[i]," in mdd_narrow is ",var_explained[1]))

print("MDD strict yes/no")
test_mdd_strict <- glm( mdd_strict ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, family=binomial(link="logit"), data=allData)
test_mdd_strict_null <- glm( mdd_strict ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre, family=binomial(link="logit"), data=allData)
print(summary(test_mdd_strict))
var_explained <- 1-(logLik(test_mdd_strict)/logLik(test_mdd_strict_null))

print(paste("variance explained by ",PT[i]," in mdd_strict is ", var_explained[1]))

print("MDD self yes/no")
test_mdd_self <- glm( depr_self ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre + score, family=binomial(link="logit"), data=allData)
test_mdd_self_null <- glm( depr_self ~ age + sex + pc1 + pc2 + pc3 +pc4 + pc5 +pc6+ pc7+ pc8 +pc9 +pc10 +gen_batch +array +assessment_centre, family=binomial(link="logit"), data=allData)
print(summary(test_mdd_self))
var_explained <- 1-(logLik(test_mdd_self)/logLik(test_mdd_self_null))

print(paste("variance explained by ",PT[i]," in mdd_self is ", var_explained[1]))

sink()
}
