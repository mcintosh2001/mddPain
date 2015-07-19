require(yaml)

# Set working directory to root of painMDD folder
setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")

# Read in files which has the absolute path to each file in it
# Includes "paindata", "pedfile", "phenotypes"
locations <- yaml.load_file("/sdata/images/projects/GENSCOT/1/andrew/painMDD/locations.yaml")

pedfile<-read.csv(locations$pedfile)
# Read pedigree file in 

load(locations$phenotypes)
totaldata$id<-factor(totaldata$id)
# Read commonly used phenotypes in

painData<-read.table(locations$paindata, header=TRUE, na.strings=-9)
names(painData)<-c("id","CPGquant","CPG2cases")  
painData$id<-factor(painData$id)

# Read pain data in and set -9s to NAs
# also change ID to id (lowercase) to match other phenotypes

# Read in MDD LDPred scores
MDDldpred <- read.table(locations$MDDpgrs, header=TRUE, na.strings=-9)
MDDldpred$mddldscore<-scale(MDDldpred$SCORE); MDDldpred$SCORE<-NULL; MDDldpred$CNT<-NULL
MDDldpred$CNT2<-NULL; MDDldpred$PHENO<-NULL; MDDldpred$FID<-NULL

BDldpred <- read.table(locations$BDpgrs, header=TRUE, na.strings=-9)
BDldpred$bdldscore<-scale(BDldpred$SCORE); BDldpred$SCORE<-NULL; BDldpred$CNT<-NULL
BDldpred$CNT2<-NULL; BDldpred$PHENO<-NULL; BDldpred$FID<-NULL

SCZldpred <- read.table(locations$SCZpgrs, header=TRUE, na.strings=-9)
SCZldpred$sczldscore<-scale(SCZldpred$SCORE); SCZldpred$SCORE<-NULL ; SCZldpred$CNT<-NULL
SCZldpred$CNT2<-NULL; SCZldpred$PHENO<-NULL; SCZldpred$FID<-NULL

allpgrs<-merge(MDDldpred,BDldpred,by="IID", all.x=TRUE, all.y=TRUE)
allpgrs<-merge(allpgrs,SCZldpred,by="IID", all.x=TRUE, all.y=TRUE)

rm(BDldpred,MDDldpred,SCZldpred)


#Merge painData with totaldata

dataframe <- merge(painData, totaldata, by="id")
dataframe <- merge(dataframe, allpgrs, by.x="id", by.y="IID") 

# Some levels of id in the dataframe are not in the pedigree?
# miss_ped <- which(!(dataframe$id %in% pedigree$id))
# Answer observations 17784 and 23596
# dataframe<-dataframe[-miss_ped]

dataframe$id<-factor(dataframe$id)
dataframe$famid<-factor(dataframe$famid)
# Remove empty factor levels in id and
# make famid a factor variable instead of integer

# Remove totaldata obejct, now replaced by dataframe
rm(totaldata)
