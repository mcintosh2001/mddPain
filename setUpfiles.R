require(yaml)

# Set working directory to root of painMDD folder
setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")

# Read in files which has the absolute path to each file in it
# Includes "paindata", "pedfile", "phenotypes"
locations <- yaml.load_file("locations.yaml")

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

#Merge painData with totaldata

dataframe <- merge(painData, totaldata, by="id")

# Some levels of id in the dataframe are not in the pedigree?
# miss_ped <- which(!(dataframe$id %in% pedigree$id))
# Answer observations 17784 and 23596
# dataframe<-dataframe[-miss_ped]

dataframe$id<-factor(dataframe$id)
# Remove empty factor levels

# Remove totaldata obejct, now replaced by dataframe
rm(totaldata)
