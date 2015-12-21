require(yaml)

# Set working directory to root of painMDD folder
setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")

# Read in files which has the absolute path to each file in it
# Includes "paindata", "pedfile", "phenotypes"
locations <- yaml.load_file("/sdata/images/projects/GENSCOT/1/andrew/painMDD/locations.yaml")

pedfile<-read.csv(locations$pedfile)
# Read pedigree file in 

relationships<-readRDS(locations$familyData)

load(locations$phenotypes)
totaldata$id<-factor(totaldata$id)
# Read commonly used phenotypes in

painData<-read.table(locations$paindata, header=TRUE, na.strings=-9)
names(painData)<-c("id","CPGquant","CPG2cases")  
painData$id<-factor(painData$id)

# Read pain data in and set -9s to NAs
# also change ID to id (lowercase) to match other phenotypes


#Create pain pgrs predscores
PAINpgrs <- readRDS(locations$PAINpgrs)

# Read in MDD PGRS scores
MDDpgrs <- readRDS(locations$MDDpgrs)

allpgrs<-merge(MDDpgrs,PAINpgrs,by="IID", all.x=TRUE, all.y=TRUE)

# Transform each pgrs to mean=0, sd=1
allpgrs[grep("pT",names(allpgrs))] <- scale(allpgrs[grep("pT",names(allpgrs))])

rm(MDDpgrs,PAINpgrs)


#Merge painData with totaldata
dataframe <- merge(painData, totaldata, by="id")

#Merge pgrs with totaldata
dataframe <- merge(dataframe, allpgrs, by.x="id", by.y="IID") 

#Add spouse and old_household and young_household variables
dataframe <- merge(dataframe, relationships, by="id")
rm(relationships)

#Add 4 MDS components. columd 1=FID, column2=IID/id
#Columns 3-22 are the 20 MDS components (PCs)
mds.obj <- read.table(locations$mds, header=FALSE)

names(mds.obj)<-c("FID","IID",paste0("C",1:20))

mds.obj$FID<-NULL; mds.obj$id<-mds.obj$IID; mds.obj$IID<-NULL

dataframe <- merge(dataframe, mds.obj, by="id")
rm(mds.obj)

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
