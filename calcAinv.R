# Calculat Ainverse file

# First make expected file with id, fatid, matid as columns
myvars<-c("volid","father.x","mother.x")

pedigree <- pedfile[myvars]
names(pedigree)<- c("id","father","mother")

pedigree$id     <- factor(pedigree$id)
pedigree$father <- factor(pedigree$father)
pedigree$mother <- factor(pedigree$mother)


# Calculate the generalise inverse and store as an object
ainv <- asreml.Ainverse(pedigree)$ginv