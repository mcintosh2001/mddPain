# Make sib variabe to identify pairs of siblings
dataframe$sib<-paste0(dataframe$father,"-",dataframe$mother)
dataframe$sib[dataframe$sib=="0-0"]<-NA
dataframe$sib[is.na(dataframe$sib)]<-dataframe$famid[is.na(dataframe$sib)]
dataframe$sib <- factor(dataframe$sib)

# Make gfactor
require(psych)
dataframe$memorytotal<-dataframe$logical_mem_1+dataframe$logical_mem_2

dataframe$zmemorytotal<-scale(dataframe$memorytotal)
dataframe$zvocabulary<-scale(dataframe$vocabulary)
dataframe$zverbal_total<-scale(dataframe$verbal_total)
dataframe$zdigit_symbol<-scale(dataframe$digit_symbol)


myvars<-c("zverbal_total","zvocabulary","zmemorytotal","zdigit_symbol")


dataframe$gfactor <- principal(dataframe[,myvars], nfactors=1, rotate="none")$scores
