rm(list=ls())
setwd("/sdata/images/projects/GENSCOT/1/andrew/painMDD/")
source("scripts/setUpfiles.R")
source("scripts/setUpvars.R")
source("scripts/calcAinv.R")
source("scripts/pin_funct.R")

# Run genetic covariance estimation with family id's included
source("scripts/calc_h2_biv_with_famid.R")
source("scripts/calc_h2_biv_with_sib.R")
