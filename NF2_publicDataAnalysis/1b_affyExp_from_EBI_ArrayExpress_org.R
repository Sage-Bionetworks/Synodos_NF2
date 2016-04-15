#script to process Expression data from EBI Array Express

#source("http://bioconductor.org/biocLite.R")
#biocLite("ArrayExpress")

library("ArrayExpress")

study <- ArrayExpress("E-MEXP-2766")

phenoData(study)