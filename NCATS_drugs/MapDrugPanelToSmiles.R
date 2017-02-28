library(synapseClient)
library(ggplot2)
library(plyr)
library(dplyr)
library(webchem)
library(ChemmineR)
synapseLogin()

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = ",", header = TRUE, fill = TRUE)

mols<-unique(drugdat$Sample.Name)

cid<-sapply(mols, function(x){
  y<-get_cid(x, from="name", first = TRUE, verbose = TRUE)
})

cid<-as.data.frame(cid)
cid$compound<-rownames(cid)
lab1<-filter(cid, !is.na(cid))
unlabeled<-filter(cid, is.na(cid))
unlabeled$compound<-sub("\\?", "", unlabeled$compound)
unlabeled$compound<-sub("\\?\\?", "", unlabeled$compound)

mols2<-unname(unlabeled$compound)

cid2<-sapply(mols2, function(x){
  y<-get_cid(x, from="name", first = TRUE, verbose = TRUE)
})

cid2<-as.data.frame(cid2)

unlabeled<-filter(cid2, is.na(cid2))
