library(synapseClient)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(VennDiagram)
synapseLogin()

##for merged table 

files <- c("syn10149605")

degenes<-lapply(files, function(x){
 bar<-synGet(x)
 foo<-read.table(bar@filePath, header = T)
  comp <- bar@fileHandle$fileName
  comp <- gsub("_edgeR_quasilikelihoodFtest-bedtoolsHtseq.txt", "", comp)
  foo$comparison <- c(rep(comp, nrow(foo)))
  foo$Hugo_Gene <- rownames(foo)
  foo
})

degenes <- ldply(degenes)

##added back ENSG ids for analysis with DAVID
map <- read.table("ensemble_hugomap.txt", sep = "\t", header = TRUE)
degenes <- left_join(degenes, map)
          
this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCwithMGHPipelinePlots-bedtoolsHtseq .R"

write.table(degenes, "schwannoma_degenes_reseq_edgeR_bedtoolsHtseq.txt", sep = "\t")
synStore(File("schwannoma_degenes_reseq_edgeR_bedtoolsHtseq.txt", parentId="syn9884455"), 
         used = files,
         executed = this.file)
