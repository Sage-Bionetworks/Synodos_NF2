#!/usr/bin/Rscript

##### REQUIRED LIBRARIES #####

library(edgeR)
library(gap)
library(ggplot2)
library(dplyr)
library(synapseClient)
synapseLogin()

##### INPUT FILES AND PARAMETERS ####

batch1 <- synGet("syn10149572")@filePath

readthreshold<-20
filter_short_genes<-"y"

shortGeneFile<- synGet("syn9884664")@filePath

table1 <- read.table(file=batch1, header = T, sep = "\t")

names(table1)[1] <- "gene"
mergedtables <- table1

shortGenes<-read.table(shortGeneFile,sep="\t")
colnames(shortGenes)=c("ensemblID","gene","genetype")

#noERCC<-subset(mergedtables, !(grepl("ERCC-", mergedtables$chr)))
noShortGenes<-subset(mergedtables,!(mergedtables$gene %in% shortGenes$gene))
if (filter_short_genes == "y") {
  noERCC<-noShortGenes
}


rownames(noERCC) <- noERCC$gene
noERCC <- noERCC[,-c(1)]

##### comparisons in Table2 
##### selecting 36 samples

counttable2 <- noERCC

##9 is HS01_DMSO_Run2_S1 - did not pass qc, so did not include
groupbatch1 = factor(c(1,3,1,1,2,2,2,2))

edgeRDEanalysis <- function(countdata,group,control,treatment,prefix,DEcontrast,readthreshold){
keep <- ((rowSums(countdata[,control] > readthreshold) == length(control) ) |  (rowSums(countdata[,treatment] > readthreshold) == length(treatment)))
cat(paste("Analyzing:",prefix,"\n",sep=""))
print(colnames(countdata[,control]))
print(colnames(countdata[,treatment]))
print(group)
cat(paste("readthreshold:",readthreshold,"\n",sep=""))
print(summary(keep))
counttableselected <- countdata[keep,]
print(colnames(counttableselected))

y <-DGEList(counttableselected,group=group)
y <- calcNormFactors(y,method="TMM")

print(y$samples)

design <- model.matrix(~0 + group,data=y$samples)

print(design)

logcpm <-cpm(y,normalized.lib.sizes=TRUE, log=TRUE)

y <- estimateDisp(y,design,robust=TRUE)
print(y$common.dispersion)

fit <- glmQLFit(y,design)
head(fit$coefficients)

qlf <- glmQLFTest(fit,contrast=DEcontrast)
qlf$table$BH = p.adjust(qlf$table$PValue,"BH")
qlf$table$bonferroni = p.adjust(qlf$table$PValue,"bonferroni")

is.de <- decideTestsDGE(qlf, p.value=0.05)
summary(is.de)


write.table(qlf$table,file=paste(prefix,"_edgeR_quasilikelihoodFtest-bedtoolsHtseq.txt",sep=""),sep="\t",quote=FALSE)

ggplot(data = qlf$table, aes(x=logFC, y=-log(BH))) + 
  geom_point(data = qlf$table %>% filter(BH <= 0.1), aes(x=logFC, y=-log(BH)), color = "red", size = 1) +
  geom_point(data = qlf$table %>% filter(BH > 0.1), aes(x=logFC, y=-log(BH)), color = "black", size = 1)
ggsave(paste(prefix,"_volcanoplot.png"))

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCWithMGHPipeline_edgeR-bedtoolsHtseq.R"

synStore(File(paste(prefix,"_edgeR_quasilikelihoodFtest-bedtoolsHtseq.txt",sep=""), parentId="syn9884467"), used = c("syn9925491","syn9884664"), executed = this.file)

}

prefix <- "HS01DMSOvsHS11DMSO"
controlCol <- c(5,6,7,8)
treatmentCol <- c(1,3,4)
contrast = c(1,-1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)






