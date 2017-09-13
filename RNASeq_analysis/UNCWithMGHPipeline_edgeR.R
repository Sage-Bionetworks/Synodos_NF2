#!/usr/bin/Rscript

##### REQUIRED LIBRARIES #####

library(edgeR)
library(gap)
library(ggplot2)
library(dplyr)
library(synapseClient)
synapseLogin()

##### INPUT FILES AND PARAMETERS ####

batch1 <- synGet("syn9925491")@filePath

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

##map for hugo genes to ensembl ids for DAVID
HugoToEnsembl <- read.table(synGet("syn10639112")@filePath, sep = "\t", comment.char = "")[-1,]
colnames(HugoToEnsembl) <- c("gene", "ensembl")
HugoToEnsembl <- filter(HugoToEnsembl, ensembl != "")

counttable2 <- noERCC

##9 is HS01_DMSO_Run2_S1 - did not pass qc, so did not include
groupbatch1 = factor(c(1,1,1,2,9,2,2,3,3,3,4,4,4,5,5,5,6,6,6,6,7,7,7,8,8,8))

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

y <- DGEList(counttableselected,group=group)
y <- calcNormFactors(y,method="TMM")


print(y$samples)

design <- model.matrix(~0 + group,data=y$samples)

print(design)

pdf(file=paste(prefix,"_MDplot.logcpm.pdf",sep=""),width=15,height=8)
plotMD(cpm(y, log=TRUE),column=1)
abline(h=0,col="red",lty=2,lwd=2)
dev.off()

logcpm <-cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
pdf(file=paste(prefix,"_MDSplot.logcpm.pdf",sep=""),width=15,height=8)
par(mfrow=c(1,2))
plotMDS(logcpm,top=40000,labels=colnames(logcpm),dim.plot=c(1,2))
plotMDS(logcpm,top=40000,labels=colnames(logcpm),dim.plot=c(2,3))
dev.off()


cpm <-cpm(y,normalized.lib.sizes=TRUE, log=FALSE)
pdf(file=paste(prefix,"_MDSplot.pdf",sep=""),width=15,height=8)
par(mfrow=c(1,2))
plotMDS(cpm,top=40000,labels=colnames(cpm),dim.plot=c(1,2))
plotMDS(cpm,top=40000,labels=colnames(cpm),dim.plot=c(2,3))

dev.off()


pdf( paste(prefix,"_samplesimilarity_edgeR.pdf",sep="") , width = 7 , height = 7 )
plotMDS.DGEList(y,labels=colnames(y$counts))
dev.off()

y <- estimateDisp(y,design,robust=TRUE)
print(y$common.dispersion)

pdf(paste(prefix,"_BCVplot.pdf",sep=""), width = 7 , height = 7 )
plotBCV(y)
dev.off()

pdf(paste(prefix,"_MeanVariance.pdf",sep=""), width = 7 , height = 7 )
meanvar <- plotMeanVar(y, show.raw.vars=TRUE,show.tagwise.vars=TRUE,show.ave.raw.vars=FALSE,NBline=TRUE,main = "Mean-Variance Plot")
plotMeanVar(y, meanvar=meanvar, show.tagwise.vars=TRUE, NBline=TRUE)
dev.off()

fit <- glmQLFit(y,design)
head(fit$coefficients)

pdf(paste(prefix,"_QLDisplot.pdf",sep=""), width = 7 , height = 7 )
plotQLDisp(fit)
dev.off()

qlf <- glmQLFTest(fit,contrast=DEcontrast)
qlf$table$BH = p.adjust(qlf$table$PValue,"BH")
qlf$table$bonferroni = p.adjust(qlf$table$PValue,"bonferroni")

is.de <- decideTestsDGE(qlf, p.value=0.05)
summary(is.de)
png(filename=paste(prefix,"_edgeR_qlf_qqplot.png",sep=""),width=800,height=600)
qqunif(qlf$table$PValue)
dev.off()

png(filename=paste(prefix,"_edgeR_qlf_plotsmear.png",sep=""),width=800,height=600)
plotSmear(qlf, de.tags=rownames(qlf)[is.de!=0])
dev.off()

templist <- qlf$table
templist$gene <- rownames(templist)

templist <- left_join(templist, HugoToEnsembl)

write.table(templist,file=paste(prefix,"_edgeR_quasilikelihoodFtest.txt",sep=""),sep="\t",quote=FALSE)

ggplot(data = qlf$table, aes(x=logFC, y=-log(BH))) + 
  geom_point(data = qlf$table %>% filter(BH <= 0.1), aes(x=logFC, y=-log(BH)), color = "red", size = 1) +
  geom_point(data = qlf$table %>% filter(BH > 0.1), aes(x=logFC, y=-log(BH)), color = "black", size = 1)
ggsave(paste(prefix,"_volcanoplot.png"))

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCWithMGHPipeline_edgeR.R"

#synStore(File(paste(prefix,"_MDplot.logcpm.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_MDSplot.logcpm.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_MDSplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_samplesimilarity_edgeR.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_MeanVariance.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_BCVplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_QLDisplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_edgeR_qlf_qqplot.png",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_edgeR_qlf_plotsmear.png",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
#synStore(File(paste(prefix,"_volcanoplot.png"), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)

#synStore(File(paste(prefix,"_edgeR_quasilikelihoodFtest.txt",sep=""), parentId="syn9884467"), used = c("syn9925491","syn9884664"), executed = this.file)

}

######
#prefix <- "schwannomaCPMs"
#controlCol <- c(14,15,16)
#treatmentCol <- c(1,2,3)
#contrast = c(1,-1,0,0,-1,1,0,0,0)
#edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)
####
#QCing

prefix <- "oldanalysis"
controlCol <- c(16,17,18,19)
treatmentCol <- c(4,6,7)
contrast = c(0,1,0,0,0,-1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "newdmsoanalysis"
controlCol <- c(17,18,19,20)
treatmentCol <- c(4,6,7)
contrast = c(0,1,0,0,0,-1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

######

prefix <- "HS01DMSOvsHS11DMSO"
controlCol <- c(17,18,19,20)
treatmentCol <- c(4,6,7)
contrast = c(0,1,0,0,0,-1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01CUDC907vsHS11CUDC907"

controlCol <- c(14,15,16)
treatmentCol <- c(1,2,3)
contrast = c(1,0,0,0,-1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS11GSK2126458"

controlCol <- c(21,22,23)
treatmentCol <- c(8,9,10)
contrast = c(0,0,1,0,0,0,-1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01panobinostatvsHS11panobinostat"

controlCol <-  c(24,25,26)
treatmentCol <- c(11,12,13)
contrast = c(0,0,0,1,0,0,0,-1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)


prefix <- "HS01CUDC907vsHS01DMSO"

controlCol <- c(4,6,7)
treatmentCol <- c(1,2,3)
contrast = c(1,-1,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS01DMSO"

controlCol <- c(4,6,7)
treatmentCol <- c(8,9,10)
contrast = c(0,-1,1,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01panobinostatvsHS01DMSO"

controlCol <- c(4,6,7)
treatmentCol <- c(11,12,13)
contrast = c(0,-1,0,1,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11CUDC907vsHS11DMSO"

controlCol <- c(17,18,19,20)
treatmentCol <- c(14,15,16)
contrast = c(0,0,0,0,1,-1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11GSK2126458vsHS11DMSO"

controlCol <- c(17,18,19,20)
treatmentCol <- c(21,22,23)
contrast = c(0,0,0,0,0,-1,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11panobinostatvsHS11DMSO"

controlCol <- c(17,18,19,20)
treatmentCol <- c(24,25,26)
contrast = c(0,0,0,0,0,-1,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

sessionInfo()


getNormCounts <- function(countdata,prefix){
  y <- DGEList(countdata)
  y <- calcNormFactors(y,method="TMM")
  cpm <-cpm(y, normalized.lib.sizes=TRUE, log=TRUE)
  this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCWithMGHPipeline_edgeR.R"
  write.table(cpm,file=paste(prefix,"_edgeR_log2_cpm.txt",sep=""),sep="\t",quote=FALSE)
  synStore(File(paste(prefix,"_edgeR_log2_cpm.txt",sep=""), parentId="syn9884467"), used = c("syn9925491","syn9884664"), executed = this.file)
  
  cpm <-cpm(y, normalized.lib.sizes=TRUE, log=TRUE)
  
  #uncomment to get cpms for all genes
  cpm <- cpm(y, normalized.lib.sizes = T, log = FALSE, prior.count = 0.25)
  write.table(cpm, paste0(prefix, "_cpm.txt"), sep = "\t")
  synStore(File("schwannomaCPMs_cpm.txt", parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
 
}

prefix <- "NormCounts"
countdata <- counttable2
getNormCounts(counttable2,prefix)

#prefix <- "IntCounts"
#countdata2 <- counttable2
#getNormCounts(counttable2,prefix)

######## not enough N for interaction terms to yield useful results but this code here for reference anyway
edgeRDEanalysisWithInteraction <- function(countdata,group,control,treatment,prefix,DEcontrast,readthreshold, condition, genotype){
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
  
  design <- model.matrix(~0+condition:genotype,data=y$samples)
  print(design)
  
  pdf(file=paste(prefix,"_MDplot.logcpm.pdf",sep=""),width=15,height=8)
  plotMD(cpm(y, log=TRUE),column=1)
  abline(h=0,col="red",lty=2,lwd=2)
  dev.off()
  
  logcpm <-cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
  pdf(file=paste(prefix,"_MDSplot.logcpm.pdf",sep=""),width=15,height=8)
  par(mfrow=c(1,2))
  plotMDS(logcpm,top=40000,labels=colnames(logcpm),dim.plot=c(1,2))
  plotMDS(logcpm,top=40000,labels=colnames(logcpm),dim.plot=c(2,3))
  dev.off()
  
  
  cpm <-cpm(y,normalized.lib.sizes=TRUE, log=FALSE)
  pdf(file=paste(prefix,"_MDSplot.pdf",sep=""),width=15,height=8)
  par(mfrow=c(1,2))
  plotMDS(cpm,top=40000,labels=colnames(cpm),dim.plot=c(1,2))
  plotMDS(cpm,top=40000,labels=colnames(cpm),dim.plot=c(2,3))
  
  dev.off()
  
  
  pdf( paste(prefix,"_samplesimilarity_edgeR.pdf",sep="") , width = 7 , height = 7 )
  plotMDS.DGEList(y,labels=colnames(y$counts))
  dev.off()
  
  y <- estimateDisp(y,design,robust=TRUE)
  print(y$common.dispersion)
  
  pdf(paste(prefix,"_BCVplot.pdf",sep=""), width = 7 , height = 7 )
  plotBCV(y)
  dev.off()
  
  pdf(paste(prefix,"_MeanVariance.pdf",sep=""), width = 7 , height = 7 )
  meanvar <- plotMeanVar(y, show.raw.vars=TRUE,show.tagwise.vars=TRUE,show.ave.raw.vars=FALSE,NBline=TRUE,main = "Mean-Variance Plot")
  plotMeanVar(y, meanvar=meanvar, show.tagwise.vars=TRUE, NBline=TRUE)
  dev.off()
  
  fit <- glmQLFit(y,design)
  head(fit$coefficients)
  
  pdf(paste(prefix,"_QLDisplot.pdf",sep=""), width = 7 , height = 7 )
  plotQLDisp(fit)
  dev.off()
  
  qlf <- glmQLFTest(fit,contrast=DEcontrast)
  qlf$table$BH = p.adjust(qlf$table$PValue,"BH")
  qlf$table$bonferroni = p.adjust(qlf$table$PValue,"bonferroni")
  
  is.de <- decideTestsDGE(qlf, p.value=0.05)
  summary(is.de)
  png(filename=paste(prefix,"_edgeR_qlf_qqplot.png",sep=""),width=800,height=600)
  qqunif(qlf$table$PValue)
  dev.off()
  
  png(filename=paste(prefix,"_edgeR_qlf_plotsmear.png",sep=""),width=800,height=600)
  plotSmear(qlf, de.tags=rownames(qlf)[is.de!=0])
  dev.off()
  
  write.table(qlf$table,file=paste(prefix,"_edgeR_quasilikelihoodFtest.txt",sep=""),sep="\t",quote=FALSE)
  
  ggplot(data = qlf$table, aes(x=logFC, y=-log(BH))) + 
    geom_point(data = qlf$table %>% filter(BH <= 0.1), aes(x=logFC, y=-log(BH)), color = "red", size = 1) +
    geom_point(data = qlf$table %>% filter(BH > 0.1), aes(x=logFC, y=-log(BH)), color = "black", size = 1)
  ggsave(paste(prefix,"_volcanoplot.png"))
  
  this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCWithMGHPipeline_edgeR.R"
  
  #synStore(File(paste(prefix,"_MDplot.logcpm.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_MDSplot.logcpm.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_MDSplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_samplesimilarity_edgeR.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_MeanVariance.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_BCVplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_QLDisplot.pdf",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_edgeR_qlf_qqplot.png",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_edgeR_qlf_plotsmear.png",sep=""), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  #synStore(File(paste(prefix,"_volcanoplot.png"), parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  
  #synStore(File(paste(prefix,"_edgeR_quasilikelihoodFtest.txt",sep=""), parentId="syn9884467"), used = c("syn9925491","syn9884664"), executed = this.file)
  
}

#list1 <- c(HS01_GSK2126458_Run1_S3,  
#           HS01_GSK2126458_Run2_S3, HS01_GSK2126458_Run3_S3, HS01_panobinostat_Run1_S4, HS01_panobinostat_Run2_S4,
#           HS01_panobinostat_Run3_S4)
#list2 <-  c(HS11_GSK2126458_Run1_S7,
#            HS11_GSK2126458_Run2_S7, HS11_GSK2126458_Run3_S7, HS11_panobinostat_Run1_S8, HS11_panobinostat_Run2_S8,
#            HS11_panobinostat_Run3_S8)

prefix <- "CUDC-HS01vsHS11-AdditiveModel"
temptable <- select(counttable2, HS01_CUDC907_Run1_S2,HS01_CUDC907_Run2_S2, HS01_CUDC907_Run3_S2, 
                    HS11_CUDC907_Run1_S6, HS11_CUDC907_Run2_S6, HS11_CUDC907_Run3_S6,
                    HS01_DMSO_Run1_S1,HS01_DMSO_Run2_S1,HS01_DMSO_Run4_S9,
                    HS11_DMSO_Run1_S5,HS11_DMSO_Run2_S5,HS11_DMSO_Run3_S5,HS11_DMSO_Run4_S9)

### 1 is drug, 0 is dmso
condition <- c(rep("CUDC", 6), rep("DMSO", 7))

### 1 is mutated, 0 is WT
genotype <-  c("WT","WT","WT","MUT","MUT","MUT",
               "WT","WT","WT","MUT","MUT","MUT","MUT")
                 
ColWt <- c(7:13)
ColMut <- c(1:6)

contrast = c(1,-1,-1,1)
groupbatch2 <- c(1,1,1,2,2,2,3,3,3,4,4,4,4)

edgeRDEanalysisWithInteraction(temptable,groupbatch2,ColWt,ColMut,prefix,contrast,readthreshold, condition, genotype)

