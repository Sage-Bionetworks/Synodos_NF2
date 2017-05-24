#!/usr/bin/Rscript

##### REQUIRED LIBRARIES #####

library(edgeR)
library(gap)
library(ggplot2)
library(dplyr)

##### INPUT FILES AND PARAMETERS ####

batch1 <- "salmon_summary.human.csv_genes.txt"

readthreshold<-20
filter_short_genes<-"y"

shortGeneFile<-"humanGRCh37_71_ERCC_shortgenes_tRNArRNA.txt"

table1 <- read.table(file=batch1, header = F, sep = "\t", colClasses = c("character", rep("numeric", 26)), skip = 1)
samp <- read.table(file=batch1, header = F, sep = "\t", colClasses = c("character"))[1,]

samp[1,1] <- "gene"
colnames(table1) <- samp
table1 <- table1[-1,]
mergedtables <- table1
shortGenes<-read.table(shortGeneFile,sep="\t")

colnames(shortGenes)=c("ensemblID","gene","genetype")

noERCC<-subset(mergedtables, !(grepl("ERCC-", mergedtables$chr)))
noShortGenes<-subset(noERCC,!(noERCC$gene %in% shortGenes$gene))
if (filter_short_genes == "y") {
  noERCC<-noShortGenes
}


rownames(noERCC) <- noERCC$gene
noERCC <- noERCC[,-c(1:4)]

##### comparisons in Table2 
##### selecting 36 samples

counttable2 <- table1
rownames(counttable2)<-table1$gene
counttable2<-counttable2[,-1]

groupbatch1 = factor(c(1,1,1,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,7,8,8,8))

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
}

prefix <- "HS01DMSOvsHS11DMSO"
controlCol <- c(10,11,12,13)
treatmentCol <- c(14,15,16,17)
contrast = c(0,0,0,-1,1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01CUDC907vsHS11CUDC907"
controlCol <- c(7,8,9)
treatmentCol <- c(18,19,20)
contrast = c(0,0,-1,0,0,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS11GSK2126458"
controlCol <- c(1,2,3)
treatmentCol <- c(21,22,23)
contrast = c(-1,0,0,0,0,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01panobinostatvsHS11panobinostat"
controlCol <- c(4,5,6)
treatmentCol <- c(24,25,26)
contrast = c(0,-1,0,0,0,0,0,1)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01CUDC907vsHS01DMSO"
controlCol <- c(14,15,16,17)
treatmentCol <- c(18,19,20)
contrast = c(0,0,0,0,-1,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS01DMSO"
controlCol <- c(14,15,16,17)
treatmentCol <- c(21,22,23)
contrast = c(0,0,0,0,-1,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01panobinostatvsHS01DMSO"
controlCol <- c(14,15,16,17)
treatmentCol <- c(24,25,26)
contrast = c(0,0,0,0,-1,0,0,1)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11CUDC907vsHS11DMSO"
controlCol <- c(10,11,12,13)
treatmentCol <- c(7,8,9)
contrast = c(0,0,1,-1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11GSK2126458vsHS11DMSO"
controlCol <- c(10,11,12,13)
treatmentCol <- c(1,2,3)
contrast = c(1,0,0,-1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11panobinostatvsHS11DMSO"
controlCol <- c(10,11,12,13)
treatmentCol <- c(4,5,6)
contrast = c(0,1,0,-1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)


sessionInfo()

