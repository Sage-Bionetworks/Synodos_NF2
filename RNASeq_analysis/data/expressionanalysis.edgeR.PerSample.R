#!/usr/bin/Rscript

##### REQUIRED LIBRARIES #####

library(edgeR)
library(gap)

##### INPUT FILES AND PARAMETERS ####


batch1 <- "batch1.cnts"
batch2 <- "batch2.cnts"
batch3 <- "batch3.cnts"

batch2idlist <- "batch2_idlist.mod.txt"


readthreshold<-20
filter_short_genes<-"y"

shortGeneFile<-"humanGRCh37_71_ERCC_shortgenes_tRNArRNA.txt"

table1 <- read.table(file=batch1, header = F, sep = ",",skip=2)
colnames(table1) <- c("chr","start","end","gene","AC027","AC028","AC029","AC034","AC7_A17","AC7_A19","AC7_A3","AC7_A4",
"BenMen1","MIN31981","MN408a","MN408b","MN408c","MN460_MN556","MN466","MN479","MN491","MN492","MN505","MN506","MN514","MN516",
"MN520","MN521","MN527","MN529","MN533","MN548","MN560","MN563","MN567","MN571","MN572","NFT_735_MN522","xT_6491_MN474")

table1 <- table1[,c(1:4,9:12)]

colnames(table1)

table2 <- read.table(file=batch2, header = F, sep = ",",skip=2)
colnames(table2) <- c("chr","start","end","gene","MTa450-10","MTa450-11","MTa450-12","MTa450-13","MTa450-14","MTa450-15",
"MTa450-16","MTa450-17","MTa450-18","MTa450-1","MTa450-2","MTa450-3","MTa450-4",
"MTa450-5","MTa450-6","MTa450-7","MTa450-8","MTa450-9","MTa451-10","MTa451-11","MTa451-12","MTa451-13",
"MTa451-14","MTa451-15","MTa451-16","MTa451-17","MTa451-18","MTa451-1","MTa451-2","MTa451-3","MTa451-4","MTa451-5","MTa451-6","MTa451-7","MTa451-8","MTa451-9")

colnames(table2)

table3 <- read.table(file=batch3, header = F, sep = ",",skip=2)
colnames(table3) <- c("chr","start","end","gene","AC029-1","AC029-2","AC030-1","AC030-2","AC033-1",
"AC033-2","AC033-3","AC6-1","AC6-2","AC7-1","AC7-2","HS01-1","HS01-2","HS01-3","HS01-4","HS01-5",
"HS01-6","HS01-7","HS01-8","HS11-1","HS11-2","HS11-3","HS11-4","HS11-5","HS11-6","HS11-7","HS11-8",
"MN491-1","MN491-2","MN527-1","MN527-2","MN571-1","MN571-2","Syn10-1","Syn10-2","Syn10-3","Syn10-4",
"Syn10-5","Syn10-6","Syn10-7","Syn10-8","Syn2-1","Syn2-2","Syn2-3","Syn3-1","Syn3-2","Syn3-3","Syn4-1","Syn4-2","Syn4-3","Syn5-1")

#data <- data[,c(1:4,14:16,22,5,6,10:12,38:40,17:19,7:9,13,34,35,23:25,20,21,28,32,33,29,36,37,30,26,27,31)]

colnames(table3)

ids <- read.table(file=batch2idlist,header=F,sep="\t")
table2 <- table2[,match(ids$V1,colnames(table2))]
colnames(table2) <- ids$V2

colnames(table2)

table2.1 <- table2[,-c(1,2,3)]
table3.1 <- table3[,-c(1,2,3)]

tempmerged <- merge(table1,table2.1,by.x="gene",by.y="gene")
mergedtables <- merge(tempmerged,table3.1,by.x="gene",by.y="gene")


shortGenes<-read.table(shortGeneFile,sep="\t")
colnames(shortGenes)=c("ensemblID","geneSymbol","genetype")
shortGenes$gene=paste(shortGenes$geneSymbol,shortGenes$ensemblID,sep="|")


noERCC<-subset(mergedtables, !(grepl("ERCC-", mergedtables$chr)))
noShortGenes<-subset(noERCC,!(noERCC$gene %in% shortGenes$gene))
if (filter_short_genes == "y") {
  noERCC<-noShortGenes
}


rownames(noERCC) <- noERCC$gene
noERCC <- noERCC[,-c(1:4)]

##### comparisons in Table2 
##### selecting 36 samples

colnames(noERCC)

counttable2 <- noERCC[,c(5:40)]

colnames(counttable2)

groupbatch1 = factor(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12))

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

}



prefix <- "Syn1CUDC907vsSyn1DMSO"
controlCol <- c(1,2,3)
treatmentCol <- c(4,5,6)
contrast = c(-1,1,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn1PanobinostatvsSyn1DMSO"
controlCol <- c(1,2,3)
treatmentCol <- c(7,8,9)
contrast = c(-1,0,1,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn1GSK2126458vsSyn1DMSO"
controlCol <- c(1,2,3)
treatmentCol <- c(10,11,12)
contrast = c(-1,0,0,1,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5CUDC907vsSyn5DMSO"
controlCol <- c(13,14,15)
treatmentCol <- c(16,17,18)
contrast = c(0,0,0,0,-1,1,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5PanobinostatvsSyn5DMSO"
controlCol <- c(13,14,15)
treatmentCol <- c(19,20,21)
contrast = c(0,0,0,0,-1,0,1,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5GSK2126458vsSyn5DMSO"
controlCol <- c(13,14,15)
treatmentCol <- c(22,23,24)
contrast = c(0,0,0,0,-1,0,0,1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6CUDC907vsSyn6DMSO"
controlCol <- c(25,26,27)
treatmentCol <- c(28,29,30)
contrast = c(0,0,0,0,0,0,0,0,-1,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6PanobinostatvsSyn6DMSO"
controlCol <- c(25,26,27)
treatmentCol <- c(31,32,33)
contrast = c(0,0,0,0,0,0,0,0,-1,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6GSK2126458vsSyn6DMSO"
controlCol <- c(25,26,27)
treatmentCol <- c(34,35,36)
contrast = c(0,0,0,0,0,0,0,0,-1,0,0,1)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5DMSOvsSyn1DMSO"
controlCol <- c(1,2,3)
treatmentCol <- c(13,14,15)
contrast = c(-1,0,0,0,1,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6DMSOvsSyn1DMSO"
controlCol <- c(1,2,3)
treatmentCol <- c(25,26,27)
contrast = c(-1,0,0,0,0,0,0,0,1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6DMSOvsSyn5DMSO"
controlCol <- c(13,14,15)
treatmentCol <- c(25,26,27)
contrast = c(0,0,0,0,-1,0,0,0,1,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5CUDC907vsSyn1CUDC907"
controlCol <- c(4,5,6)
treatmentCol <- c(16,17,18)
contrast = c(0,-1,0,0,0,1,0,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6CUDC907vsSyn1CUDC907"
controlCol <- c(4,5,6)
treatmentCol <- c(28,29,30)
contrast = c(0,-1,0,0,0,0,0,0,0,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6CUDC907vsSyn5CUDC907"
controlCol <- c(16,17,18)
treatmentCol <- c(28,29,30)
contrast = c(0,0,0,0,0,-1,0,0,0,1,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5PanobinostatvsSyn1Panobinostat"
controlCol <- c(7,8,9)
treatmentCol <- c(19,20,21)
contrast = c(0,0,-1,0,0,0,1,0,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6PanobinostatvsSyn1Panobinostat"
controlCol <- c(7,8,9)
treatmentCol <- c(31,32,33)
contrast = c(0,0,-1,0,0,0,0,0,0,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6PanobinostatvsSyn5Panobinostat"
controlCol <- c(19,20,21)
treatmentCol <- c(31,32,33)
contrast = c(0,0,0,0,0,0,-1,0,0,0,1,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn5GSK2126458vsSyn1GSK2126458"
controlCol <- c(10,11,12)
treatmentCol <- c(22,23,24)
contrast = c(0,0,0,-1,0,0,0,1,0,0,0,0)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6GSK2126458vsSyn1GSK2126458"
controlCol <- c(10,11,12)
treatmentCol <- c(34,35,36)
contrast = c(0,0,0,-1,0,0,0,0,0,0,0,1)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn6GSK2126458vsSyn5GSK2126458"
controlCol <- c(22,23,24)
treatmentCol <- c(34,35,36)
contrast = c(0,0,0,0,0,0,0,-1,0,0,0,1)

edgeRDEanalysis(counttable2,groupbatch1,controlCol,treatmentCol,prefix,contrast,readthreshold)


# [1] "AC7-17-syn4" 1        "AC7_A19-syn5" 1       "AC7_A3-syn2" 1             
# [4] "AC7_A4-syn3" 1        "AC029-1"  2           "AC029-2" 2            
# [7] "AC030-1" 2           "AC030-2" 2            "AC033-1" 2           
#[10] "AC033-2" 2            "AC033-3" 3            "AC6-1" 4             
#[13] "AC6-2" 4              "AC7-1" 4              "AC7-2" 4             
#[16] "HS01-1-DMSO" 5        "HS01-2-DMSO" 5        "HS01-3-CUDC" 6           
#[19] "HS01-4-CUDC" 6        "HS01-5-Pen" 7          "HS01-6-Pen" 7             
##[22] "HS01-7-GSK" 8         "HS01-8-GSK" 8         "HS11-1-DMSO" 9             
#[25] "HS11-2-DMSO" 9         "HS11-3-CUDC" 10        "HS11-4-CUDC" 10             
#[28] "HS11-5-Pen" 11         "HS11-6-Pen" 11         "HS11-7-GSK" 12             
#[31] "HS11-8-GSK" 12         "MN491-1" 13            "MN491-2" 13           
#[34] "MN527-1" 13            "MN527-2" 13             "MN571-1" 13            
#[37] "MN571-2" 13            "Syn10-1-DMSO" 14        "Syn10-2-DMSO" 14            
#[40] "Syn10-3-CUDC"15        "Syn10-4-CUDC" 15        "Syn10-5-Pen" 16            
#[43] "Syn10-6-Pen" 16        "Syn10-7-GSK" 17         "Syn10-8-GSK" 17           
#[46] "Syn2-1-F" 18           "Syn2-2-F" 18           "Syn2-3" 19            
#[49] "Syn3-1-F" 20           "Syn3-2-F" 20           "Syn3-3" 19            
#[52] "Syn4-1-F" 20           "Syn4-2-F" 20           "Syn4-3" 19            
#[55] "Syn5-1" 19

groupbatch2 <- factor(c(1,1,1,1,2,2,2,2,2,2,3,4,4,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,13,13,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20,19,20,20,19,19))

counttable3 <- noERCC[,-c(5:40)]
colnames(counttable3)

prefix <- "Syn10CUDC907vsSyn10DMSO"
controlCol <- c(38,39)
treatmentCol <- c(40,41)
contrast <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0)


edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn10PanobinostatvsSyn10DMSO"
controlCol <- c(38,39)
treatmentCol <- c(42,43)
contrast <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0)


edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "Syn10GSK2126458vsSyn10DMSO"
controlCol <- c(38,39)
treatmentCol <- c(44,45)
contrast <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)





prefix <- "HS01CUDC907vsHS01DMSO"
controlCol <- c(16,17)
treatmentCol <- c(18,19)
contrast <-c(0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)


prefix <- "HS01PanobinostatvsHS01DMSO"
controlCol <- c(16,17)
treatmentCol <- c(20,21)
contrast <-c(0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS01DMSO"
controlCol <- c(16,17)
treatmentCol <- c(22,23)
contrast <-c(0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01DMSOvsHS11DMSO"
controlCol <- c(24,25)
treatmentCol <- c(16,17)
contrast <-c(0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11CUDC907vsHS11DMSO"
controlCol <- c(24,25)
treatmentCol <- c(26,27)
contrast <-c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11PanobinostatvsHS11DMSO"
controlCol <- c(24,25)
treatmentCol <- c(28,29)
contrast <-c(0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS11GSK2126458vsHS11DMSO"
controlCol <- c(24,25)
treatmentCol <- c(30,31)
contrast <-c(0,0,0,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01CUDC907vsHS11CUDC907"
controlCol <- c(26,27)
treatmentCol <- c(18,19)
contrast <-c(0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01PanobinostatvsHS11Panobinostat"
controlCol <- c(28,29)
treatmentCol <- c(20,21)
contrast <-c(0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "HS01GSK2126458vsHS11GSK2126458"
controlCol <- c(30,31)
treatmentCol <- c(22,23)
contrast <-c(0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

prefix <- "FreshSyn3andSyn4vsFreshSyn2"
controlCol <- c(46,47)
treatmentCol <- c(49,50,52,53)

contrast <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)


prefix <- "AC6AC7vsAC029AC030AC033"
controlCol <- c(5,6,7,8,9,10)
treatmentCol <- c(12,13,14,15)
contrast <-c(0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)


prefix <- "MN491MN527MN572vsAC029AC030AC033"
controlCol <- c(5,6,7,8,9,10)
treatmentCol <- c(32,33,34,35,36,37)
contrast <-c(0,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0)


edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

groupbatch2 <- factor(c(1,1,1,1,2,2,2,2,2,2,3,4,4,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,13,13,13,13,14,14,15,15,16,16,17,17,18,18,18,19,19,19,19,19,19,19))


prefix <- "Batch3Syn3andSyn4andSyn5vsBatchSyn2"
controlCol <- c(46,47,48)
treatmentCol <- c(49,50,51,52,53,54,55)
contrast <-c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1,1)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

groupbatch2 <- factor(c(1,1,18,1,2,2,2,2,2,2,3,4,4,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,13,13,13,13,14,14,15,15,16,16,17,17,18,18,18,1,1,1,1,1,1,1))

prefix <- "Batch1.3Syn3andSyn4andSyn5vsBatchSyn2"
controlCol <- c(3,46,47,48)
treatmentCol <- c(1,2,4,49,50,51,52,53,54,55)
contrast <-c(1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1)

edgeRDEanalysis(counttable3,groupbatch2,controlCol,treatmentCol,prefix,contrast,readthreshold)

sessionInfo()
