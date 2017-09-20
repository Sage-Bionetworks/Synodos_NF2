#!/usr/bin/Rscript

##### REQUIRED LIBRARIES #####

library(edgeR)
library(gap)

##### INPUT FILES AND PARAMETERS ####

## all necessary files can be found in syn6038272
batch1 <- "data/batch1.cnts"
batch2 <- "data/batch2.cnts"
batch3 <- "data/batch3.cnts"

batch2idlist <- "data/batch2_idlist.mod.txt"

#args<-commandArgs(TRUE)
#countfile1 <- args[1]

readthreshold<-20
filter_short_genes<-"y"

shortGeneFile <- "data/humanGRCh37_71_ERCC_shortgenes_tRNArRNA.txt"

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

counttable2$gene <- rownames(counttable2)

counttable2 <- counttable2 %>% 
  separate(gene, c("Hugo_Gene", "ensembl"), sep = '\\|') %>% 
  dplyr::select(-ensembl)

counttable2 <- aggregate(. ~ Hugo_Gene, data = counttable2, sum)

rownames(counttable2) <- counttable2$Hugo_Gene
counttable2 <- select(counttable2, -Hugo_Gene)

getNormCounts <- function(countdata,prefix){
  y <- DGEList(countdata)
  y <- calcNormFactors(y,method="TMM")
  cpm <-cpm(y, normalized.lib.sizes=TRUE, log=TRUE)
  #this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCWithMGHPipeline_edgeR.R"
  #write.table(cpm,file=paste(prefix,"_edgeR_log2_cpm.txt",sep=""),sep="\t",quote=FALSE)
  #synStore(File(paste(prefix,"_edgeR_log2_cpm.txt",sep=""), parentId="syn9884467"), used = c("syn9925491","syn9884664"), executed = this.file)
  
  cpm <-cpm(y, normalized.lib.sizes=TRUE, log=TRUE)
  
  #uncomment to get cpms for all genes
  cpm <- cpm(y, normalized.lib.sizes = T, log = TRUE, prior.count = 0.25)
  write.table(cpm, paste0(prefix, "_cpm.txt"), sep = "\t")
  #synStore(File("schwannomaCPMs_cpm.txt", parentId="syn9884466"), used = c("syn9925491","syn9884664"), executed = this.file)
  
}

countdata <- counttable2
prefix <- "meningiomaCPM"

getNormCounts(counttable2,prefix)

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/meningiomaCPMs.R"
synStore(File("meningiomaCPM_cpm.txt", parentId="syn9884466"), executed = this.file)


