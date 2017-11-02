library(synapseClient)
library(tidyr)
library(dplyr)

synapseLogin()

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

counttable2 <- noERCC

##map for hugo genes to ensembl ids for DAVID
genenames <- rownames(counttable2)

library(refGenome)
gtf <- ensemblGenome()
read.gtf(gtf, filename="Homo_sapiens_GRCh37_71_ERCC_noPatch.gtf")
genes <- gtf@ev$gtf[,c("gene_id","gene_name")] %>% distinct() %>% filter(gene_name %in% genenames)

repgenes <- table(genes$gene_name) %>% as.data.frame() %>% filter(Freq>1) 
repgenes <- repgenes$Var1

genes <- filter(genes, !gene_name %in% repgenes)

##some genes have multiple ensembl ids per gene. Extract these to get single consensus ID from hgnc downloa

map<-read.table("HugotoEnsembl.txt", sep ="\t", header = T) %>% filter(!Approved.Symbol %in% genes$gene_name) %>% 
  dplyr::rename("gene_name" = Approved.Symbol, "gene_id" = Ensembl.Gene.ID ) %>% filter(gene_id != "")
genes <- bind_rows(genes, map)

##remaining map with biomart export db
leftover <- as.data.frame(genenames) %>% filter(!genenames %in% genes$gene_name) 
bm <- read.table("mart_export.txt", sep = "\t") %>% dplyr::rename("gene_id" = V1, "gene_name" = V2) %>% 
  filter(gene_name %in% leftover$genenames)
bm$count <- c(1:nrow(bm))
bm <- bm %>% group_by(gene_name) %>% top_n(1, count) %>% ungroup() %>% select(-count)
genes <- bind_rows(genes, bm)

write.table(genes, "HugoToEnsembl.txt", sep = "\t", row.names = F)
synStore(File("HugoToEnsembl.txt", parentId = "syn9884455"), executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/HugotoEnsemblMap.R")
