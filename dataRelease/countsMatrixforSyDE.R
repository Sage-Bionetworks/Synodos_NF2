library(synapseClient)
library(tidyr)
library(dplyr)
library(Biobase)
synapseLogin()

sch <- read.table(synGet("syn9926791")@filePath, sep = "\t", header = T) 
sch$gene <- rownames(sch)

mn <- read.table(synGet("syn10595790")@filePath, sep = "\t", header = T) 
mn$gene <- rownames(mn)

cts<-full_join(sch, mn)

rownames(cts) <- cts$gene
cts<-cts %>% dplyr::select(-gene)
cts<-as.matrix(cts)

meta<-read.table(synGet("syn10831413")@filePath, sep = "\t", header = T)[1:62,]
rownames(meta) <- meta$name
meta <- meta %>% select(-name)

colnames(meta) <- c("Genotype", 'Treatment', "Cell_Type", "Cell_Line")

df<-AnnotatedDataFrame(meta)

eset<-Biobase::ExpressionSet(assayData = cts, phenoData = df)

saveRDS(eset, "NF2_eset.rds")

this.file = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/RNASeq_analysis/countsMatrixforSyDE.R"
synStore(File("NF2_eset.rds", parentId = "syn10845573"), used = c("syn9926791", "syn10595790", "syn10831413"), executed = this.file)
