library(synapseClient)
library(tidyr)
library(dplyr)
library(Biobase)
synapseLogin()

meta<-read.table(synGet("syn10831413")@filePath, sep = "\t", header = T)
rownames(meta) <- gsub("\\_", "\\.",meta$name)
meta <- meta %>% select(-name)

colnames(meta) <- c("Genotype", 'Treatment', "Cell_Type", "Cell_Line")

sch <- read.table(synGet("syn10882987")@filePath, sep = "\t", header = T) %>% 
  rename("gene" = "Gene")
sch <- sch %>% aggregate(. ~ gene, data = ., sum)

mn <- read.table(synGet("syn10882903")@filePath, sep = "\t", header = T)[,c(1,11:46)]
mn$gene <- sapply(mn$gene, function(x) unlist(strsplit(x, split = "\\|"))[1])
mn <- mn %>% aggregate(. ~ gene, data = ., sum)

cts<-inner_join(sch, mn)
colnames(cts) <- gsub("\\_", "\\.", colnames(cts))
colnames(cts) <- as.character(gsub("CUDC\\.907", "CUDC907", colnames(cts)))

sampls <- rownames(meta)
genes <- cts$gene
cts<-cts %>% dplyr::select(one_of(c(sampls))) %>% mutate_all(funs(min_rank(.)))

cts<-as.matrix(cts) %>% scale(center = TRUE, scale = TRUE)
rownames(cts) <- genes

meta<-Biobase::AnnotatedDataFrame(meta)

eset<-Biobase::ExpressionSet(assayData = cts, phenoData = meta)

saveRDS(eset, "NF2_eset.rds")

this.file = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/RNASeq_analysis/countsMatrixforSyDE.R"
synStore(File("NF2_eset.rds", parentId = "syn10845573"), used = c("syn9926791", "syn10595790", "syn10831413"), executed = this.file)
