library(synapseClient)
library(dplyr)
library(ggplot2)
library(biomaRt)
synapseLogin()

files<- c("syn9925317",
          "syn9925318",
          "syn9925319",
          "syn9925323",
          "syn9925324",
          "syn9925325",
          "syn9925326",
          "syn9925343",
          "syn9925346",
          "syn9925348",
          "syn9925354",
          "syn9925356",
          "syn9925357",
          "syn9925358",
          "syn9925359",
          "syn9925360",
          "syn9925361",
          "syn9925362",
          "syn9925363",
          "syn9925364",
          "syn9925365",
          "syn9925366",
          "syn9925367",
          "syn9925368",
          "syn9925369",
          "syn9925371")

dat<-lapply(files, function(x){
  foo <- synGet(x)
  nam <- gsub("_sorted.bam.count", "", foo@fileHandle$fileName)
  bar <- read.table(foo@filePath, header = FALSE)
  colnames(bar) <- c("ensembl", nam)
  bar
})

mat <- dat %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ensembl"), .)

mat$ensembl <- gsub("\\..+$","",mat$ensembl)

###plot failed counts

counts<-rbind(c("total",colSums(mat[,2:27])),mat[58175:58179,])
rownames(counts) <- counts[,1]
counts <- as.data.frame(t(counts[,-1]))
names <- rownames(counts)
counts<- data.frame(apply(counts, 2, function(x) as.numeric(as.character(x))))
counts$sample<-names
ggplot(counts, aes(x=sample, y =X__no_feature/total*100)) +
  geom_bar(stat = "identity") +
  labs(x = "sample", y = "no feature (%)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("no_feature_percentage.png", width = 8, height = 5)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 mart = mart)
colnames(results) <- c("ensembl", "hugo_gene")

mat <- left_join(mat, results)
hugo_gene <- mat$hugo_gene
mat <- cbind(hugo_gene, mat[2:27])

mat <- aggregate(. ~ hugo_gene, data = mat, sum)
mat <- filter(mat, hugo_gene != "")

write.table(mat, "schwannoma_reseq_raw_counts.txt", sep = "\t")

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCExpressionMatrix.R"
synStore(File("schwannoma_reseq_raw_counts.txt", parentId = "syn9925306"), used = files, executed = this.file)
