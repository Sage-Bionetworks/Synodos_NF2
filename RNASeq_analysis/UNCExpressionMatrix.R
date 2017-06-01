library(synapseClient)
library(dplyr)
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

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 mart = mart)
colnames(results) <- c("ensembl", "hugo_gene")

mat <- left_join(mat, results)
hugo_gene <- mat$hugo_gene
mat <- cbind(hugo_gene, mat[2:27])

mat <- aggregate(. ~ hugo_gene, data = mat, sum)
mat <- filter(mat, hugo_gene != "")

write.table