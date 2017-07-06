library(synapseClient)
library(dplyr)
library(ggplot2)
library(biomaRt)
synapseLogin()

files<- c("syn10140183",
             "syn10140184",
             "syn10140185",
             "syn10140186",
             "syn10140188",
             "syn10140189",
             "syn10140190",
             "syn10140192",
             "syn10140193",
             "syn10140195",
             "syn10140196",
             "syn10140197",
             "syn10140199",
             "syn10140200",
             "syn10140202",
             "syn10140203",
             "syn10140205",
             "syn10140206",
             "syn10140207",
             "syn10140208",
             "syn10140209",
             "syn10140210",
             "syn10140211",
             "syn10140212",
             "syn10140213",
             "syn10140214")

dat<-lapply(files, function(x){
  foo <- synGet(x)
  nam <- gsub("_R1_001.uniq.sorted.txt", "", foo@fileHandle$fileName)
  bar <- read.table(foo@filePath, header = FALSE) %>% dplyr::select(V4, V13)
  colnames(bar) <- c("ensembl", nam)
  bar
})

mat <- dat %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ensembl"), .)

mat <- aggregate(. ~ ensembl, data = mat, sum)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 mart = mart)
colnames(results) <- c("ensembl", "Hugo_Gene")

mat <- left_join(mat, results)

#map <- distinct(dplyr::select(mat, Hugo_Gene, ensembl)) %>% 
#  dplyr::filter(Hugo_Gene != "" & Hugo_Gene != "NA") %>% 
#  group_by(Hugo_Gene) %>% 
#  top_n(1, ensembl) %>% 
#  ungroup()

#write.table(map, "ensemble_hugomap.txt", sep = "\t", quote = F)

Hugo_Gene <- mat$Hugo_Gene
mat <- cbind(Hugo_Gene, mat[2:27])

mat <- filter(mat, Hugo_Gene != "")
mat <- aggregate(. ~ Hugo_Gene, data = mat, sum)

write.table(mat, "schwannoma_reseq_raw_counts.txt", sep = "\t")

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCExpressionMatrix.R"
synStore(File("schwannoma_reseq_raw_counts.txt", parentId = "syn9925306"), used = files, executed = this.file)


