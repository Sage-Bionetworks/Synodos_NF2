library(synapseClient)
library(dplyr)
library(ggplot2)
library(biomaRt)
synapseLogin()

files<- c("syn10149538",
             "syn10149539",
             "syn10149540",
             "syn10149542",
             "syn10149543",
             "syn10149544",
             "syn10149545",
             "syn10149546")
             

dat<-lapply(files, function(x){
  foo <- synGet(x)
  nam <- gsub("_R1_001.uniq.sorted.bam.txt", "", foo@fileHandle$fileName)
  bar <- read.table(foo@filePath, header = FALSE)
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
mat <- cbind(Hugo_Gene, mat[2:9])

mat <- filter(mat, Hugo_Gene != "")
mat <- aggregate(. ~ Hugo_Gene, data = mat, sum)

write.table(mat, "schwannoma_reseq_raw_counts_bedtoolsHtseq.txt", sep = "\t")

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCExpressionMatrix-bedtoolsHtseq.R"
synStore(File("schwannoma_reseq_raw_counts_bedtoolsHtseq.txt", parentId = "syn9925306"), used = files, executed = this.file)


