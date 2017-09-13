library(biomaRt)

ensembl = useEnsembl(biomart="ensembl",GRCh=37,dataset = "hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

dat <- read.table("HS11panobinostatvsHS11DMSO.txt", sep = "\t", header = TRUE, quote = "", comment.char = "")
dat$Genes <- as.character(dat$Genes)

for(i in 1:nrow(dat)){
  genelist <- as.character(dat[[i,4]])
  genelist <- (unlist(strsplit(genelist, split=", ")))
  genelist <- gsub("\"", "", genelist)
  print(length(genelist))
  hugo <- getBM(mart = ensembl, filters = "ensembl_gene_id", attributes = c("hgnc_symbol"), values = genelist)
  print(nrow(hugo))
  dat[[i,4]] <- paste(hugo$hgnc_symbol, sep = "", collapse = ", ")
}

write.table(dat, "HS11panobinostatvsHS11DMSO_2.txt", sep = "\t")
