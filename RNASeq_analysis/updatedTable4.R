library(dplyr)
library(tidyr)

dat<-read.table("Supplemental_Table_5.txt", sep = "\t", header = T, skip = 1L)[,1:4]
colnames(dat) <- c("GeneID", "logFC.m", "logCPM.m", "BH.m")

sch<-read.table(synGet("syn9884491")@filePath, sep = "\t", header = T) %>% 
  unite(GeneID, gene, ensembl, sep = "|") %>% 
  select(GeneID, logFC, logCPM, BH) 

dat <- full_join(dat, sch) %>% 
  arrange(BH.m)
 
write.table(dat, "updatedtable4.txt", sep = "\t")
