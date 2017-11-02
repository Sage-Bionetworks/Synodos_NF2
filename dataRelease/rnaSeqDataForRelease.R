library(synapseClient)
library(dplyr)
library(tidyr)
library(stringr)
synapseLogin()
this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/dataRelease/rnaSeqDataForRelease.R"

###human meningioma differential expression data
mendiff<-read.table(synGet("syn6166525")@filePath, header = T, sep = "\t") %>% 
  filter(cellLine1 != "HS01", cellLine1 != "HS11") %>% 
  filter(cellLine2 != "HS01", cellLine2 != "HS11")

write.table(mendiff, "meningioma_differential_expression.txt", sep = "\t")
synStore(File("meningioma_differential_expression.txt", parentId = "syn10881879"), used = "syn6166525", executed = this.file)
###

### human schwannoma differential expression data
schdiff<-read.table(synGet("syn9884855")@filePath, header = T, sep = "\t") %>%
  unite(geneName, ensembl, Hugo_Gene, sep = "|") %>%
  mutate(diffExptest = comparison) %>% 
  separate(comparison, c("cond1", "cond2"), sep = "vs") %>% 
  separate(cond1, c("cellLine1", "treatment1"), sep = 4) %>% 
  separate(cond2, c("cellLine2", "treatment2"), sep = 4) %>% 
  select(logFC, logCPM, `F`, PValue, BH, bonferroni, diffExptest, geneName, cellLine1, treatment1, cellLine2, treatment2)

write.table(schdiff, "schwannoma_differential_expression.txt", sep = "\t")
synStore(File("schwannoma_differential_expression.txt", parentId = "syn10881879"), used = "syn9884855", executed = this.file)
###

### mouse schwannoma data directly copied from internal project

###human meningioma raw counts matrix - remove non-Synodos samples and outdated schwannoma
mencount <- read.table(synGet("syn6045994")@filePath, header = T, sep = "\t") %>% select(
  gene,AC027_untreated,AC028_untreated,AC029_untreated,
  AC034_untreated,Syn4_untreated,Syn5_untreated,Syn2_untreated,
  Syn3_untreated,Syn6_untreated,
  Syn1.2_CUDC.907,Syn1.3_CUDC.907,Syn5.1_CUDC.907,Syn5.2_CUDC.907,
  Syn5.3_CUDC.907,Syn1.1_Panobinostat,Syn1.2_Panobinostat,Syn1.3_Panobinostat,
  Syn5.1_Panobinostat,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,
  Syn5.1_DMSO,Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,
  Syn6.2_DMSO,Syn1.1_CUDC.907,Syn5.1_GSK2126458,Syn5.2_GSK2126458,
  Syn5.3_GSK2126458,Syn6.1_GSK2126458,Syn6.2_GSK2126458,Syn6.3_DMSO,
  Syn6.3_CUDC.907,Syn6.3_Panobinostat,Syn6.3_GSK2126458,Syn6.1_CUDC.907,
  Syn6.2_CUDC.907,Syn5.2_Panobinostat,Syn5.3_Panobinostat,Syn6.1_Panobinostat,
  Syn6.2_Panobinostat,Syn1.1_GSK2126458,Syn1.2_GSK2126458,Syn1.3_GSK2126458,
  AC029.1_untreated,AC029.2_untreated,AC030.1_untreated,AC030.2_untreated,
  AC033.1_untreated,AC033.2_untreated,AC033.3_untreated,AC6.1_untreated,
  AC6.2_untreated,AC7.1_untreated,AC7.2_untreated,Syn10.1_DMSO,Syn10.2_DMSO,Syn10.3_CUDC.907,
  Syn10.4_CUDC.907,Syn10.5_Panobinostat,Syn10.6_Panobinostat,Syn10.7_GSK2126458,
  Syn10.8_GSK2126458,Syn2.1_untreated,Syn2.2_untreated,Syn2.3_untreated,
  Syn3.1_untreated,Syn3.2_untreated,Syn3.3_untreated,Syn4.1_untreated,
  Syn4.2_untreated,Syn4.3_untreated,Syn5.1_untreated)

write.table(mencount, "meningioma_raw_counts.txt", sep = "\t", row.names = F, quote = F)
synStore(File("meningioma_raw_counts.txt", parentId = "syn10881882"), used = "syn6045994", executed = this.file)

###human schwannoma raw counts matrix - renamed
schcount <- read.table(synGet("syn9779230")@filePath, header = T, sep = "\t") %>% rename("Gene" = X)
write.table(schcount, "schwannoma_raw_counts.txt", sep = "\t", row.names = F, quote = F)
synStore(File("schwannoma_raw_counts.txt", parentId = "syn10881882"), used = "syn9779230", executed = this.file)

### mouse schwannoma raw counts matrix - batch 5
mscount<-read.table(synGet("syn7467413")@filePath, header = T, sep = ",", comment.char = "") %>% 
  rename("chr" = X..chr) %>% 
  select(1:4,29:48)
write.table(mscount, "mouse_schwannoma_raw_counts_batch5.txt", sep = "\t", row.names = F, quote = F)
synStore(File("mouse_schwannoma_raw_counts_batch5.txt", parentId = "syn10881882"), used = "syn7467413", executed = this.file)

### mouse schwannoma raw counts matrix - batch 4
mscount2<-read.table(synGet("syn7467412")@filePath, header = T, sep = ",", comment.char = "") %>% rename("chr" = X..chr)
write.table(mscount2, "mouse_schwannoma_raw_counts_batch4.txt", sep = "\t", row.names = F, quote = F)
synStore(File("mouse_schwannoma_raw_counts_batch4.txt", parentId = "syn10881882"), used = "syn7467412", executed = this.file)

