library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
library(tidyr)
synapseLogin()
source("NCATS_helpers.R")
library(parallel)


this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/NCATS_drugs/listofDrugTargets.R"

##use polypharmacology db app list of mols
fp.evo <- readRDS("fpevo.rds")
evo <- readRDS("evotec_dgidb.rds")
evo <- evo %>% filter(N_quantitative > N_inactive | N_qualitative > N_inactive | N_DGIDB > 0)

##map with modifed helper functions from app
x<-synGet("syn11559906")@filePath
drugstruct<-read.table(x, sep = '\t', header = TRUE, comment.char = "", quote = "")

targs <- mclapply(as.character(drugstruct$smiles2), function(x){
  x<- getSimMols(x, 0.95)
  
  if(nrow(x)>0){
    x
  }else{
    as.data.frame(c("NothingFound"))
  }
  }, mc.cores = detectCores())

names(targs) <- as.character(drugstruct$smiles2) 
targs<-ldply(targs)

map<- read.table(synGet("syn11591280")@filePath, sep = "\t", quote = "", comment.char = "", header = T) 
colnames(map)[3] <- ".id"
targs2 <- left_join(targs, map) %>% filter(`Tanimoto Similarity` > 0.99) %>% dplyr::select(-6)

targs2 <- targs2 %>% left_join(evo, by = "Original_molecule_SMILES")
summary <- as.data.frame(table(targs2$Hugo_Gene)) %>% filter(Freq >= 5)
targs3 <- filter(targs2, Hugo_Gene %in% summary$Var1)

ncats <- read.table(synGet("syn8314523")@filePath, header = T, sep = "\t") %>% filter(Cell.line %in% c("HS01", "HS12", "Syn5", "Syn1"))


res<-lapply(unique(targs3$Hugo_Gene), function(x){
  foo <- targs3 %>% dplyr::filter(Hugo_Gene == x)syn8395266
  drugs.ncats <- ncats %>% filter(Sample.ID %in% foo$ncgc)
  null <- drugs.ncats %>% filter(NF2.status != "NF2 expressing")
  wt <- drugs.ncats %>% filter(NF2.status != "NF2 null")
  bar<-wilcox.test(wt$AUC, null$AUC)
  c(bar$statistic, "pval" = bar$p.value)
})

### nothing signficant yielded!!!

names(res) <- unique(targs3$Hugo_Gene)
res <- ldply(res)

targs4 <- dplyr::select(targs3, .id, Common_Name.x, Original_molecule_SMILES, `Tanimoto Similarity`, ncgc, Hugo_Gene) %>% 
  add_column("true" = rep(TRUE, nrow(targs3))) %>% distinct() %>% spread(Hugo_Gene, "true", fill = FALSE)

write.table(targs3, "NCATSTargetsGather.txt", sep ="\t", row.names = F)
synStore(File("NCATSTargetsGather.txt", parentId = "syn11609815"), used = c("syn11559906", "syn11591280", "https://github.com/Sage-Bionetworks/polypharmacology-db/blob/master/Data/fpevo.rds", "https://github.com/Sage-Bionetworkspolypharmacology-db/Data/fpevo.rds"),
         executed = this.file)

write.table(targs4, "NCATSTargetsSpread.txt", sep ="\t", row.names = F)
synStore(File("NCATSTargetsSpread.txt", parentId = "syn11609815"), used = c("syn11559906", "syn11591280", "https://github.com/Sage-Bionetworkspolypharmacology-db/Data/evotec_dgidb.RDS", "https://github.com/Sage-Bionetworkspolypharmacology-db/Data/fpevo.rds"),
         executed = this.file)

# 
# ## map uniprot targets to hugo genes - generate list of dataframes, each  
# ## data frame lists targets of drug
# targets <-
#   left_join(targets, bm, by = "Uniprot_accession_numbers")
# 
# targets.ncats<-filter(targets,  Structure_ID %in% alldata$Structure_ID)
# targets.ncats<-full_join(targets.ncats, alldata, "Structure_ID")
# targets.ncats<-targets.ncats %>% dplyr::distinct() %>% filter(!is.na(Uniprot_accession_numbers))
# targets.ncats<-dplyr::select(targets.ncats, Drug, Structure_ID, CID, CSID, Uniprot_accession_numbers, Hugo_Gene, Protein_names, 
#                       Min_activity_operator, MinActivity_nM, MeanActivity_nM, MedianActivity_nM, Geom_meanActivity_nM, 
#                       Max_activity_operator, MaxActivity_nM, N_quantitative, N_qualitative, N_inactive, Supplier , 
#                       Supplier_Molname, Supplier_Data_1, Supplier_Data_2, Supplier_Data_3)       
# targets.ncats.simple<-dplyr::select(targets.ncats, Drug, Structure_ID, CID, CSID, Uniprot_accession_numbers, Hugo_Gene, Protein_names)   
# 
# write.table(targets.ncats, "Full_NCATS_EVOTEC_Target_Map.txt", sep = "\t")
# synStore(File("Full_NCATS_EVOTEC_Target_Map.txt", parentId = "syn8398700"), used = c("syn8314523", "syn8361855", "syn7341038"), executed = this.file)
# write.table(targets.ncats.simple, "Simple_NCATS_EVOTEC_Target_Map.txt", sep = "\t")
# synStore(File("Simple_NCATS_EVOTEC_Target_Map.txt", parentId = "syn8398700"), used = c("syn8314523", "syn8361855", "syn7341038"), executed = this.file)
