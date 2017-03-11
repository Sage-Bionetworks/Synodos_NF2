library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
library(tidyr)
synapseLogin()

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/NCATS_drugs/listofDrugTargets.R"

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)

colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Compound.Name', 'AUC', 'AUC.Fit', 'Gene.Name')

#lookup pubchem IDs with https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
druggies<-as.data.frame(unique(drugdat$Compound.Name))
write.table(druggies, "NCATSdrugs.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)               
cids<-read.table("NCATS_drugs_to_CID.txt", sep ="\t", head = FALSE, quote = "")
names(cids) <- c("Drug", "CID")

#only one was not mapped to a CID, which is nice
length(intersect(cids$Drug, druggies$`unique(drugdat$Compound.Name)`))

##map drugs to EVOTEC ids using CID
evo<-synGet("syn8361855")@filePath
evodat<-read.table(evo)
evodat<-filter(evodat, !is.na(CID))
alldat<-left_join(evodat, cids, by = "CID")
alldat<-filter(alldat, !is.na(Drug))
alldata<-alldat[!duplicated(alldat$Structure_ID),]

##get targets using EVOTEC IDs and biomart
targets <- synTableQuery("SELECT * FROM syn7341038")
targets <- as.data.frame(targets@values)
targets <- filter(targets, Organism == 'Homo sapiens')

## obtain data to map Uniprot id with Hugo Genes
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"),
            mart = mart)
bm <- filter(bm, uniprot_swissprot != "")
colnames(bm) <- c("Hugo_Gene", "Uniprot_accession_numbers")

## map uniprot targets to hugo genes - generate list of dataframes, each  
## data frame lists targets of drug
targets <-
  left_join(targets, bm, by = "Uniprot_accession_numbers")

targets.ncats<-filter(targets,  Structure_ID %in% alldata$Structure_ID)
targets.ncats<-full_join(targets.ncats, alldata, "Structure_ID")
targets.ncats<-targets.ncats %>% dplyr::distinct() %>% filter(!is.na(Uniprot_accession_numbers))
targets.ncats<-dplyr::select(targets.ncats, Drug, Structure_ID, CID, CSID, Uniprot_accession_numbers, Hugo_Gene, Protein_names, 
                      Min_activity_operator, MinActivity_nM, MeanActivity_nM, MedianActivity_nM, Geom_meanActivity_nM, 
                      Max_activity_operator, MaxActivity_nM, N_quantitative, N_qualitative, N_inactive, Supplier , 
                      Supplier_Molname, Supplier_Data_1, Supplier_Data_2, Supplier_Data_3)       
targets.ncats.simple<-dplyr::select(targets.ncats, Drug, Structure_ID, CID, CSID, Uniprot_accession_numbers, Hugo_Gene, Protein_names)   

write.table(targets.ncats, "Full_NCATS_EVOTEC_Target_Map.txt", sep = "\t")
synStore(File("Full_NCATS_EVOTEC_Target_Map.txt", parentId = "syn8398700"), used = c("syn8314523", "syn8361855", "syn7341038"), executed = this.file)
write.table(targets.ncats.simple, "Simple_NCATS_EVOTEC_Target_Map.txt", sep = "\t")
synStore(File("Simple_NCATS_EVOTEC_Target_Map.txt", parentId = "syn8398700"), used = c("syn8314523", "syn8361855", "syn7341038"), executed = this.file)
