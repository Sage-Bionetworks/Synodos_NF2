library(synapseClient)
synapseLogin()
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)

source("kinomeAnalysis/kinome_data_functions.R")

###### Meningioma Data
meningioma_peptide_ids = synapseQuery('select id,name from file where parentId == "syn5870242"')
meningioma_peptide_ids <- meningioma_peptide_ids[grepl('peptides', meningioma_peptide_ids$file.name),]
meningioma_peptide_data <- lapply(meningioma_peptide_ids$file.id, get_peptide_data)
#melting
meningioma_peptide_data <- lapply(meningioma_peptide_data, function(x){
  melt(x, id.vars=c('Protein_Group_Accessions', 'fileId'))
})
#combining
meningioma_peptide_data <- do.call(rbind, meningioma_peptide_data)

##### Schwannoma
schwannoma_peptide_ids = synapseQuery('select id,name from file where parentId == "syn5870281"')
schwannoma_peptide_ids <- schwannoma_peptide_ids[grepl('peptides', schwannoma_peptide_ids$file.name),]
schwannoma_peptide_data <- lapply(schwannoma_peptide_ids$file.id, get_peptide_data)
#melting
schwannoma_peptide_data <- lapply(schwannoma_peptide_data, function(x){
  melt(x, id.vars=c('Protein_Group_Accessions', 'fileId'))
})
schwannoma_peptide_data <- do.call(rbind, schwannoma_peptide_data)


#combining
kinome_peptide_data <- rbind(schwannoma_peptide_data, meningioma_peptide_data)
colnames(kinome_peptide_data) <- c('Accession', 'fileId', 'condition','log2ratio')
kinome_peptide_data$log2ratio <- log2(kinome_peptide_data$log2ratio)
#split the meta data col
kinome_peptide_data <- separate(data=kinome_peptide_data, col=condition, 
                                into=c('cellLine', 'replicate','drug', 'time'),
                                sep='_')


#normalizing by median kinase ratio from protein
protein_derived_normFactors <- read.csv(synGet('syn4951080')@filePath, header=T, sep="\t")
protein_derived_normFactors <- protein_derived_normFactors %>% select(cellLine, drug, replicate, time, median_kinase_ratio)
protein_derived_normFactors <- protein_derived_normFactors[!duplicated(protein_derived_normFactors),]
kinome_peptide_data <- merge(kinome_peptide_data, protein_derived_normFactors)
kinome_peptide_data['log2NormRatio'] <- kinome_peptide_data$log2ratio - log2(kinome_peptide_data$median_kinase_ratio)

# Upload the data to synapse
parentId = 'syn4259360'
outfile <- "Synodos_kinome_data_peptide_level.csv"
write.table(kinome_peptide_data, file=outfile, sep="\t")
used_ids <- c(schwannoma_peptide_ids$file.id, meningioma_peptide_ids$file.id)
tmp <- synStore(File(outfile, parentId = parentId), used=used_ids)
unlink(outfile)
