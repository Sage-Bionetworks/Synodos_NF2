library(synapseClient)
synapseLogin()
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)

get_peptide_data <- function(id, filter=TRUE){
  peptide_file <- synGet(id)
  peptides <- fread(peptide_file@filePath, header=T, data.table=F)
  cat(paste("processing", id, "\n"))
  #remove unwanted cols
  unwanted_cols <- c('Search ID', 'Processing Node No', 'Sequence', 'Unique Sequence ID',
                     'Protein Descriptions' ,'Modifications', 'Activation Type', 
                     'DeltaScore', 'DeltaCn', 'Rank', 'Search Engine Rank',
                     'Precursor Area', 'QuanResultID', 'Decoy Peptides Matched',
                     'Peptides Matched', 'XCorr', '# Missed Cleavages',
                     'Ion Inject Time [ms]', 'Charge',  'm/z [Da]',  'MH+ [Da]',  'Delta Mass [Da]',
                     'Delta Mass [PPM]', 'RT [min]', 'First Scan', 'Last Scan',  'MS Order',  'Ions Matched',
                     'Matched Ions', 'Total Ions',  'Spectrum File', 'Annotation', 'Confidence_Level',
                     'PSM Ambiguity', 'q-Value', 'PEP', 'Intensity',
                     'Confidence Level', '# Proteins', '# Protein Groups')
  
  peptides <- peptides[,!colnames(peptides) %in% unwanted_cols]
  colnames(peptides) <-  gsub('|\\[|\\]|\\%|', '',colnames(peptides))
  colnames(peptides) <-  gsub('\\s+$', '',colnames(peptides), perl=T)
  colnames(peptides) <-  gsub('# ', '',colnames(peptides))
  colnames(peptides) <-  gsub(' ', '_',colnames(peptides))
  
  if (filter == TRUE){
    peptides <- peptides %>% 
      filter( Quan_Info == 'Unique' & Quan_Usage == 'Used' & Isolation_Interference <= 75 )
    #remove further cols
    peptides <- peptides[,!colnames(peptides)  %in% c('Quan_Info', 'Quan_Usage', 'Isolation_Interference')]
  }
  peptides['fileId'] = id
  peptides
}

######################
# Meningioma Data
meningioma_peptide_ids = synapseQuery('select id,name from file where parentId == "syn4214447"')
meningioma_peptide_ids <- meningioma_peptide_ids[grepl('peptides', meningioma_peptide_ids$file.name),]
meningioma_peptide_data <- lapply(meningioma_peptide_ids$file.id, get_peptide_data)
#melting
meningioma_peptide_data <- lapply(meningioma_peptide_data, function(x){
  melt(x, id.vars=c('Protein_Group_Accessions', 'fileId'))
})
#combining
meningioma_peptide_data <- do.call(rbind, meningioma_peptide_data)

##### Schwannoma
schwannoma_peptide_ids = synapseQuery('select id,name from file where parentId == "syn4214452"')
schwannoma_peptide_ids <- schwannoma_peptide_ids[grepl('peptides', schwannoma_peptide_ids$file.name),]
schwannoma_peptide_data <- lapply(schwannoma_peptide_ids$file.id, get_peptide_data)
#melting
schwannoma_peptide_data <- lapply(schwannoma_peptide_data, function(x){
  melt(x, id.vars=c('Protein_Group_Accessions', 'fileId'))
})

#combining
schwannoma_peptide_data <- do.call(rbind, schwannoma_peptide_data)
kinome_peptide_data <- rbind(schwannoma_peptide_data, meningioma_peptide_data)
colnames(kinome_peptide_data) <- c('Accession', 'fileId', 'condition','log2ratio')
kinome_peptide_data$log2ratio <- log2(kinome_peptide_data$log2ratio)
#split the meta data col
kinome_peptide_data <- separate(data=kinome_peptide_data, col=condition, 
                                into=c('cellLine', 'replicate','drug', 'time'),
                                sep='_')

# Upload the data to synapse
parentId = 'syn4259360'
outfile <- "Synodos_kinome_peptide_data.csv"
write.table(kinome_peptide_data, file=outfile, sep="\t")
used_ids <- c(schwannoma_peptide_ids$file.id, meningioma_peptide_ids$file.id)
tmp <- synStore(File(outfile, parentId = parentId), used=used_ids)
unlink(outfile)
