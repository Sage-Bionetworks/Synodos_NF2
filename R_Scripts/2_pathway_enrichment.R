library("synapseClient")
library("GSVA")

synapseLogin()

#get the MsigDB object
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
cat('..Done\n\n')



studyName <- 'GSE16581'
upload_bucket <- 'syn2482470'
geneExp_synId <- 'syn2504985'
geneExp <- synGet(geneExp_synId)
geneExp <- read.table(geneExp@filePath, sep="\t", header=T)


res_ssGSEA <- gsva(as.matrix(geneExp), MSigDB$C2.CP.KEGG, method="ssgsea")
res_gsva <- gsva(as.matrix(geneExp), MSigDB$C2.CP.KEGG)


#get the phenotype data
geo_study <- getGEO(geo_id)
geo_study <- geo_study[[1]]
phenotype_data <- get_formatted_phenoData(geo_study)
phenotype_data_file <- paste0(base_dir, geo_id, '/', geo_id, '_phenotype_data.tsv')
write.table(phenotype_data, phenotype_data_file, col.names=T, row.names=F, quote=FALSE, sep="\t")
phenotype_data_file <- synStore(File(phenotype_data_file, parentId =study_folder$properties$id ), used= source_url$properties$id)



tempdir()
(res_ssGSEA)

?gsva

heatmap(res_ssGSEA,scale="none",center="none")
heatmap(res_gsva$es.obs)

head(res_ssGSEA)