code_base_dir <- paste0(path.expand("~"),"/dev/Synodos_NF2/R_Scripts/" )

#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
# biocLite("frma")
# biocLite("hugene10stv1frmavecs")
# biocLite("oligo")
# biocLite("hugene10stv1.r3cdf")
# biocLite("hgu133plus2cdf")
# biocLite("hgu133plus2frmavecs")
# biocLite("hgu133afrmavecs")
# biocLite("hugene10stv1cdf")
# biocLite("hgu133plus2frmavecs")
#biocLite("hgu133afrmavecs")
#source("http://bioconductor.org/biocLite.R")
#biocLite("hugene.1.0.st.v1frmavecs")
#biocLite("hugene10stprobeset.db")


#call the GEO functions
source(paste0(code_base_dir,"/lib/GEO_functions.R"))


library("GEOquery")
library("tools")
library("frma")
library("affy")
library("synapseClient")

#login to synapse
synapseLogin()

#synid of the folder where the expression data will be stored
meningioma_studies_bucket <- "syn2537450"
schwannoma_studies_bucket <- "syn2537451"


#this script github URL
thisscript_URL <- "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/1_publicData_org.R"

##get the best probeSet per gene
HGU133Plus2.0 <- "http://www.cbs.dtu.dk/biotools/jetset/current/jetset.scores.hgu133plus2_2.14.0.csv.zip"
HGU133A <- "http://www.cbs.dtu.dk/biotools/jetset/current/jetset.scores.hgu133a_2.14.0.csv.zip"
U133_X3P <- "http://www.cbs.dtu.dk/biotools/jetset/current/jetset.scores.u133x3p_2.14.0.csv.zip"

get_probeset_to_geneMapping <- function(zipFile){
  temp <- tempfile()
  download.file(zipFile,temp)
  data <- read.table(unzip(temp),sep=",", header=T) 
  data  <- subset(data, best == "TRUE")
  #could be improved
  #currently using only the uniq gene symbols
  data <- data[!duplicated(data$symbol),] 
}

get_probeset_to_geneMapping_hugene10st <- function(){
  library("hugene10sttranscriptcluster.db")
  df <- toTable(hugene10sttranscriptclusterSYMBOL)
  df$probeset <- df$probe_id  
  df$probe_id <- NULL
  df
}

get_probeset_to_geneMapping_hgu133a <- function(){
  library("hgu133a.db")
  df <- toTable(hgu133aSYMBOL)
  df$probeset <- df$probe_id  
  df$probe_id <- NULL  
  df
}


#function
#to map the probeid to genenames AND
#avg the geneExp values for probes mapping to same gene 
process_geneExp <- function(geneExp, probeSet_to_genes ){
  dim(geneExp)
  geneExp <- geneExp[rownames(geneExp) %in% probeSet_to_genes$probeset, ]
  temp_df <- data.frame('probeset' = rownames(geneExp))
  temp_df <- merge(temp_df, probeSet_to_genes, by.x='probeset', by.y='probeset')
  rownames(geneExp) <- temp_df$symbol
  
  #take avg of the rows for the same geneNames
  #multiple probes mapping to a same gene are summarized (mean)
  avg_genes <- function(gene){
    e <- geneExp[rownames(geneExp) %in% gene,]
    apply(e,2,mean)
  }
  
  duplicate_genes <- names(which(table(rownames(geneExp)) > 1))
  if(length(duplicate_genes) >= 1){
    avg_geneExp_duplicate_genes <- sapply(duplicate_genes,avg_genes)
    avg_geneExp_duplicate_genes <- t(avg_geneExp_duplicate_genes)
    #keep only the genes which have one unique probe
    geneExp <- geneExp[(!rownames(geneExp) %in% duplicate_genes),]
    #now add the avg expression values for genes with multiprobes
    geneExp <- rbind(geneExp, avg_geneExp_duplicate_genes)  
  }
  geneExp
}

HGU133Plus2.0_probeSet_to_genes <- get_probeset_to_geneMapping(HGU133Plus2.0)
HGU133A_probeSet_to_genes <- get_probeset_to_geneMapping(HGU133)
Aff_HUGENE10_ST_probeSet_to_genes <- get_probeset_to_geneMapping_hugene10st()
Affy_U133_XP <- get_probeset_to_geneMapping(U133_X3P)
Affy_HGU133A <- get_probeset_to_geneMapping_hgu133a()


######################
##MAIN
######################
#specify the ID of the GEO study to download and analyze
geo_id = 'GSE43290'

#create a folder 
## choose the parentId based on study type (meningioma & schwannoma)
parentId <- meningioma_studies_bucket
study_folder <- synStore(Folder(name=geo_id, parentId = parentId ))
#push the study page as a reference
source_url = get_GEO_ftpLink(geo_id)
source_url = synStore(File(source_url, name=geo_id, synapseStore=F, parentId=study_folder$properties$id))

orig_working_dir <- getwd()
base_dir <- "/external-data/DAT_123__Synodos_NF2/publicData//expressionData/"
setwd(base_dir)

#download the raw data
df <- getGEOSuppFiles(geo_id)
setwd(orig_working_dir)
rownames(df)
sapply(rownames(df),extractGEOdata) 



#get the phenotype data
geo_study <- getGEO(geo_id)
geo_study <- geo_study[[1]]



phenotype_data <- get_formatted_phenoData(geo_study)
phenotype_data_file <- paste0(base_dir, geo_id, '/', geo_id, '_phenotype_data.tsv')
write.table(phenotype_data, phenotype_data_file, col.names=T, row.names=F, quote=FALSE, sep="\t")
phenotype_data_file <- synStore(File(phenotype_data_file, parentId =study_folder$properties$id ), 
                                used= source_url$properties$id)



#read the expression data
#choose the right method based on a method below

#Method 1
affyBatch_object <- ReadAffy(celfile.path=dirname(rownames(df)[1]))

#Method2
#create a vector of files names we want to read
#this is needed in the case where we dont want to read all the CEL files
#sometimes a GEO study has a mix of platform types and versions and each grp needs specific processing
#cel_files <- paste0(base_dir,'/', geo_id, '/',rownames(pData(geo_study)), '.cel.gz')
#affyBatch_object <- ReadAffy(filenames=cel_files)

#Method 3
#reading Affy 1.0 data
#read the pre-normalized data
geneExp <- exprs(geo_study)


#optional if normalization needed
#normalize the expression data
normData <- frma(affyBatch_object, summarize="random_effect")


#get the expression matrix
geneExp <- exprs(normData)

#or 
#if directly using data from study
#geneExp <- exprs(geo_study)

#filter the probeset ids to keep only the best probeset / gene 
#map the porbeset to gene symbol
#CHOOSE the right probeset based on version of the chip

# HGU133Plus2.0_probeSet_to_genes OR HGU133_probeSet_to_genes
# OR Aff_HUGENE10_ST_probeSet_to_genes

geneExp <- process_geneExp(geneExp, HGU133A_probeSet_to_genes )


exp_file <- paste0(base_dir, geo_id, '/', geo_id, '_expvals.tsv')
write.table(geneExp, file=exp_file, col.names=T, row.names=T, quote=F, sep="\t")
exp_file <- synStore(File(exp_file, parentId=study_folder$properties$id),
                     used=source_url$properties$id,
                     executed=thisscript_URL)

#boxplot
boxplot_file <- paste0(base_dir, geo_id, '/', geo_id, '_expBoxplot.png')
png(boxplot_file, width=8, height=6, units="in", res=200)
boxplot(geneExp, xlab="CEL File", ylab="expr", main=paste0('GEO Study ', geo_id, ' expression boxplot' ))
dev.off()
boxplot_file <- synStore(File(boxplot_file, parentId=study_folder$properties$id), 
                         used=source_url$properties$id,
                         executed=thisscript_URL)


