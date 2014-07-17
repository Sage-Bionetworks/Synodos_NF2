source("http://bioconductor.org/biocLite.R")
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


library("GEOquery")
library("tools")
library("frma")
library("affy")
library("synapseClient")
synapseLogin()

#synid of the folder where the expression data will be stored
expressionData_synfolder_id <- "syn2481212"


#this script github URL
thisscript_URL <- "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/1_publicData_org.R"


##get the best probeSet per gene
HGU133Plus2.0 <- "http://www.cbs.dtu.dk/biotools/jetset/current/jetset.scores.hgu133plus2_2.14.0.csv.zip"
HGU133 <- "http://www.cbs.dtu.dk/biotools/jetset/current/jetset.scores.hgu133a_2.14.0.csv.zip"

get_probeset_to_geneMapping <- function(zipFile){
  temp <- tempfile()
  download.file(zipFile,temp)
  data <- read.table(unzip(temp),sep=",", header=T) 
  data  <- subset(data, best == "TRUE")
  #could be improved
  #currently using only the uniq gene symbols
  data <- data[!duplicated(data$symbol),] 
}

HGU133Plus2.0_probeSet_to_genes <- get_probeset_to_geneMapping(HGU133Plus2.0)
HGU133_probeSet_to_genes <- get_probeset_to_geneMapping(HGU133)
  
extractGEOdata <- function(f){
  if(file_ext(f) == "tar"){
    untar(f, exdir=dirname(f))
  }
}

get_GEO_ftpLink <- function(GEO){
  geotype <- toupper(substr(GEO, 1, 3))
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  if (geotype == "GSM") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/", 
                   stub, GEO)
  } else if (geotype == "GSE") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", 
                   stub, GEO)
  } else if (geotype == "GPL") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/", 
                   stub, GEO)
  } else {
    url <- None
  }
}

get_formatted_phenoData <- function(geo_study){
  df <- pData(geo_study)
  idx <- grep("characteristics", colnames(df), fixed=T, value=T)
  df <- df[,idx]
  for(i in idx){
    tmp <- strsplit(as.character(df[[i]]), ": ", fixed=T)
    nm <- sapply(tmp, "[[", 1)
    df[[unique(nm)]] <- sapply(tmp, "[[", 2)
    df[[i]] <- NULL
  }
  df['sampleId'] <- rownames(df)
  return(df)
}


######################
##MAIN
######################

#specify the ID of the GEO study to download and analyze
geo_id = 'GSE39645'

#create a folder 
study_folder <- synStore(Folder(name=geo_id, parentId = expressionData_synfolder_id ))
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
#Method 1
affyBatch_object <- ReadAffy(celfile.path=dirname(rownames(df)[1]))
#Method2
#create a vector of files names we want to read
#this is needed in the case where we dont want to read all the CEL files
#sometimes a GEO study has a mix of platform types and versions and each grp needs specific processing
cel_files <- paste0(base_dir,'/', geo_id, '/',rownames(pData(geo_study)), '.cel.gz')
affyBatch_object <- ReadAffy(filenames=cel_files)


#normalize the expression data
normData <- frma(affyBatch_object, summarize="random_effect")

#get the expression matrix
geneExp <- exprs(normData)


#filter the probeset ids to keep only the best probeset / gene 
#map the porbeset to gene symbol

#CHOOSE the right probeset based on version of the chip
probeSet_to_genes <- HGU133Plus2.0_probeSet_to_genes
probeSet_to_genes <- HGU133_probeSet_to_genes

geneExp <- geneExp[rownames(geneExp) %in% probeSet_to_genes$probeset, ]
temp_df <- data.frame('probeset' = rownames(geneExp))
temp_df <- merge(temp_df, probeSet_to_genes, by.x='probeset', by.y='probeset')

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


boxplot(log2(geneExp_raw), xlab="CEL File", ylab="expr", 
        main=paste0('GEO Study ', geo_id, ' expression boxplot' ))


dim(geneExp)
dim(phenotype_data)
phenotype_data
boxplot(geneExp, xlab="CEL File", ylab="expr", 
        col = factor(phenotype_data['who grade']),
        main=paste0('GEO Study ', geo_id, ' expression boxplot' ))




