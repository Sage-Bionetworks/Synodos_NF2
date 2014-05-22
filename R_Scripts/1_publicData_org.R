#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("frma")
biocLite("hugene10stv1frmavecs")
biocLite("oligo")
biocLite("hugene10stv1.r3cdf")
biocLite("hgu133plus2cdf")
biocLite("hgu133plus2frmavecs")


library("GEOquery")
library("tools")
library("frma")
library("affy")
library("synapseClient")
synapseLogin()

#
expressionData_synfolder_id <- "syn2481212"


#this script github URL
thisscript_URL <- 

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
  return(df)
}


######################
##MAIN
######################

#specify the ID of the GEO study to download and analyze
geo_id = 'GSE16581'

#create a folder 
study_folder <- synStore(Folder(name=geo_id, parentId = expressionData_synfolder_id ))
#push the study page as a reference
source_url = get_GEO_ftpLink(geo_id)
source_url = synStore(File(source_url, name=geo_id, synapseStore=F, parentId=study_folder$properties$id))

orig_working_dir <- getwd()
base_dir <- "/external-data/DAT_123__Synodos_NF2/publicData//expressionData/"
setwd(base_dir)

#get the phenotype data
geo_study <- getGEO(geo_id)[[1]]
phenotype_data <- get_formatted_phenoData(geo_study)
phenotype_data_file <- paste0(base_dir, geo_id, '/', geo_id, '_phenotype_data.txt')
write.table(phenotype_data, file=phenotype_data_file, col.names=T, row.names=T, quote=FALSE)
phenotype_data_file <- synStore(File(phenotype_data_file, parentId =study_folder$properties$id ), used= source_url$properties$id)

#download the raw data
df <- getGEOSuppFiles(geo_id)
setwd(orig_working_dir)
sapply(rownames(df),extractGEOdata) 

#process the expression data
affyBatch_object <- ReadAffy(celfile.path=dirname(rownames(df)[1]))
normData <- frma(affyBatch_object, summarize="random_effect")
geneExp <- exprs(normData)


#boxplot
boxplot_file <- paste0(base_dir, geo_id, '/', geo_id, '_expBoxplot.png')
png(boxplot_file)
boxplot(geneExp, col=factor(phenotype_data$gender), xlab="CEL File", ylab="expr")
?boxplot
dev.off()
boxplot_file <- synStore(File(boxp))






