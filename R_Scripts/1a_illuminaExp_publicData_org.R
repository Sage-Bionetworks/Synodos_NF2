source("http://bioconductor.org/biocLite.R")
biocLite (c("beadarray", "limma", "GEOquery", "illuminaHumanv1.db",
             "illuminaHumanv2.db", "illuminaHumanv3.db", "BeadArrayUseCases" ))
biocLite("illuminaHumanv4.db")
biocLite (c( "GOstats", "GenomicRanges", "Biostrings" ))



library("GEOquery")
library("tools")
library("frma")
library("affy")
library("synapseClient")
synapseLogin()

#
expressionData_synfolder_id <- "syn2481212"

#this script github URL
thisscript_URL <- "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/1_publicData_org.R"

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
geo_id = 'GSE58037'

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
df
setwd(orig_working_dir)
rownames(df)
sapply(rownames(df),extractGEOdata) 


#get the phenotype data
geo_study <- getGEO(geo_id)
geo_study <- geo_study[[1]]
phenotype_data <- get_formatted_phenoData(geo_study)
phenotype_data_file <- paste0(base_dir, geo_id, '/', geo_id, '_phenotype_data.tsv')
write.table(phenotype_data, phenotype_data_file, col.names=T, row.names=F, quote=FALSE, sep="\t")
phenotype_data_file <- synStore(File(phenotype_data_file, parentId =study_folder$properties$id ), used= source_url$properties$id)

exprs(geo_study)

#process the expression data
affyBatch_object <- ReadAffy(celfile.path=dirname(rownames(df)[1]))
#num probes
#dim(pm(affyBatch_object))
normData <- frma(affyBatch_object, summarize="random_effect")

geneExp_raw <- exprs(affyBatch_object)
geneExp <- exprs(normData)

exp_file <- paste0(base_dir, geo_id, '/', geo_id, '_expvals.txt')
write.table(geneExp, file=exp_file, col.names=T, row.names=T, quote=F)
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
        col = factor(phenotype_data$tissue),
        main=paste0('GEO Study ', geo_id, ' expression boxplot' ))




