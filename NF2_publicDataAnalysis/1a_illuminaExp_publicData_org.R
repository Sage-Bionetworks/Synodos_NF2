code_base_dir <- paste0(path.expand("~"),"/dev/Synodos_NF2/R_Scripts/" )

#source("http://bioconductor.org/biocLite.R")
#biocLite (c("beadarray", "limma", "GEOquery", "illuminaHumanv1.db",
#             "illuminaHumanv2.db", "illuminaHumanv3.db", "BeadArrayUseCases" ))
#biocLite("illuminaHumanv4.db")
#biocLite (c( "GOstats", "GenomicRanges", "Biostrings" ))


library("GEOquery")
library("tools")
library("synapseClient")
library("beadarray")
library("illuminaHumanv4.db")
library("sva")
synapseLogin()

#call the GEO functions
source(paste0(code_base_dir,"/lib/GEO_functions.R"))

#synid of the folder where the expression data will be stored
meningioma_studies_bucket <- "syn2537450"
schwannoma_studies_bucket <- "syn2537451"
synodos_code_dir <- "syn2572412"


#this script github URL
thisscript_URL <- "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/1a_illuminaExp_publicData_org.R"


#####################
##MAIN
######################

#specify the ID of the GEO study to download and analyze
geo_id = 'GSE58037'

#create a folder 
study_folder <- synStore(Folder(name=geo_id, parentId = meningioma_studies_bucket ))
#push the study page as a reference
source_url = get_GEO_ftpLink(geo_id)
source_url = synStore(File(source_url, name=geo_id, synapseStore=F, parentId=study_folder$properties$id))

base_dir <- "/external-data/DAT_123__Synodos_NF2/publicData//expressionData/"
setwd(base_dir)
#download the raw data
#df <- getGEOSuppFiles(geo_id)
# setwd(orig_working_dir)
# getwd()
# rownames(df)
# sapply(rownames(df),extractGEOdata) 


#get the phenotype data
geo_study <- getGEO(geo_id)
geo_study <- geo_study[[1]]
phenotype_data <- get_formatted_phenoData(geo_study)
phenotype_data_file <- paste0(base_dir, geo_id, '/', geo_id, '_phenotype_data.tsv')
write.table(phenotype_data, phenotype_data_file, col.names=T, row.names=F, quote=FALSE, sep="\t")
phenotype_data_file <- synStore(File(phenotype_data_file, parentId =study_folder$properties$id ), used= source_url$properties$id)



geo_study <- addFeatureData(geo_study, 
                            toAdd = c("SYMBOL", "PROBEQUALITY", "CODINGZONE", "PROBESEQUENCE", "GENOMICLOCATION"),
                            annotation="Humanv4")

#get the feature names of the illumina array
ids <- as.character(featureNames(geo_study))
probe_quality <- unlist(mget(ids, illuminaHumanv4PROBEQUALITY, ifnotfound=NA))

#remove low quality probes
to_remove <- probe_quality == "No match" | probe_quality == "Bad" | is.na(probe_quality)
#remove the bad probes from the downstream analysis
geo_study.flt <- geo_study[!to_remove,]

#create  a matrix of gene expression
geneExp <- exprs(geo_study.flt)

#convert to illumina probe id's to geneNames(symbols)
geneNames <- unlist(mget(rownames(geneExp),illuminaHumanv4SYMBOL, ifnotfound=NA))
rownames(geneExp) <- geneNames

#remove the rows where the genename is NA
rows_to_keep <- !is.na(rownames(geneExp))
geneExp <- geneExp[rows_to_keep,]

#take avg of the rows for the same geneNames
#multiple probes mapping to a same gene are summarized (mean)
avg_genes <- function(gene){
  e <- geneExp[rownames(geneExp) %in% gene,]
  apply(e,2,mean)
}

duplicate_genes <- names(which(table(rownames(geneExp)) > 1))
avg_geneExp_duplicate_genes <- sapply(duplicate_genes,avg_genes)
avg_geneExp_duplicate_genes <- t(avg_geneExp_duplicate_genes)

#keep only the genes which have one unique probe
geneExp <- geneExp[(!rownames(geneExp) %in% duplicate_genes),]

#now add the avg expression values for genes with multiprobes
geneExp <- rbind(geneExp, avg_geneExp_duplicate_genes)


exp_file <- paste0(base_dir, geo_id, '/', geo_id, '_expvals.txt')
write.table(geneExp, file=exp_file, col.names=T, row.names=T, quote=F, sep="\t")
exp_file <- synStore(File(exp_file, parentId=study_folder$properties$id),
                     used=source_url$properties$id,
                     executed=thisscript_URL)

#boxplot
boxplot_file <- paste0(base_dir, geo_id, '/', geo_id, '_expBoxplot.png')
png(boxplot_file, width=8, height=6, units="in", res=200)
boxplot(geneExp, xlab="sample", ylab="expr", main=paste0('GEO Study ', geo_id, ' expression boxplot' ))
dev.off()
boxplot_file <- synStore(File(boxplot_file, parentId=study_folder$properties$id), 
                         used=source_url$properties$id,
                         executed=thisscript_URL)


#####
# Data Exploration Code
#####
boxplot(geneExp, xlab="CEL File", ylab="expr", main=paste0('GEO Study ', geo_id, ' expression boxplot' ))

#top 1000 genes that vary
genes_var <- apply(geneExp,1,var)
top_genes <- geneExp[order(-genes_var),][1:500,]
colnames(top_genes) <- gsub('.cel.gz','',colnames(top_genes))
library(limma)

r1 <- eBayes(lmFit(geneExp, model.matrix(~ phenotype_data$batch + phenotype_data$center)))
hist(r1$p.value[,-1])

top_genes_scaled <- t(scale(t(top_genes)))
memoised_pheatmap(top_genes_scaled, annotation=phenotype_data[,c('center','batch')])


batch <- phenotype_data[,c('batch', 'center'),drop=F]
mod.matrix <- model.matrix(~1, data=phenotype_data)  
mod.matrix
x <- ComBat(geneExp, batch=batch)
r2 <- eBayes(lmFit(x, model.matrix(~ phenotype_data$batch + phenotype_data$center)))
hist(r2$p.value[,-1])

?ComBat

top <- x[rownames(x) %in% rownames(top_genes),]
top <- t(scale(t(top)))
memoised_pheatmap(top, annotation=phenotype_data[,c('center','batch')])


boxplot(x)

?model.matrix