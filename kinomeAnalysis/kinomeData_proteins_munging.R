library(synapseClient)
library("gdata")
library("tidyr")
library("dplyr")
library("reshape2")
require("parallel")
library("plyr")
library("doMC")
registerDoMC(8)
synapseLogin()

source("kinomeAnalysis/kinome_data_functions.R")

#OLD format data
# #process KINOME RUN 1 : FullSerum
# kinomeRun1 <- synGet('syn2679231')
# kinomeRun1 <- read.xls(kinomeRun1@filePath, sheet=2, header=T, check.names=F)
# #keep only selected columns (75% co-isolation)
# kinomeRun1 <- kinomeRun1[,c(1:14)]
# kinomeRun1_iTRAQ_to_cellLine <- list('116'='Syn2','117'='Syn3','114'='Syn5_1','115'='Syn5_2')                                                    
# kinomeRun1_condition <- 'FullSerum'
# kinomeRun1 <- get_old_kinome_proteinData(kinomeRun1, kinomeRun1_iTRAQ_to_cellLine, kinomeRun1_condition)
# 
# #process KINOME RUN 2 : Serum Free
# kinomeRun2 <- synGet('syn2679230')
# kinomeRun2 <- read.xls(kinomeRun2@filePath, sheet=2, header=T, check.names=F)
# #keep only selected columns (75% co-isolation)
# kinomeRun2 <- kinomeRun2[,c(1:14)]
# kinomeRun2_iTRAQ_to_cellLine <- list('114'='Syn2','115'='Syn3','116'='Syn5_1','117'='Syn5_2')                                                    
# kinomeRun2_condition <- 'SerumFree'
# kinomeRun2 <- get_old_kinome_proteinData(kinomeRun2, kinomeRun2_iTRAQ_to_cellLine, kinomeRun2_condition)


###############
# New format data (post April 12, 2016 - this is the non-Normalized data)
# we will do custom normalization - median kinase 
##############
#get all Single drug treatment kinome data 
kinome_data_ids = synapseQuery('select id,name from file where parentId == "syn6178256"')
kinome_data_ids <- kinome_data_ids[!grepl('peptides', kinome_data_ids$file.name),]
kinome_data <- lapply(kinome_data_ids$file.id, get_new_kinome_proteinData)
kinome_data <- do.call(rbind, kinome_data)
kinome_data <- kinome_data %>% separate(sample, into=c('cellLine', 'replicate', 'drug', 'time'), sep="_")


##########
# create normalization factor
##########
normFactors <- kinomeData %>% dplyr::group_by(cellLine, drug, replicate, time) %>% dplyr::summarise(median_kinase_ratio=median(ratio))
kinomeData <- merge(normFactors, kinomeData)


kinomeData['log2ratio'] = log2(as.numeric(kinomeData$ratio))
kinomeData['log2NormRatio'] = log2(kinomeData$ratio) - log2(kinomeData$median_kinase_ratio)
kinomeData['uniq_peptides'] =  kinomeData['# Unique Peptides']
kinomeData['# Unique Peptides'] <- NULL

#creating a unique col to identify each condition of run
kinomeData['condition'] = apply(kinomeData[,c('cellLine', 'drug', 'replicate', 'time')], 1,function(x) paste(x, collapse='_') )

#fix colnames
cnames <- gsub('# ','num_',colnames(kinomeData))
cnames <- gsub(' ','',cnames)
colnames(kinomeData) <- cnames

# Upload the data to synapse
parentId = 'syn4259360'
outfile <- "Synodos_kinome_data_protein_level.csv"
write.table(kinomeData, file=outfile, sep="\t")
used_ids <- c(meningioma_ids$file.id,schwannoma_ids$file.id)
tmp <- synStore(File(outfile, parentId = parentId),
                used=used_ids)
unlink(outfile)
