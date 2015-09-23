temp_kinomeData_dataProcessing <- function(df, iTRAQ_to_cellLine,runCondition ){
  
  #get the colnames which are other than the first 5 cols and not containing the variability and count
  ratio_cols_match_pattern = paste0("Variability|Count|", paste(colnames(df)[1:5], collapse="|"))
  ratio_cols <- colnames(df)[!grepl(ratio_cols_match_pattern,colnames(df),perl=T)]
  variability_cols <- colnames(df)[grepl("Variability",colnames(df),perl=T)]
  count_cols <- colnames(df)[grepl("Count",colnames(df),perl=T)]
  
  ratios <- melt(df, id.vars=colnames(df)[1:5], measure.vars = ratio_cols, variable.name = 'sample' , value.name = 'ratio')
  ratio_variability <- melt(df,id.vars=colnames(df)[1:5], measure.vars = variability_cols, variable.name = 'sample' , value.name = 'variability')
  ratio_variability$sample <- gsub(' Var.*','', ratio_variability$sample,perl=T)
  psm_counts <-   melt(df,id.vars=colnames(df)[1:5], measure.vars = count_cols, variable.name = 'sample' , value.name = 'count')
  psm_counts$sample <- gsub(' Count.*','', psm_counts$sample,perl=T)
  
  #create the final data frame
  merge.all <- function(x,y){ merge(x,y) }
  kinomeData <-Reduce(merge.all, list(ratios, ratio_variability, psm_counts))
  condition <- lapply(strsplit(as.character(kinomeData$sample), split='/'),
                      function(x) paste0(iTRAQ_to_cellLine[[x[[1]]]], '/', iTRAQ_to_cellLine[[x[[2]]]]) )
  condition <- paste0(condition, ' ', runCondition)
  kinomeData['condition'] <- condition
  kinomeData
}


process_newKinomeData_processing <- function(df){
  reqd_cols <- c('Accession', 'Description', 'Gene', 'Uniprot', 'Family', 'Score', 'Coverage',
                 '# Proteins', '# Unique Peptides',  '# Peptides',  '# PSMs', 'Area',
                 '# AAs', 'MW [kDa]', 'calc. pI')
  measured_cols <- df[,!colnames(df) %in% reqd_cols]
  ratio_cols <- colnames(measured_cols)[!grepl("Variability|Count",colnames(measured_cols),perl=T)]
  variability_cols <- colnames(measured_cols)[grepl("Variability",colnames(measured_cols),perl=T)]
  count_cols <- colnames(measured_cols)[grepl("Count",colnames(measured_cols),perl=T)]
  
  ratios <- melt(df, id.vars=reqd_cols, measure.vars = ratio_cols, variable.name = 'sample' , value.name = 'ratio')
  ratio_variability <- melt(df,id.vars=reqd_cols, measure.vars = variability_cols, variable.name = 'sample' , value.name = 'variability')
  ratio_variability$sample <- gsub(' Var.*','', ratio_variability$sample,perl=T)
  psm_counts <-   melt(df,id.vars=reqd_cols, measure.vars = count_cols, variable.name = 'sample' , value.name = 'count')
  psm_counts$sample <- gsub(' Count.*','', psm_counts$sample,perl=T)
  merge.all <- function(x,y){ merge(x,y) }
  kinomeData <-Reduce(merge.all, list(ratios, ratio_variability, psm_counts))
  
  kinomeData <- kinomeData %>%
                  separate(sample, into=c('cellLine', 'replicate', 'drug', 'time'), sep="_")
  kinomeData
}



temp_get_newkinomeData <- function(synid){
  print(paste("processing", synid))
  f <- synGet(synid)
  df <- read.xls(f@filePath, sheet="kinases", check.names=F)
  process_newKinomeData_processing(df)
}


library(synapseClient)
library("gdata")
library("tidyr")
library("dplyr")
library("reshape2")
require("parallel")
library("plyr")
library("doMC")
registerDoMC(4)
synapseLogin()


############
#OLD format data
############
#process KINOME RUN 1 : FullSerum
kinomeRun1 <- synGet('syn2679231')
kinomeRun1 <- read.xls(kinomeRun1@filePath, sheet=2, header=T, check.names=F)
#keep only selected columns (75% co-isolation)
kinomeRun1 <- kinomeRun1[,c(1:14)]
kinomeRun1_iTRAQ_to_cellLine <- list('116'='Syn2','117'='Syn3','114'='Syn5_1','115'='Syn5_2')                                                    
kinomeRun1_condition <- 'FullSerum'
kinomeRun1 <- temp_kinomeData_dataProcessing(kinomeRun1, kinomeRun1_iTRAQ_to_cellLine, kinomeRun1_condition)

#process KINOME RUN 2 : Serum Free
kinomeRun2 <- synGet('syn2679230')
kinomeRun2 <- read.xls(kinomeRun2@filePath, sheet=2, header=T, check.names=F)
#keep only selected columns (75% co-isolation)
kinomeRun2 <- kinomeRun2[,c(1:14)]
kinomeRun2_iTRAQ_to_cellLine <- list('114'='Syn2','115'='Syn3','116'='Syn5_1','117'='Syn5_2')                                                    
kinomeRun2_condition <- 'SerumFree'
kinomeRun2 <- temp_kinomeData_dataProcessing(kinomeRun2, kinomeRun2_iTRAQ_to_cellLine, kinomeRun2_condition)



###############
# New format data (post May 15, 2015)
##############
#get all Single drug treatment for meningioma 
meningioma_ids = synapseQuery('select id,name from file where parentId == "syn4214447"')
meningioma_ids <- meningioma_ids[!grepl('peptides', meningioma_ids$file.name),]
meningioma_kinome_data <- lapply(meningioma_ids$file.id, temp_get_newkinomeData)
meningioma_kinome_data <- do.call(rbind, meningioma_kinome_data)

#get all Single drug treatment for schwanoma
schwannoma_ids = synapseQuery('select id,name from file where parentId == "syn4214452"')
schwannoma_ids <- schwannoma_ids[!grepl('peptides', schwannoma_ids$file.name),]
schwannoma_kinome_data <- lapply(schwannoma_ids$file.id, temp_get_newkinomeData)
#schwannoma_kinome_data <- ldply(schwannoma_ids$file.id, .fun=temp_get_newkinomeData, .parallel = T)
schwannoma_kinome_data <- do.call(rbind, schwannoma_kinome_data)


#########
# Final Kinome Data
#########
kinomeData <- rbind(meningioma_kinome_data, schwannoma_kinome_data)
kinomeData['log2ratio'] = log2(kinomeData$ratio)
kinomeData['log2ratio_max'] =  kinomeData$log2ratio    + kinomeData$log2ratio*(kinomeData$variability/100)
kinomeData['log2ratio_min'] =  kinomeData$log2ratio    - kinomeData$log2ratio*(kinomeData$variability/100)
kinomeData['uniq_peptides'] =  kinomeData['# Unique Peptides']
kinomeData['# Unique Peptides'] <- NULL

#creating a unique col to identify each condition of run
kinomeData['condition'] = apply(kinomeData[,c('cellLine', 'drug', 'replicate', 'time')], 1,function(x) paste(x, collapse='_') )

# Upload the data to synapse
parentId = 'syn4259360'
write.table(kinomeData, file="Synodos_compiled_kinomeData.csv", sep="\t")
used_ids <- c(meningioma_ids$file.id,schwannoma_ids$file.id)
tmp <- synStore(File("Synodos_compiled_kinomeData.csv", parentId = parentId),
                used=used_ids)
unlink("Synodos_compiled_kinomeData.csv")
