#install.packages(c("WGCNA", "cluster"))
library(WGCNA) 
#Allowing 6 threads for WGCNA
allowWGCNAThreads(6)
library(cluster) 
library(synapseClient)
options(stringsAsFactors = FALSE) 


setwd("/external-data/DAT_123__Synodos_NF2//publicData//expressionData/networkAnalysis/")


synapseLogin()

#get all the entities under a project
synQuery_results <- synQuery("select id,name from entity where benefactorId == 'syn2347420'")

process_exp_file <- function(x){
  geo_study_id <- as.character(x[3])
  expFile <- synGet(x[2])
  expFile <- read.table(expFile@filePath,sep="\t",header=T)
  colnames(expFile) <- paste0(geo_study_id,'_',colnames(expFile))
  expFile
}


#subset to select only the GSVA enrichment files
filter_df <- function(df,pattern){
  temp_df <- df[with(df,grepl(pattern,entity.name)),]  
  temp_df['geo_study_id'] = gsub(paste0('_',pattern,'.*'),"",temp_df$entity.name)
  temp_df
}

#get all the expression files
expVals_files <- filter_df(synQuery_results,'expvals')

schwannoma_studies = c('GSE30563','GSE39645')

#get only the meningioma studies
meningioma_expVals_files <- expVals_files[!expVals_files$geo_study_id %in% schwannoma_studies,]

#keep only two studies for consensus analysis
selectedStudies <- meningioma_expVals_files[meningioma_expVals_files$geo_study_id %in% c('GSE43290', 'GSE16581'),]


#filter out any study < 18k genes
# studies_toKeep <- unlist(lapply(meningioma_exp_mats, function(x) nrow(x) > 18000))
# meningioma_exp_mats_flt <- meningioma_exp_mats[studies_toKeep]
# studyLabels <- studyLabels[studies_toKeep]


studyLabels <- selectedStudies$geo_study_id


#list of the exp mats across all the studies
expMat <- apply(selectedStudies, 1, process_exp_file)


#find the list of common genes across all the filtered studies
common_genes <- (Reduce(intersect, lapply(expMat, function(x) rownames(x))))


#get only the common genes in all the studies and create a list as required
prep_data_for_WCGNA <- function(x){
  rows_to_keep <- rownames(x) %in% common_genes
  x <- x[rows_to_keep,]
  list(data=t(x))
}

multiExpr <- lapply(expMat, prep_data_for_WCGNA)

# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)
exprSize


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK



#cluster samples
sampleTrees = list()
for (set in 1:length(multiExpr))
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "complete")
}

#plot dendograms
pdf(file = "Meningioma_SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(length(multiExpr),1))
for (set in 1:length(multiExpr))
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", studyLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();




consensus_network = blockwiseConsensusModules(multiExpr,
                                 corType="bicor", 
                                 maxBlockSize=13000,
                                 networkType="signed",
                                 power=14,
                                 minModuleSize=30,
                                 numericLabels=TRUE,
                                 saveTOMs=TRUE,
                                 mergeCutHeight=0.25, 
                                 verbose=5,
                                 deepSplit = 2,
                                 pamRespectsDendro=FALSE)



save(consensus_network,"consensus_networkAnalysis.RData")


load("consensus_networkAnalysis.rda")

consMEs = consensus_network$multiMEs
moduleLabels = consensus_network$colors
moduleColors = labels2colors(moduleLabels)
consTree = consensus_network$dendrograms[[1]]

plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")


load("")





# 
# moduleLabels = consensus_network$colors
# table(moduleLabels)
# moduleColors = labels2colors(moduleLabels)
# consTree = net$dendrograms[[1]]
# datME = net$MEs
# 
# 
# plotDendroAndColors(consTree, moduleColors, "Module Colors", dendroLabels = F, hang = 0.03, addGuide=TRUE,
#                     guideHang = 0.05)
# 
# 
# cor(datME, use="p")
# 
# 
# dissimME=(1-t(cor(datME, method="p")))/2
# hclustdatME=hclust(as.dist(dissimME), method="average" )
# # Plot the eigengene dendrogram
# par(mfrow=c(1,1))
# plot(hclustdatME, main="Clustering tree based of the module eigengenes")
# 
# 
# 
# 
# 
