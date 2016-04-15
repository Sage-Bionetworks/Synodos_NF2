#combined pathway heatmap

library("synapseClient")
library("GSVA")
library("reshape")

#source the heatmap code
source("~/dev//apRs/expression_heatmap.R")  


thisScriptURL = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/2_combined_studies_heatmap.R"
  

#subset to select only the GSVA enrichment files
filter_df <- function(df,pattern){
  temp_df <- df[with(df,grepl(pattern,entity.name)),]  
  temp_df['geo_study_id'] = gsub(paste0('_',pattern,'.*'),"",temp_df$entity.name)
  temp_df
}

process_enrichment_file <- function(x){
  geo_study_id <- as.character(x[3])
  enrichment_file <- synGet(x[2])
  enrichment_file <- read.table(enrichment_file@filePath,sep="\t",header=T)
  colnames(enrichment_file) <- paste0(geo_study_id,'_',colnames(enrichment_file))
  enrichment_file['pathway'] <- rownames(enrichment_file)
  enrichment_file
}

process_exp_file <- function(x){
  geo_study_id <- as.character(x[3])
  expFile <- synGet(x[2])
  expFile <- read.table(expFile@filePath,sep="\t",header=T)
  colnames(expFile) <- paste0(geo_study_id,'_',colnames(expFile))
  expFile['genes'] <- rownames(expFile)
  expFile
}


rank_enrichment_file <- function(x){
  geo_study_id <- as.character(x[3])
  study_file <- synGet(as.character(x[2]))
  study_file <- read.table(study_file@filePath,sep="\t",header=T)
  geometric_mean <- function(x) {exp(mean(log(x+1)))-1}
  study_rank <- as.data.frame(rank(apply(study_file,1,geometric_mean)))
  colnames(study_rank)=c(geo_study_id)
  study_rank['pathway'] = rownames(study_rank)
  study_rank
}

merge_studies <- function(studies_list,by){
  merge.all <- function(x,y){ merge(x,y,by=by,all=T) }
  studies_merged <-Reduce(merge.all, studies_list)
  #setup the row names
  rownames(studies_merged) <- studies_merged[[by]]
  studies_merged[by] <- NULL
  studies_merged
}


custom_heatmap <- function(mat,annotation){
  memoised_pheatmap(mat,
                    scale="none",
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    clustering_method = "average",
                    border_color = NA,
                    fontsize_row=3,
                    fontsize_col=2,
                    cor_method = "spearman",
                    annotation = annotation
  )
}



##########
#Main
##########

#login
synapseLogin()

#get all the entities under a project
synQuery_results <- synQuery("select id,name from entity where benefactorId == 'syn2347420'")

combined_studies_bucket <- "syn2536538"

pathways_enrichment_ssGSEA_files <- filter_df(synQuery_results,'pathways_GSVA_enrichment.tsv')
pathways_enrichment_GSVA_files <- filter_df(synQuery_results,'pathways_GSVA_enrichment.tsv')

setwd("/external-data//DAT_123__Synodos_NF2/publicData//expressionData/")

#1.
pathway_ssGSEA_enrichment_scores <- apply(pathways_enrichment_ssGSEA_files,1,process_enrichment_file)
pathway_ssGSEA_enrichment_scores <- merge_studies(pathway_ssGSEA_enrichment_scores,'pathway')
pathway_ssGSEA_enrichment_scores_scaled <- t(scale(t(pathway_ssGSEA_enrichment_scores)))


#prepare annotation
study_ids <- sapply(strsplit(colnames(pathway_ssGSEA_enrichment_scores),'_',fixed=T),"[[",1)
annotation <- data.frame('GEO_ids'= study_ids)
annotation['type'] = 'meningioma'
schwannoma_studies = c('GSE30563','GSE39645')
annotation$type[annotation$GEO_ids %in% schwannoma_studies] = 'schwannoma'
rownames(annotation) <- colnames(pathway_ssGSEA_enrichment_scores)

png("pathway_ssGSEA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(pathway_ssGSEA_enrichment_scores_scaled, annotation)
dev.off()

#2.
pathway_GSVA_enrichment_scores <- apply(pathways_enrichment_GSVA_files,1,process_enrichment_file)
pathway_GSVA_enrichment_scores <- merge_studies(pathway_GSVA_enrichment_scores, 'pathway')
write.table(pathway_GSVA_enrichment_scores, "pathways_GSVA_enrichment_scores.tsv", sep="\t",
            col.names=T, row.names=T )
temp <- synStore( File("pathways_GSVA_enrichment_scores.tsv", parentId=combined_studies_bucket),
                   used = pathways_enrichment_GSVA_files$entity.id,
                   executed = thisScriptURL
                   )

pathway_GSVA_enrichment_scores_scaled <- t(scale(t(pathway_GSVA_enrichment_scores)))
png("pathway_GSVA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(pathway_GSVA_enrichment_scores_scaled, annotation)
dev.off()


###################
# combined expression matrix
###################
expVals_files <- filter_df(synQuery_results,'expvals')
expVals <- apply(expVals_files, 1, process_exp_file)
expVals <- merge_studies(expVals, 'genes')

write.table(expVals, "combined_studies_expressionVals.tsv", sep="\t", col.names=T, row.names=T )
temp <- synStore( File("combined_studies_expressionVals.tsv", parentId=combined_studies_bucket),
                  used = expVals_files$entity.id,
                  executed = thisScriptURL
)



################
## based on pathway ranking
################

#1
BioCarta_ssGSEA_enrichment_ranks <- apply(BioCarta_ssGSEA_enrichment_files,1,rank_enrichment_file)
BioCarta_ssGSEA_enrichment_ranks <- merge_studies(BioCarta_ssGSEA_enrichment_ranks)
BioCarta_ssGSEA_enrichment_ranks_scaled <- t(scale(t(BioCarta_ssGSEA_enrichment_ranks)))

#prepare annotation
annotation_ranks <- data.frame('GEO_ids'= colnames(BioCarta_ssGSEA_enrichment_ranks))
annotation_ranks['type'] = 'meningioma'
schwannoma_studies = c('GSE30563','GSE39645')
annotation_ranks$type[annotation_ranks$GEO_ids %in% schwannoma_studies] = 'schwannoma'
rownames(annotation_ranks) <- colnames(BioCarta_ssGSEA_enrichment_ranks)
png("BioCarta_ssGSEA_pathway_ranks_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(BioCarta_ssGSEA_enrichment_ranks_scaled, annotation_ranks)
dev.off()

#2
BioCarta_GSVA_enrichment_ranks <- apply(BioCarta_GSVA_enrichment_files,1,rank_enrichment_file)
BioCarta_GSVA_enrichment_ranks <- merge_studies(BioCarta_GSVA_enrichment_ranks)
BioCarta_GSVA_enrichment_ranks_scaled <- t(scale(t(BioCarta_GSVA_enrichment_ranks)))
png("BioCarta_GSVA_pathway_ranks_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(BioCarta_GSVA_enrichment_ranks_scaled, annotation_ranks)
dev.off()







dim(KEGG_GSVA_enrichment_scores)
dim(BioCarta_GSVA_enrichment_scores)
x <- rbind(KEGG_GSVA_enrichment_scores_scaled,BioCarta_GSVA_enrichment_scores_scaled)

cols = as.factor(annotation$type)
levels(cols)[levels(cols)=="meningioma"] <- "blue"
levels(cols)[levels(cols)=="schwannoma"] <- "red"
plotMDS(x,col=as.character(cols))



y <- rbind(KEGG_GSVA_enrichment_ranks_scaled, BioCarta_GSVA_enrichment_ranks_scaled)
cols <- as.factor(annotation_ranks$type)
levels(cols)[levels(cols)=="meningioma"] <- "blue"
levels(cols)[levels(cols)=="schwannoma"] <- "red"
plotMDS(y,col=as.character(cols))
