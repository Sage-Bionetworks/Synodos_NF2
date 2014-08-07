#combined pathway heatmap

library("synapseClient")
library("GSVA")
library("reshape")

#source the heatmap code
source("~/dev//apRs/expression_heatmap.R")  

#login
synapseLogin()

combined_studies_bucket <- "syn2536538"

thisScriptURL = ""


#get all the entities under a project
synQuery_results <- synQuery("select id,name from entity where benefactorId == 'syn2347420'")

#subset to select only the GSVA enrichment files
filter_df <- function(df,pattern){
  temp_df <- df[with(df,grepl(pattern,entity.name)),]  
  temp_df['geo_study_id'] = gsub(paste0('_',pattern),"",temp_df$entity.name)
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

merge_studies <- function(studies_list){
  merge.all <- function(x,y){ merge(x,y,by="pathway") }
  studies_merged <-Reduce(merge.all, studies_list)
  #setup the row names
  rownames(studies_merged) <- studies_merged$pathway
  studies_merged$pathway <- NULL
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


#KEGG pathway files
KEGG_ssGSEA_enrichment_files <- filter_df(synQuery_results,'KEGG_ssGSEA_enrichment.tsv')
KEGG_GSVA_enrichment_files <- filter_df(synQuery_results,'KEGG_GSVA_enrichment.tsv')

#BioCarta enrichment files
BioCarta_ssGSEA_enrichment_files <- filter_df(synQuery_results,'BioCarta_ssGSEA_enrichment.tsv')
BioCarta_GSVA_enrichment_files <- filter_df(synQuery_results,'BioCarta_GSVA_enrichment.tsv')


setwd("/external-data//DAT_123__Synodos_NF2/publicData//expressionData/")

#1.
KEGG_ssGSEA_enrichment_scores <- apply(KEGG_ssGSEA_enrichment_files,1,process_enrichment_file)
KEGG_ssGSEA_enrichment_scores <- merge_studies(KEGG_ssGSEA_enrichment_scores)
KEGG_ssGSEA_enrichment_scores_scaled <- t(scale(t(KEGG_ssGSEA_enrichment_scores)))


#prepare annotation
study_ids <- sapply(strsplit(colnames(KEGG_ssGSEA_enrichment_scores),'_',fixed=T),"[[",1)
annotation <- data.frame('GEO_ids'= study_ids)
annotation['type'] = 'meningioma'
schwannoma_studies = c('GSE30563','GSE39645')
annotation$type[annotation$GEO_ids %in% schwannoma_studies] = 'schwannoma'
rownames(annotation) <- colnames(KEGG_ssGSEA_enrichment_scores)

png("KEGG_ssGSEA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(KEGG_ssGSEA_enrichment_scores_scaled, annotation)
dev.off()

#2.
KEGG_GSVA_enrichment_scores <- apply(KEGG_GSVA_enrichment_files,1,process_enrichment_file)
KEGG_GSVA_enrichment_scores <- merge_studies(KEGG_GSVA_enrichment_scores)
write.table(KEGG_GSVA_enrichment_scores, "KEGG_GSVA_enrichment_scores.tsv", sep="\t",
            col.names=T, row.names=T )
_temp <- synStore( File("KEGG_GSVA_enrichment_scores.tsv", parentId=combined_studies_bucket),
                   used = KEGG_GSVA_enrichment_files$geo_study_id
                   executed = thisScriptURL
                   )

KEGG_GSVA_enrichment_scores_scaled <- t(scale(t(KEGG_GSVA_enrichment_scores)))
png("KEGG_GSVA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(KEGG_GSVA_enrichment_scores_scaled, annotation)
dev.off()


KEGG_GSVA_enrichment_files

#3. 
BioCarta_GSVA_enrichment_scores <- apply(BioCarta_GSVA_enrichment_files,1,process_enrichment_file)
BioCarta_GSVA_enrichment_scores <- merge_studies(BioCarta_GSVA_enrichment_scores)
BioCarta_GSVA_enrichment_scores_scaled <- t(scale(t(BioCarta_GSVA_enrichment_scores)))
png("BioCarta_GSVA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(BioCarta_GSVA_enrichment_scores_scaled, annotation)
dev.off()


#4. 
BioCarta_ssGSEA_enrichment_scores <- apply(BioCarta_ssGSEA_enrichment_files,1,process_enrichment_file)
BioCarta_ssGSEA_enrichment_scores <- merge_studies(BioCarta_ssGSEA_enrichment_scores)
BioCarta_ssGSEA_enrichment_scores_scaled <- t(scale(t(BioCarta_ssGSEA_enrichment_scores)))
png("BioCarta_ssGSEA_enrichment_scores_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(BioCarta_ssGSEA_enrichment_scores_scaled,annotation)
dev.off()



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

#3
KEGG_GSVA_enrichment_ranks <- apply(KEGG_GSVA_enrichment_files,1,rank_enrichment_file)
KEGG_GSVA_enrichment_ranks <- merge_studies(KEGG_GSVA_enrichment_ranks)
KEGG_GSVA_enrichment_ranks_scaled <- t(scale(t(KEGG_GSVA_enrichment_ranks)))
png("KEGG_GSVA_pathway_ranks_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(KEGG_GSVA_enrichment_ranks_scaled, annotation_ranks)
dev.off()

#4
KEGG_ssGSEA_enrichment_ranks <- apply(KEGG_ssGSEA_enrichment_files,1,rank_enrichment_file)
KEGG_ssGSEA_enrichment_ranks <- merge_studies(KEGG_ssGSEA_enrichment_ranks)
#fitler out low variance
KEGG_ssGSEA_enrichment_ranks <- KEGG_ssGSEA_enrichment_ranks[!apply(KEGG_ssGSEA_enrichment_ranks,1,sd) == 0,]
KEGG_ssGSEA_enrichment_ranks_scaled <- t(scale(t(KEGG_ssGSEA_enrichment_ranks)))
png("KEGG_ssGSEA_pathway_ranks_heatmap.png",width=10,height=8,units="in",res=300)
custom_heatmap(KEGG_ssGSEA_enrichment_ranks_scaled, annotation_ranks)
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
