library("synapseClient")
library("GSVA")

#####
#custom function: to do pathway enrichment analysis 
pathway_enrichment <- function(expMat, pathways, outfile_prefix, annotation=NULL, method="gsva"){
  res <- gsva(as.matrix(expMat),pathways , method=method)
  if(method == "gsva"){
    res <- res$es.obs
  }
  outfile <- paste0(base_dir, geo_id, '/', geo_id, outfile_prefix, '_enrichment.tsv')
  write.table(res, outfile, col.names=T, row.names=T, quote=FALSE, sep="\t")
  pathway_analysis <- synStore(Folder(name='pathway_analysis', parentId=study_synId))
  
  #store the pathway enrichment matrices in synapse
  outfile <- synStore(File(outfile, parentId =pathway_analysis$properties$id ),
                      used = geneExp_synId,
                      executed = thisScript )
  
  #plot the heatmap
  source("~/dev/apRs/expression_heatmap.R")
  res_scaled <- t(scale(t(res)))
  res_plot <- paste0(base_dir, geo_id, '/', geo_id, outfile_prefix, '_enrichment_heatmap.png')
  png(res_plot,width=10,height=8,units="in",res=300)
  memoised_pheatmap(res_scaled,
                    scale="none",
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    clustering_method = "complete",
                    border_color = NA,
                    fontsize_row=0,
                    fontsize_col=0,
                    annotation=annotation,
                    corr_method="spearman"
  )
  dev.off()  
  #store the heatmap in synapse with provenance
  res_plot <- synStore(File(res_plot, parentId =pathway_analysis$properties$id ),
                       used = outfile$properties$id,
                       executed = thisScript )
}

synapseLogin()

#get the MsigDB object
cat('Reading the MSIGDB object from synapse...')
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
cat('..Done\n\n')

pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)


thisScript <- "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/R_Scripts/2_pathway_enrichment.R"
base_dir <- "/external-data/DAT_123__Synodos_NF2/publicData//expressionData/"

#replace the following with the relevant study file id's from synapse
geo_id <- 'GSE39645'
study_synId <- 'syn2550097'
geneExp_synId <- 'syn2554578'
phenotype_synId <- 'syn2550100'


#download the normalized expression study
geneExp <- synGet(geneExp_synId)
geneExp <- read.table(geneExp@filePath, sep="\t", header=T)

phenotype <- synGet(phenotype_synId)
phenotype <- read.table(phenotype@filePath, sep="\t",header=T)


apply(phenotype,2,table)
#create an annotation dataframe for coloring
annotation <- phenotype[,c('batch'),drop=F]
rownames(annotation) <- phenotype$sampleId

#modify the column name to match the rowname of the annotation dataframe
colnames(geneExp)
colnames(geneExp) <- gsub('_.*.CEL.gz','',colnames(geneExp))


# pathway enrichment using GSVA
pathway_enrichment(geneExp, pathways_list, '_pathways_GSVA', annotation, method='gsva')

# Biocarta pathway enrichment using ssGSEA
pathway_enrichment(geneExp, pathways_list, '_pathway_ssGSEA', annotation, method='ssgsea')


