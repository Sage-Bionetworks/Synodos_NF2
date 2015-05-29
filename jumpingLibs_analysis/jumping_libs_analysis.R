# Synodos : Jumping library analysis 
#library("devtools")
#devtools::install_github("andrewhzau/intansv")
library("intansv")
library("synapseClient")
#read in the mod 
source("intansv_mod_readDelly_code.R")

data_dir <- "~/testing/synodos_jumping_lib_delly_SV/"
#dir.create(data_dir)


#download the SV data
synapseLogin()
dels <- synGet('syn3917437', downloadLocation=data_dir)
dups <- synGet('syn3917464', downloadLocation=data_dir)
invs <- synGet('syn3917446', downloadLocation=data_dir)

delly <- readdelly_local(dataDir=data_dir, regSizeLowerCutoff=100, 
          regSizeUpperCutoff=10000000, readsSupport=5,
          method="DELLY", pass=TRUE, minMappingQuality=30)

#combine all the delly SVs
delly_sv <- do.call(rbind,delly)
delly_sv$type <- gsub('\\..*','',rownames(delly_sv))
rownames(delly_sv) <- 1:nrow(delly_sv)


library(GenomicRanges)
library(Homo.sapiens)
geneRanges <- function(db, column="ENTREZID"){
  g <- genes(db, columns=column)
  col <- mcols(g)[[column]]
  genes <- granges(g)[rep(seq_along(g), elementLengths(col))]
  mcols(genes)[[column]] <- as.character(unlist(col))
  genes
}

splitColumnByOverlap <-function(query, subject, column="ENTREZID", ...){
  olaps <- findOverlaps(query, subject, ...)
  f1 <- factor(subjectHits(olaps),
               levels=seq_len(subjectLength(olaps)))
  splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
}


delly_gr <- makeGRangesFromDataFrame(delly_sv, keep.extra.columns = T,
                                     ignore.strand = FALSE,
                                     seqnames.field = 'chromosome',
                                     start.field = 'pos1',
                                     end.field = 'pos2')


hg19_genes <- geneRanges(Homo.sapiens, column="SYMBOL")
hg19_genes <- keepStandardChromosomes(hg19_genes)
seqlevelsStyle(hg19_genes) <- "NCBI"

genes_in_sv = splitColumnByOverlap(hg19_genes, delly_gr, "SYMBOL")
genes_in_sv <- as.vector(genes_in_sv)
delly_sv['genes'] = as.character(lapply(genes_in_sv,function(x) paste(x,collapse=", ")))

paste(unique(delly_sv$genes), collapse=',')

#generate HTMLreport
library(ReportingTools)
htmlRep <- HTMLReport(shortName = "Synodos Jumping Library preliminary analysis",
                      reportDirectory = "./reports")
publish(delly_sv,htmlRep)


#chr plots
temp_plot <- function(i){
  chr <- seqlevels(hg19_genes)[i]
  size <- seqlengths(hg19_genes)[i]
  gr  <- GRanges(seqnames = Rle(chr),ranges   = IRanges(1,size))
  plotname <- paste0('SV_chr',chr,'.png')
  plotChromosome(gr, delly, windowSize=1000000)
  ggsave(filename=paste0("reports/plots/",plotname), width=5, height=4,dpi=200,
         units="in")
  return(plotname)
}

plots <- lapply(1:length(seqlengths(hg19_genes)),temp_plot)
plots <- unlist(plots)
plots <- paste0("plots/",plots)

library(hwriter)
html_plots <- hwriteImage(matrix(plots,ncol=5),link=matrix(plots,ncol=5),
                          width=200, height=200)
publish("<br><br>",htmlRep)
publish("<b>Structural Variation plots per chromosome</b>", htmlRep)
publish("<i>Blue: deletions, Red: duplications, Green:inversions </i>", htmlRep)
publish(html_plots,htmlRep)
finish(htmlRep)



source("http://bioconductor.org/biocLite.R")
biocLite("ReportingTools")



#get list of all the genes 
# MSIGDB_syn<-synGet("syn2227979")
# load(MSIGDB_syn@filePath) #available as MSigDB R object
# pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
# pathways_list

###annotation
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# hg19_genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
# #hg19_genes <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
# seqlevelsStyle(hg19_genes) <- "NCBI"
# hg19_genes <- keepStandardChromosomes(hg19_genes)

# library(org.Hs.eg.db)
# entrezid_key<- mcols(hg19_genes)$gene_id
# hg19_gene_names <- select(org.Hs.eg.db, keys=entrezid_key, columns=c("SYMBOL","GENENAME"), 
#                           keytype="ENTREZID")
# values(hg19_genes) <- cbind(values(hg19_genes), extra_annot_cols)


#intersect with genes
# library(biomaRt)
# mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# regions <- paste(delly_sv$chromosome,delly_sv$pos1,delly_sv$pos2,sep=":")
# delly_sv_annot <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
#                                        "start_position", "end_position"),
#                         filters = c("chromosomal_region"),
#                         values = regions,
#                         mart=mart)
# 
# 
# delly.anno <- llply(delly, svAnnotation, genomeAnnotation=hg19_genes)
# svAnnotation
# 



  

chr <- hg19_genome_size[23]
chr
chr_name <- as.character(seqnames(chr)@values)[1]
plotname <- paste0('SV_chr',chr_name,'.png')
png(filename=plotname)
plotChromosome(chr, delly, windowSize=1000000)
chr
dev.off()
levels(seqnames(chr)) <- chr_name



getwd()



x <- seqnames(gr[1])
lapply(gr,function(x) colnames(x))



# library(biomaRt)
# mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
#                                 "start_position", "end_position"),
#                  filters = c("chromosome_name", "start", "end"), 
#                  values=list(1, 94312388, 96000000),
#                  mart=mart)