library(NMF)
library('stringr')
library('ggplot2')
library(synapseClient)
library(plyr)
library(dplyr)
library(data.table)
library('rbokeh')
library(reshape2)
library("combinat")


#login to synapse
synapseLogin()


#get the kinome protein data
kinome_protein_data <- synGet('syn4259365')
kinome_protein_data <- read.table(kinome_protein_data@filePath, sep="\t",
                                  check.name=F, header=T)

#get the kinome peptide data
kinome_peptide_data <- synGet('syn4378913')
kinome_peptide_data <- read.table(kinome_peptide_data@filePath,
                                  sep="\t",  check.name=F, header=T)
kinome_protein_data$Gene <- toupper(kinome_protein_data$Gene)
unique(kinome_protein_data$Gene)

#summary of data
summary <- kinome_protein_data %>%
  group_by(cellLine, drug, time ) %>%
  summarise(replicates=length(unique(replicate)))
p <- ggplot(data=summary, aes(x=cellLine, y=replicates)) + geom_point()
p <- p + facet_grid(time ~ drug) + scale_fill_brewer(palette = "Set2") + theme_bw()
p <- p + ggtitle("Kinome Data Summary") 
png("plots/kinomeData_summary.png", width=12, height=8, units="in", res=200)
p + ylab('number of replciates') +  scale_y_discrete(breaks=c(1,2,3,4)) + expand_limits(y=c(0,4))
dev.off()


# Heatmap of log2ratio
df <- filter(kinome_protein_data, cellLine %in% c('MS02','MS03', 'MS12', 'HS11', 'HS01'))
df <- dcast(df, Gene  ~ cellLine + drug + time + replicate , value.var = "log2ratio" )
df <- df[apply(df,1, function(x) sum(is.na(x)) == 0 ),]
rownames(df) <- as.character(df$Gene)
df$Gene <- NULL
m <- t(scale(t(df)))
plot(hclust(dist(t(m))))
annotCol = data.frame(matrix(unlist(strsplit(colnames(m), split="_")), ncol=4, byrow=T))
colnames(annotCol) = c('cellLine', 'drug', 'time' , 'replicate')
png('plots/heatmap_schwannoma_cellLines_kinome.png', width=12, height=8, units="in", res=200)
aheatmap(m,scale="none",cexCol = 4,cexRow=0, hclustfun="complete",annCol = annotCol[,c('time')])
dev.off()


#correlation between replicates
get_correlations <- function(df,...){
  gene_lists <- lapply(unique(df$replicate), function(x) filter(df, replicate == x) %>% select(Gene) %>% transform(Gene = unique(as.character(Gene))) %>%.$Gene)
  names(gene_lists) <- unique(df$replicate)
  combinations <- t(combn(unique(df$replicate),2))
  combinations <- split(combinations, 1:nrow(combinations))
  rep_cor <- sapply(combinations, function(x) {
    selected_gene_lists <- gene_lists[names(gene_lists) %in% x]
    common_genes <- Reduce(intersect, selected_gene_lists)
    tmp <- df %>% filter(Gene %in% common_genes & replicate %in% x) %>%
      dcast(Gene ~ replicate, value.var = "log2ratio") 
    cor(tmp[,2],tmp[,3], method="spearman")
  })
  results <- data.frame(combinations = sapply(combinations,function(x) paste(x, collapse='+')),
                        rep_cor = rep_cor)  
  results  
}
df <- filter(kinome_protein_data, cellLine %in% c('MS02', 'MS03', 'MS12', 'HS01', 'HS11'))
replicates_corr <- ddply(.data=df, .variables = c('cellLine', 'drug', 'time'),
                             .fun = get_correlations)
p <- ggplot(data=replicates_corr, aes(x=cellLine, y=rep_cor, color=combinations)) + geom_point()
p <- p + facet_grid(time ~ drug) + scale_fill_brewer(palette = "Set2") + theme_bw()
p <- p + ggtitle("between replicates correlation") + scale_color_brewer(palette = "Set1")
png('plots/corr_schwannoma_cellLines_kinome.png',width=12, height=8, units="in", res=200)
p + ylab('spearman correlation') + geom_hline(yintercept=0)
dev.off()


#finding correlations using proteins which have > 10 peptides in all 3 replicates
df <- df[df['# Peptides'] > 10,]
replicates_corr_1 <- ddply(.data=df, .variables = c('cellLine', 'drug', 'time'),
                         .fun = get_correlations)
p <- ggplot(data=replicates_corr_1, aes(x=cellLine, y=rep_cor, color=combinations)) + geom_point()
p <- p + facet_grid(time ~ drug) + scale_fill_brewer(palette = "Set2") + theme_bw()
p <- p + ggtitle("between replicates correlation") + scale_color_brewer(palette = "Set1")
png('plots/corr_peptides_gt10_schwannoma_cellLines_kinome.png',width=12, height=8, units="in", res=200)
p + ylab('spearman correlation') + geom_hline(yintercept=0)
dev.off()


#histogram of log2ratio
df <- filter(kinome_protein_data, cellLine %in% c('MS03', 'MS02', 'MS12', 'HS11', 'HS01'))
p <- ggplot(data=df, aes(x=log2ratio,fill=replicate))
p <- p + facet_grid(drug + time  ~ cellLine) + theme_bw() + scale_fill_brewer(palette = "Set1")
png('plots/hist_log2ratios.png',width=12, height=8, units="in", res=200)
p + geom_density(binwidth=.05,alpha=.5)
dev.off()

# percent overlaps between proteins for a cellLine
#compare replicates for a cellLine
get_overlaps <- function(df){
  gene_lists <- lapply(unique(df$replicate), function(x) filter(df, replicate == x) %>% select(Gene) %>% transform(Gene = unique(as.character(Gene))) %>%.$Gene)
  names(gene_lists) <- unique(df$replicate)
  
  combinations <- t(combn(unique(df$replicate),2))
  combinations <- split(combinations, 1:nrow(combinations))
  combinations <- c(combinations,list(names(gene_lists)))
  
  percent_overlap <- sapply(combinations, function(x) {
    selected_gene_lists <- gene_lists[names(gene_lists) %in% x]
    length(Reduce(intersect, selected_gene_lists)) / length(Reduce(union, selected_gene_lists))  
  })
  percent_overlap <- as.numeric(percent_overlap)
  results <- data.frame(combinations = sapply(combinations,function(x) paste(x, collapse='+')),
                        overlap = percent_overlap * 100)  
  results  
}
df <- filter(kinome_protein_data, cellLine %in% c('MS02', 'MS03', 'MS12', 'HS01', 'HS11'))
replicates_overlaps <- ddply(.data=df, .variables = c('cellLine', 'drug', 'time'),
      .fun = get_overlaps)
replicates_overlaps$combinations <- factor(replicates_overlaps$combinations, 
                                          levels = c("run1+run2", "run1+run3", "run2+run3", "run1+run2+run3"))
p <- ggplot(data=replicates_overlaps, aes(x=cellLine, y=overlap, fill=combinations) ) + geom_bar(stat="identity", position="dodge")
p <- p + facet_grid(time ~ drug) + scale_fill_brewer(palette = "Set2") + theme_bw()
png('plots/proteins_overlap.png',width=12, height=8, units="in", res=200)
p + ggtitle("common proteins % across replicates") + ylab('overlap %')
dev.off()




calc_correlation <- function(df){
  df <- na.omit(df)
  data.frame('spearman_cor' = cor(df$run1, df$run2, method="spearman"))
}
cellLine_corr <- ddply(.data = df, 
      .variables = c('cellLine', 'time', 'drug'),
      .fun = calc_correlation)


View(cellLine_corr)

library(ggplot2)
p <- ggplot(data=cellLine_corr, aes(x=cellLine, y=spearman_cor, color=drug)) + geom_point() 
p  + facet_grid( . ~ time) + theme_bw() +  geom_hline(aes(yintercept=0),color="darkred")
geom_point() + facet_wrap(time ~.) +


## overall cellLines / drugs and / timepoints
x <- dcast(kinome_protein_data, Accession ~ cellLine + replicate + drug + time ,
           value.var = 'log2ratio')
png('1.png')
hist(apply(x,1, function(y) sum(is.na(y))/ncol(x)), breaks=10,
     xlab='percent samples', ylab = '#kinases',
     main='hist coverage of kinome across 48 datasets')
dev.off()


### 2. pathway present in differentially exp proteins
MSIGDB_syn<-synGet("syn2227979")
load(MSIGDB_syn@filePath) #available as MSigDB R object
#pathways_list <- c(MSigDB$C2.CP.BIOCARTA, MSigDB$C2.CP.KEGG, MSigDB$C2.CP.REACTOME)
pathways_list <- c(MSigDB$C2.CP.KEGG)
pathway_gene_map <- ldply(pathways_list, function(x) {
  data.frame(gene=as.character(x))})
colnames(pathway_gene_map) <- c('pathway', 'gene')




