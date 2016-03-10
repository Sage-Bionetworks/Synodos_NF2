#diff expression analysis
library(synapseClient)
library("DESeq2")
synapseLogin()
library(plyr)
library(dplyr)
library(limma)
library(edgeR)
library(gplots)
library(gdata)
library(tidyr)
library(reshape2)
library(WGCNA)

source("~/dev/apRs/plotting_helpers.R")
source("~/dev/apRs/expression_heatmap.R")
source("~/dev/apRs/dataExploration_plots.R")

prepareData <- function(selected_cellLines, selected_treatments){
  flt_metadata <- mdata %>% filter(cellLine %in% selected_cellLines & Treatment %in% selected_treatments)
  flt_data <- data[,colnames(data) %in% flt_metadata$SampleID]
  flt_metadata['replicates'] = paste0(flt_metadata$cellLine,'_',flt_metadata$Treatment)
  #filter low counts and zero counts genes
  flt_data <- flt_data[complete.cases(flt_data),]
  flt_data <- flt_data[rowSums(flt_data) / ncol(flt_data) > 1, ]
  flt_data <- flt_data[rowSums(flt_data) > 30,] #30 selected arbitrarily
  df1=data.frame(SampleID = colnames(flt_data))
  df2=merge(df1, flt_metadata)
  colnames(flt_data) <- df2$sampleAlias
  list(flt_data, flt_metadata)
}
runLimma <- function(flt_data,flt_metadata){
  dge <- DGEList(counts=flt_data)
  dge <- dge[rowSums(dge$counts) > 30, keep.lib.size = FALSE]
  dge <- calcNormFactors(dge)
  #model matrix 
  flt_metadata$cellLine <- as.factor(flt_metadata$cellLine)
  flt_metadata$Treatment <- as.factor(flt_metadata$Treatment)
  treatment <- factor(paste(flt_metadata$cellLine, flt_metadata$Treatment, sep="."))
  #design = model.matrix( ~ 0 + cellLine + Treatment  , data=flt_metadata)
  design = model.matrix( ~0+treatment)
  #data norm using voom
  dge_norm <- voom(dge, design=design,plot=F)
  #fit the model
  fit <- lmFit(dge_norm, design, weights = dge_norm$weights)
  fit
}

processLimmaResults <- function(limma.fit.results){
  diffExpResults <- lapply(colnames(limma.fit.results$contrasts), function(x){
    df <- topTable(fit.new, coef=x, number = Inf)
    df['diffExptest'] = x
    df$geneName <- rownames(df)
    rownames(df) <- NULL
    df
  })
  diffExpResults <- do.call(rbind, diffExpResults)
  diffExpResults$diffExptest <- gsub('treatment', '', diffExpResults$diffExptest)
  diffExpResults <- diffExpResults %>% mutate(temp = diffExptest) %>% separate(temp, into=c('cellLine1', 'cellLine2'), sep='-') %>% 
    separate(cellLine1, into=c('cellLine1', 'treatment1')) %>% separate(cellLine2, into=c('cellLine2', 'treatment2'))
}

getLimmaSummary <- function(df, pval){
  summary <- df %>% group_by(cellLine1, cellLine2, treatment1, treatment2) %>%
    summarise(diffExpGenes_pVals = sum(adj.P.Val <= pval),
              UP_pVals = sum(logFC > 0 & adj.P.Val <= pval), 
              DOWN_pVals = sum(logFC < 0 & adj.P.Val <= pval),
              diffExpGenes_foldChange = sum(abs(logFC) > 1 ),
              UP_foldChange = sum(logFC > 1), 
              DOWN_foldChange = sum(logFC < -1 ))
}

#download mdata
mdata = synTableQuery("select * from syn5569975", loadResult = T)
mdata = mdata@values
mdata$Treatment = gsub('-','',mdata$Treatment)

#downlaod data
data = read.table(synGet("syn5568205")@filePath, header=T, sep="\t", as.is = T, check.names = F )
rownames(data) <- data[,1]
data[,1] <- NULL  

#remove low quality samples
low_quality_samples <- c('AC027', 'AC028', 'AC029', 'AC034', 'MN529', 'MN533', 'MN567')
mdata <- mdata[!mdata$cellLine %in% low_quality_samples,]
data <- data[ , colnames(data) %in% mdata$SampleID]

##################
#diff expression
##################
#1. HS01 vs HS11
selected_cellLines <- c('HS11', 'HS01')
selected_treatments <- c('CUDC907', 'Panobinostat', 'GSK2126458', 'DMSO')
prepd_data <- prepareData(selected_cellLines, selected_treatments)
flt_data <- prepd_data[[1]]
flt_metadata <- prepd_data[[2]]  
fit <- runLimma(flt_data, flt_metadata)
contrasts <- makeContrasts('treatmentHS01.Panobinostat-treatmentHS01.DMSO',
                           'treatmentHS01.GSK2126458-treatmentHS01.DMSO',
                           'treatmentHS01.CUDC907-treatmentHS01.DMSO',
                           'treatmentHS11.Panobinostat-treatmentHS11.DMSO',
                           'treatmentHS11.GSK2126458-treatmentHS11.DMSO',
                           'treatmentHS11.CUDC907-treatmentHS11.DMSO',
                           'treatmentHS01.Panobinostat-treatmentHS11.Panobinostat',
                           'treatmentHS01.GSK2126458-treatmentHS11.GSK2126458',
                           'treatmentHS01.CUDC907-treatmentHS11.CUDC907',
                           'treatmentHS01.DMSO-treatmentHS11.DMSO',
                           levels=fit$design)
fit.new <- contrasts.fit(fit, contrasts)
fit.new <- eBayes(fit.new)
HS01_vs_HS11 <- processLimmaResults(fit.new)
summary_HS01_vs_HS11 <- getLimmaSummary(HS01_vs_HS11, pval=.10)
png("plots/table.png",width = 10, height = 3, units = "in", res=200)
grid.arrange(tableGrob(summary_HS01_vs_HS11[,-c(8,9,10)]))
dev.off()
p <- lapply(unique(HS01_vs_HS11$treatment1), 
            function(x){
              length(x)
              d <- HS01_vs_HS11 %>% filter(treatment1 == x )
              ggplot(data=d, aes(x=logFC, fill=diffExptest)) + geom_density(alpha=.5) + facet_grid(~treatment1) + theme_bw()
})
p[[1]]
ggsave(filename = "plots/temp_plot.png", width = 10, height = 8, units = "in", dpi=200)
png('plots/temp_plot.png', width=15, height=12, units="in", res=300)
multiplot(plotlist = p, cols=2)
dev.off()

# volacno plot
p <- ggplot(data=HS01_vs_HS11 %>% filter(cellLine1 == cellLine2), aes(x=logFC, y=-log10(adj.P.Val))) 
p <- p + geom_jitter(alpha=.3, color="grey50", size=.5, fill="grey") + facet_grid(cellLine1~treatment1)
p <- p + theme_bw()
ggsave(filename = "plots/temp_plot.png", plot = p, width = 10, height = 8, units = "in", dpi=200)

# PCA plot // by samples
m <- t(scale(t(flt_data)))
colorBy = flt_metadata$Treatment
names(colorBy) = flt_metadata$sampleAlias
doPCA(t(m), colorBy = colorBy)

# correlation heatmap
png('plots/corPlot_2.png', width=15, height=12, units="in", res=300)
expHeatMap(cor(flt_data), fontsize=12)
dev.off()

#clustering
png("plots/temp_plot.png",width = 12, height = 8, units = "in", res=200)
doClustering(m,cex=1)
dev.off()

#################################################

#2. Syn1 vs Syn5 
selected_cellLines <- c('Syn1', 'Syn5')
selected_treatments <- c('CUDC907', 'Panobinostat', 'GSK2126458', 'DMSO')
prepd_data <- prepareData(selected_cellLines, selected_treatments)
flt_data <- prepd_data[[1]]
flt_metadata <- prepd_data[[2]]  
fit <- runLimma(flt_data, flt_metadata)
contrasts <- makeContrasts(#Syn5 vs Syn1
                          'treatmentSyn1.Panobinostat-treatmentSyn1.DMSO',
                          'treatmentSyn1.GSK2126458-treatmentSyn1.DMSO',
                          'treatmentSyn1.CUDC907-treatmentSyn1.DMSO',
                          'treatmentSyn5.Panobinostat-treatmentSyn5.DMSO',
                          'treatmentSyn5.GSK2126458-treatmentSyn5.DMSO',
                          'treatmentSyn5.CUDC907-treatmentSyn5.DMSO',
                          'treatmentSyn5.Panobinostat-treatmentSyn1.Panobinostat',
                          'treatmentSyn5.GSK2126458-treatmentSyn1.GSK2126458',
                          'treatmentSyn5.CUDC907-treatmentSyn1.CUDC907',
                          'treatmentSyn5.DMSO-treatmentSyn1.DMSO',
                          levels=fit$design)
fit.new <- contrasts.fit(fit, contrasts)
fit.new <- eBayes(fit.new)
Syn5_vs_Syn1 <- processLimmaResults(fit.new)
summary_Syn5_vs_Syn1 <- getLimmaSummary(Syn5_vs_Syn1,pval=.10)
p <- lapply(unique(Syn5_vs_Syn1$treatment1), 
            function(x){
              length(x)
              d <- Syn_lines %>% filter(treatment1 == x )
              ggplot(data=d, aes(x=logFC, fill=diffExptest)) + geom_density(alpha=.5) + facet_grid(~treatment1) + theme_bw()
            })
p[[1]]
ggsave(filename = "plots/temp_plot.png", width = 10, height = 8, units = "in", dpi=200)
png('plots/temp_plot.png', width=10, height=8, units="in", res=300)
multiplot(plotlist = p, cols=2)
dev.off()

#2. Syn1 vs Syn6 
selected_cellLines <- c('Syn1', 'Syn6')
selected_treatments <- c('CUDC907', 'Panobinostat', 'GSK2126458', 'DMSO')
prepd_data <- prepareData(selected_cellLines, selected_treatments)
flt_data <- prepd_data[[1]]
flt_metadata <- prepd_data[[2]]  
fit <- runLimma(flt_data, flt_metadata)
contrasts <- makeContrasts(#Syn6 vs Syn1
                          'treatmentSyn6.Panobinostat-treatmentSyn6.DMSO',
                          'treatmentSyn6.GSK2126458-treatmentSyn6.DMSO',
                          'treatmentSyn6.CUDC907-treatmentSyn6.DMSO',
                          'treatmentSyn6.Panobinostat-treatmentSyn1.Panobinostat',
                          'treatmentSyn6.GSK2126458-treatmentSyn1.GSK2126458',
                          'treatmentSyn6.CUDC907-treatmentSyn1.CUDC907',
                          'treatmentSyn6.DMSO-treatmentSyn1.DMSO',
                          levels=fit$design)
fit.new <- contrasts.fit(fit, contrasts)
fit.new <- eBayes(fit.new)
Syn6_vs_Syn1 <- processLimmaResults(fit.new)
summary_Syn6_vs_Syn1 <- getLimmaSummary(Syn6_vs_Syn1,pval=.10)


#now doing the QC plots together for 3 cell lines
selected_cellLines <- c('Syn1', 'Syn5', 'Syn6')
selected_treatments <- c('CUDC907', 'Panobinostat', 'GSK2126458', 'DMSO')
prepd_data <- prepareData(selected_cellLines, selected_treatments)
flt_data <- prepd_data[[1]]
flt_metadata <- prepd_data[[2]]  

# PCA plot // by samples
m <- t(scale(t(flt_data)))
colorBy = flt_metadata$cellLine
names(colorBy) = flt_metadata$sampleAlias
doPCA(t(m), colorBy = colorBy)
ggsave(filename = "plots/temp_plot.png", width = 4, height = 4, units = "in", dpi=100)

png('plots/corPlot_2.png', width=15, height=12, units="in", res=300)
expHeatMap(cor(flt_data), fontsize=12)
dev.off()

#clustering
png("plots/temp_plot.png",width = 12, height = 8, units = "in", res=200)
doClustering(m,cex=1)
dev.off()


#################
#Upload Diff exp data to synapse
################
diff_exp_genes <- rbind(HS01_vs_HS11, Syn5_vs_Syn1, Syn6_vs_Syn1)

outfile <- "Synodos_RNASeq_diffExp_data.tsv"
write.table(diff_exp_genes, file=outfile, sep="\t", col.names=T, row.names = F)
synStore(File(outfile, parentId = 'syn4259360'),
         used = c('syn5568205','syn5569975'),
         executed = "")
unlink(outfile)

#create a synapse table
tab <-as.tableColumns(diff_exp_genes)
cols<-tab$tableColumns
fileHandleId <- tab$fileHandleId
projectId <-"syn2347420"
schema<-TableSchema(name=outfile, parent=projectId, columns=cols)
table<-Table(schema, fileHandleId)
table<-synStore(table, retrieveData=TRUE)


