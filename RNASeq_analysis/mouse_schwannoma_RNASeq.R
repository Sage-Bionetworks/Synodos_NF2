#!/usr/bin/env Rscript



library("edgeR")
library("gap")
library("gdata")
library(data.table)
library("synapseClient")
library("plyr")
library("dplyr")
library("tidyr")

synapseLogin()

#GLOBAL VARS
readthreshold <- 30
filter_short_genes<-"y"

edgeRDEanalysis <- function(countdata, group, control, treatment, prefix, DEcontrast,
                            readthreshold, outputPath){
  
  keep <- ((rowSums(countdata[,control] > readthreshold) == length(control) ) |  (rowSums(countdata[,treatment] > readthreshold) == length(treatment)))
  cat(paste("Analyzing:",prefix,"\n",sep=""))
  print(colnames(countdata[,control]))
  print(colnames(countdata[,treatment]))
  print(group)
  cat(paste("readthreshold:",readthreshold,"\n",sep=""))
  print(summary(keep))
  counttableselected <- countdata[keep,]
  print(colnames(counttableselected))
  
  y <-DGEList(counttableselected,group=group)
  y <- calcNormFactors(y,method="TMM")
  
  design <- model.matrix(~ 0 + group,data=y$samples)
  
  print(design)
  
  pdf(file=paste(outputPath, '/', prefix,"_MDplot.logcpm.pdf",sep=""),width=15,height=8)
  plotMD(cpm(y, log=TRUE),column=1)
  abline(h=0,col="red",lty=2,lwd=2)
  dev.off()
  
  logcpm <-cpm(y,normalized.lib.sizes=TRUE, log=TRUE)
  pdf(file=paste(outputPath, '/', prefix,"_MDSplot.logcpm.pdf", sep=""),width=15,height=8)
  par(mfrow=c(1,2))
  
  plotMDS(logcpm,top=10000,labels=groups,dim.plot=c(1,2))
  plotMDS(logcpm,top=10000,labels=groups,dim.plot=c(2,3))
  dev.off()
  
  
  pdf( paste(prefix,"_samplesimilarity_edgeR.pdf",sep="") , width = 7 , height = 7 )
  plotMDS.DGEList(y,labels=colnames(y$counts))
  dev.off()
  
  y <- estimateDisp(y,design,robust=TRUE)
  print(y$common.dispersion)
  
  pdf(paste(outputPath, '/',prefix,"_BCVplot.pdf",sep=""), width = 7 , height = 7 )
  plotBCV(y)
  dev.off()
  
  pdf(paste(outputPath, '/',prefix,"_MeanVariance.pdf",sep=""), width = 7 , height = 7 )
  meanvar <- plotMeanVar(y, show.raw.vars=TRUE,show.tagwise.vars=TRUE,show.ave.raw.vars=FALSE,NBline=TRUE,main = "Mean-Variance Plot")
  plotMeanVar(y, meanvar=meanvar, show.tagwise.vars=TRUE, NBline=TRUE)
  dev.off()
  
  fit <- glmQLFit(y,design)
  head(fit$coefficients)
  
  pdf(paste(outputPath, '/',prefix,"_QLDisplot.pdf",sep=""), width = 7 , height = 7 )
  plotQLDisp(fit)
  dev.off()
  
  qlf <- glmQLFTest(fit,contrast=DEcontrast)
  qlf$table$BH = p.adjust(qlf$table$PValue,"BH")
  qlf$table$bonferroni = p.adjust(qlf$table$PValue,"bonferroni")
  
  is.de <- decideTestsDGE(qlf, p.value=0.05)
  summary(is.de)
  png(filename=paste(outputPath, '/',prefix,"_edgeR_qlf_qqplot.png",sep=""),width=800,height=600)
  qqunif(qlf$table$PValue)
  dev.off()
  
  png(filename=paste(outputPath, '/',prefix,"_edgeR_qlf_plotsmear.png",sep=""),width=800,height=600)
  plotSmear(qlf, de.tags=rownames(qlf)[is.de!=0])
  dev.off()
  
  #write.table(qlf$table,file=paste(prefix,"_edgeR_quasilikelihoodFtest.txt",sep=""),sep="\t",quote=FALSE)
  qlf$table['gene'] = rownames(qlf$table)
  qlf$table
}


#batch 5
b5_metaData <- read.xls(synGet("syn7436870")@filePath)
b5_metaData <- b5_metaData[1:24,]
b5_metaData <- b5_metaData[, c(1:4)]
b5_metaData <- b5_metaData %>% mutate(Sample.ID = as.character(b5_metaData$Sample.ID)) %>%
  mutate(treatment = gsub('CUDC-907' ,'CUDC907', treatment)) %>%
  tidyr::separate(treatment, into=c('treatment','replicate'), sep='-') %>%
  mutate(tmp = Sample.ID) %>%
  tidyr::separate(tmp, into=c('sampleName','tmp'), sep='-')


#batch5 counts
b5_counts <- fread(synGet("syn7467412")@filePath, data.table=F)
b5_counts <- b5_counts[,-c(1:3)]
rownames(b5_counts) <- b5_counts$gene
b5_counts$gene <- NULL
colnames(b5_counts) <- gsub('_','-',colnames(b5_counts))
#b5_counts <- b5_counts[2:nrow(b5_counts),]

head(b5_counts)

#group variable
b5_counts <- b5_counts[,b5_metaData$Sample.ID]

groups <- paste0(b5_metaData$sampleName, '-',b5_metaData$treatment)
#available contrasts
colnames(model.matrix( ~ 0 + groups))


############################
### MS03 DMSO vs MS12 DMSO
############################
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'DMSO') %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'DMSO') %>% .$Sample.ID
prefix <- "MS03-DMSO_vs_MS12-DMSO"
contrast = c(0,1,0,0,0,-1,0,0)
MS03DMSO_vs_MS12DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                          outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03DMSO_vs_MS12DMSO <- MS03DMSO_vs_MS12DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS12',
         treatment1 = 'DMSO',
         treatment2 = 'DMSO')

############################
### MS03  vs MS12  for CUDC
############################
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'CUDC907') %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'CUDC907') %>% .$Sample.ID
prefix <- "MS03-CUDC_vs_MS12-CUDC"
colnames(model.matrix( ~ 0 + groups))

contrast = c(1,0,0,0,-1,0,0,0)
MS03CUDC_vs_MS12CUDC <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                        outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03CUDC_vs_MS12CUDC <- MS03CUDC_vs_MS12CUDC %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>% 
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS12',
         treatment1 = 'CUDC',
         treatment2 = 'CUDC')


############################
### MS03  vs MS12  for GSK2126458
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'GSK2126458') %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'GSK2126458') %>% .$Sample.ID
prefix <- "MS03-GSK458_vs_MS12-GSK458"

contrast = c(0,0,1,0,0,0,-1,0)
MS03GSK458_vs_MS12GSK458 <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                            outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03GSK458_vs_MS12GSK458 <- MS03GSK458_vs_MS12GSK458 %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS12',
         treatment1 = 'GSK458',
         treatment2 = 'GSK458')


############################
### MS03  vs MS12  for Pano
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'Panobin') %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'Panobin') %>% .$Sample.ID
prefix <- "MS03Pano_vs_MS12Pano"

contrast = c(0,0,0,1,0,0,0,-1)
MS03Pano_vs_MS12Pano <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                            outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03Pano_vs_MS12Pano <- MS03Pano_vs_MS12Pano %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS12',
         treatment1 = 'Pano',
         treatment2 = 'Pano')


############################
### MS03 for DMSO vs CUDC
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'CUDC907') %>% .$Sample.ID
prefix <- "MS03CUDC_vs_MS03DMSO"

contrast = c(1,-1,0,0,0,0,0,0)
MS03CUDC_vs_MS03DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                        outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03CUDC_vs_MS03DMSO <- MS03CUDC_vs_MS03DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS03',
         treatment1 = 'CUDC',
         treatment2 = 'DMSO')


############################
### MS03 for DMSO vs GSK
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'GSK2126458') %>% .$Sample.ID
prefix <- "MS03GSK458_vs_MS03DMSO"

contrast = c(0,-1,1,0,0,0,0,0)
MS03GSK458_vs_MS03DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                        outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03GSK458_vs_MS03DMSO <- MS03GSK458_vs_MS03DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS03',
         treatment1 = 'GSK458',
         treatment2 = 'DMSO')

############################
### MS03 for DMSO vs Panobinostat
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS03' & treatment == 'Panobin') %>% .$Sample.ID
prefix <- "MS03Panobin_vs_MS03DMSO"

contrast = c(0,-1,0,1,0,0,0,0)
MS03Panobin_vs_MS03DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                          outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS03Panobin_vs_MS03DMSO <- MS03Panobin_vs_MS03DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS03',
         cellLine2 = 'MS03',
         treatment1 = 'Pano',
         treatment2 = 'DMSO')

############################
### MS12 for DMSO vs CUDC
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'CUDC907') %>% .$Sample.ID
prefix <- "MS12CUDC_vs_MS12DMSO"

contrast = c(0,0,0,0,1,-1,0,0)
MS12CUDC_vs_MS12DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                        outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS12CUDC_vs_MS12DMSO <- MS12CUDC_vs_MS12DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS12',
         cellLine2 = 'MS12',
         treatment1 = 'CUDC',
         treatment2 = 'DMSO')


############################
### MS12 for DMSO vs GSK
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'GSK2126458') %>% .$Sample.ID
prefix <- "MS12GSK458_vs_MS12DMSO"

contrast = c(0,0,0,0,0,-1,1,0)
MS12GSK458_vs_MS12DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                          outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS12GSK458_vs_MS12DMSO <- MS12GSK458_vs_MS12DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS12',
         cellLine2 = 'MS12',
         treatment1 = 'GSK458',
         treatment2 = 'DMSO')

############################
### MS12 for DMSO vs Panobinostat
############################
colnames(model.matrix( ~ 0 + groups))
controlSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == "DMSO") %>% .$Sample.ID
caseSamples <- filter(b5_metaData, sampleName == 'MS12' & treatment == 'Panobin') %>% .$Sample.ID
prefix <- "MS12Panobin_vs_MS12DMSO"

contrast = c(0,0,0,0,0,-1,0,1)
MS12Panobin_vs_MS12DMSO <- edgeRDEanalysis(b5_counts,groups,controlSamples,caseSamples,prefix,contrast,readthreshold,
                                           outputPath="~/Documents/Github/Synodos_NF2/RNASeq_analysis/plots")
MS12Panobin_vs_MS12DMSO <- MS12Panobin_vs_MS12DMSO %>%
  tidyr::separate(gene, into=c('ensemblId', 'geneName'), sep="\\|") %>%
  mutate(cellLine1 = 'MS12',
         cellLine2 = 'MS12',
         treatment1 = 'Pano',
         treatment2 = 'DMSO')


mouse_schwannoma_diffExpgenes <- rbind(MS03Pano_vs_MS12Pano, MS03GSK458_vs_MS12GSK458,
                                       MS03CUDC_vs_MS12CUDC, MS03DMSO_vs_MS12DMSO, 
                                       MS03CUDC_vs_MS03DMSO, MS03GSK458_vs_MS03DMSO,
                                       MS03Panobin_vs_MS03DMSO, MS12CUDC_vs_MS12DMSO,
                                       MS12GSK458_vs_MS12DMSO, MS12Panobin_vs_MS12DMSO)



########
# upload the data
########
outFile = "Synodos_mouseSchwannoma_RNASeq_diffExp_genes.tsv"
write.table(mouse_schwannoma_diffExpgenes, file=outFile, sep="\t", col.names=T, row.names = F)
synStore(File(outFile, parentId = 'syn6004175'),
         used = c('syn7436870','syn7467412'),
         executed = "https://Github.com/Sage-Bionetworks/Synodos_NF2/blob/master/RNASeq_analysis/mouse_schwannoma_RNASeq.R")
unlink(outFile)

