options(stringsAsFactors = F)

library(nplr)
library(dplyr)
library(rCharts)
library(ggvis)
library(ggplot2)
library(grid)
library(reshape2)
library(gridExtra)
library("synapseClient")
library('sjPlot')
synapseLogin()

source("~/dev/apRs/plotting_helpers.R")

UCF_DS_cleanData <- 'syn2773796'
UCF_DS <- synGet(UCF_DS_cleanData)
UCF_DS <- read.delim(UCF_DS@filePath, sep="\t", check.names=F)
#filter out any rows where DRUG = ''
UCF_DS <- UCF_DS[!(is.na(UCF_DS$drug) | UCF_DS$drug == ""),]
#modify colnames
UCF_DS <- melt(UCF_DS, id.vars=c('name', '0.1% DMSO', 'Ctrl', 'Untreated', '0.1% H2O', 'row',
                                 'drug', 'no cells', 'synid', 'experiment', 'cellLine'))
UCF_DS$experiment <- gsub(' ', '', UCF_DS$experiment)

#create a numeric conc column
UCF_DS['conc'] <- as.numeric(gsub(' uM','',UCF_DS$variable)) * 1e-6
UCF_DS$variable <- NULL
#create a numeric column for T0
UCF_DS['T0'] = as.numeric(UCF_DS[,'0.1% DMSO'])
#where .1% DMSO is NA which mainly for the Drug Perifosine which is soluble in water
UCF_DS[is.na(UCF_DS$T0),'T0'] <- UCF_DS[is.na(UCF_DS$T0),'0.1% H2O']

#convert to numeric
UCF_DS$Untreated <- as.numeric(UCF_DS$Untreated)
UCF_DS['viability'] <- as.numeric(UCF_DS$value)
UCF_DS$value <- NULL
#remove any rows where cell viability is not defined
UCF_DS <- UCF_DS[!is.na(UCF_DS$viability),]

#create a column for plate
plates <- factor(UCF_DS$name)
levels(plates)  <- paste0('plate-', 1:length(unique(UCF_DS$name)))
UCF_DS['plate'] <- plates

#fix drug names
UCF_DS$drug[UCF_DS$drug %in% c('AR-42')] = 'AR42'
UCF_DS$drug[UCF_DS$drug %in% c('GDC 0941', 'GCD0941')] = 'GDC0941'
UCF_DS$drug[UCF_DS$drug %in% c('OSU-03012')] = 'OSU03012'
UCF_DS$drug[UCF_DS$drug %in% c('Ganetspib')] = 'Ganetespib'


#######
#Normalization
#######

#convert vialbility to Zscore by using negative control like Untreated or DMSO/H20/T0
norm_by_untreated <- function(df){
  m <- mean(df$Untreated)
  stdev <- sd(df$Untreated)  
  df$viability <- (df$viability - m )/stdev
  df$Untreated <- (df$Untreated - m )/stdev
  df$T0 <- (df$T0 - m )/stdev
  df
}

#taking untreated as 100% viability
norm_by_untreated2 <- function(df){
  m <- mean(df$Untreated)  
  df$normViability <- df$viability/m
  df$norMUntreated <- df$Untreated/m
  df$T0 <- df$T0/m
  df$Ctrl <- df$Ctrl/m
  df
}

#taking DMSO as 100% viability
norm_by_meanDMSO <- function(df){
  m <- mean(df$T0)
  df['normViability'] <- df$viability/m
  df['meanDMSO'] <- m
  df$Untreated <- df$Untreated/m
  df$T0 <- df$T0/m
  df$Ctrl <- df$Ctrl/m
  df
}

UCF_DS_norm_by_meanDMSO <- ddply(.data = UCF_DS, .variables = c('plate', 'drug'),
                                 .fun = norm_by_meanDMSO)
colnames(UCF_DS_norm_by_meanDMSO)
#reshaping the data frame to match a common structure between MGH and UCF data
#dropping the unwanted cols
drop_cols <- c('name', '0.1% DMSO', '0.1% H2O', 'Ctrl', 'Untreated', 'T0', 'no cells', 'synid')
UCF_DS_norm_by_meanDMSO['replicate'] <- UCF_DS_norm_by_meanDMSO$row
UCF_DS_norm_by_meanDMSO$row <- NULL
new_col_order <- c('drug', 'conc', 'replicate', 'viability', 'cellLine', 'experiment', 'meanDMSO', 'normViability', 'plate')
UCF_DS_norm_by_meanDMSO <- UCF_DS_norm_by_meanDMSO[ , !colnames(UCF_DS_norm_by_meanDMSO) %in% drop_cols]
UCF_DS_norm_by_meanDMSO <- UCF_DS_norm_by_meanDMSO[, new_col_order]
write.table(ds_norm_by_meanDMSO, file="UCF_DrugScreen_DMSONorm_data.tsv", col.names=T, 
            sep="\t", quote=F, row.names=F)
synStore(File("UCF_DrugScreen_DMSONorm_data.tsv", parentId="syn2773788"), 
         used = UCF_DS_cleanData,
         executed = )











##################
# QC
##################
ds_run1 <- filter(ds, run=='Run1')
ds_run2 <- filter(ds, run=='Run2')
run1_drugDist <- table(ds_run1$drug, ds_run1$cellLine) 
run2_drugDist <- table(ds_run2$drug, ds_run2$cellLine) 
#convert to #replicates / drug (each drug has 9 concentration time points)
run1_drugDist <- round(run1_drugDist/9,0)
run2_drugDist <- round(run2_drugDist/9,0)
tableToPlot(run1_drugDist)
ggsave("plots/UCF_DrugScreen_run1_drugs.png", width=5, height=8, units="in")
tableToPlot(run2_drugDist)
ggsave("plots/UCF_DrugScreen_run2_drugs.png", width=5, height=8, units="in")

ds_run1 <- filter(ds, run=='Run1')
ds_run2 <- filter(ds, run=='Run2')
run1_drugDist <- table(ds_run1$drug, ds_run1$cellLine) 
run2_drugDist <- table(ds_run2$drug, ds_run2$cellLine) 
#convert to #replicates / drug (each drug has 9 concentration time points)
run1_drugDist <- round(run1_drugDist/9,0)
run2_drugDist <- round(run2_drugDist/9,0)
tableToPlot(run1_drugDist)
ggsave("plots/UCF_DrugScreen_run1_corrected_drugs.png", width=5, height=8, units="in")
tableToPlot(run2_drugDist)
ggsave("plots/UCF_DrugScreen_run2_corrected_drugs.png", width=5, height=8, units="in")

#get the order of plates arranged by increasing median DMSO viability
med_T0 <- tapply(ds$T0, as.character(ds$plate), median)
new_plate_levels <- names(sort(med_T0))

temp_transform_func <- function(x,logged) { if(logged == 'TRUE'){log10(x)} else {x} }

drugScreen_QC_plot <- function(df,plate_levels,ylab, logy=T){
  df$plate <- factor(df$plate,plate_levels)
  
  df$T0 <- temp_transform_func(df$T0,logy)
  df$Untreated <- temp_transform_func(df$Untreated,logy)
  df$Ctrl <- temp_transform_func(df$Ctrl,logy)
  
  p1 <- ggplot(df, aes(x=plate, y=T0, fill=cellLine)) + geom_boxplot()
  p1 <- p1 + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) 
  p1 <- p1 + ylab(ylab) + xlab('') + ggtitle('DMSO/H20')
  
  #boxplot of Positive Control variability across the plates using the same levels as of T0 for fair comparison
  p2 <- ggplot(df, aes(x=plate, y=Ctrl, fill=cellLine)) + geom_boxplot()
  p2 <- p2 + theme(axis.text.x=element_text(angle=90, hjust=1))
  p2 <- p2 + ylab(ylab) + xlab('plate number') + ggtitle('Positive Control')
  
  #boxplot of Untreated(negative control)  variability across the plates using the same levels as of T0 for fair comparison
  p3 <- ggplot(df, aes(x=plate, y=Untreated, fill=cellLine)) + geom_boxplot()
  p3 <- p3 + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
  p3 <- p3 +  ylab(ylab) + ggtitle('Negative Control(untreated)') + xlab('')
  multiplot(p3,p1,p2)
}



#drug screen QC post Normalization
ds_norm_by_meanUntreated <- ddply(.data = ds, .variables = c('plate', 'drug'),
                                  .fun = norm_by_untreated2)
drugScreen_QC_plot(ds_norm_by_meanUntreated, new_plate_levels,'viability(%)', logy=F)


drugScreen_QC_plot(ds, new_plate_levels, 'viability(log10)', logy=TRUE)

#cell viability across drug concentrations
ds$conc <- round(ds$conc,3)
new_levels <- as.character(sort(unique(ds$conc)))
p <- ggplot(ds, aes(x=factor(ds$conc, new_levels), y=log10(viability) , fill=cellLine))
p <- p + geom_boxplot() + xlab('conc(log10)') + ggtitle('cell viability v/s dosage')
p + theme_bw()


######
#between replicates variaiton across each plates
######

temp_CalcVar <- function(df){
  viab_var <- tapply(log10(df$viability), df$conc,sd)
  cols <- rownames(viab_var)
  viab_var <- data.frame(t(viab_var))
  colnames(viab_var) <- cols
  viab_var['untreated'] <- sd(log(df$Untreated))
  viab_var['T0'] <- sd(log10(df$T0))
  viab_var  
}

calc_meanViab <- function(df){
  mean_viab <- tapply(log10(df$viability), df$conc,mean)
  cols <- rownames(mean_viab)
  mean_viab <- data.frame(t(mean_viab))
  colnames(mean_viab) <- cols
  mean_viab['untreated'] <- mean(log(df$Untreated))
  mean_viab['T0'] <- mean(log10(df$T0))
  mean_viab  
}


var_between_replicates <- ddply(ds, .variables =c('drug', 'run', 'cellLine', 'run'), .fun = temp_CalcVar )
scaled_var_btwn_replicates <- apply(var_between_replicates[,-c(1:3)],2,scale)
scaled_var_btwn_replicates <- cbind(var_between_replicates[,c(1:3)],scaled_var_btwn_replicates)


meanViab <- ddply(ds, .variables =c('drug', 'run', 'cellLine', 'run'), .fun = calc_meanViab )
meanViab$untreated <- NULL
meanViab$T0 <- NULL
#meanViab <- cbind(meanViab[,c(1:4)] , scale(meanViab[,-c(1:4)]))
meanViab <- melt(meanViab, id.vars = c('drug', 'run', 'cellLine', 'run'))
ggplot(meanViab, aes(x=variable, y=drug,)) + geom_tile(aes(fill=value)) + facet_grid(~ cellLine)

#Edge variation
plate_untreated <- ds[,c('plate', 'drug', 'Untreated','T0', 'row', 'run', 'cellLine')]
plate_untreated <- plate_untreated[!duplicated(plate_untreated),]
temp_func <- function(x){
    unique(x)[1]
}
plate_untreated_cols <- dcast(plate_untreated, row ~ plate + drug , 
                              value.var="Untreated",
                              fun.aggregate=temp_func)
#remove cols which have all NA as some plates don't have untreated cells viability
plate_untreated_cols <- plate_untreated_cols[,!apply(plate_untreated_cols,2,function(x) all(is.na(x)))]
rownames(plate_untreated_cols) <- plate_untreated_cols$row
plate_untreated_cols$row <- NULL
source("~/dev/apRs/expression_heatmap.R")
m <- log10(plate_untreated_cols)
m.scaled <- t(scale(t(m)))
annotation <- plate_untreated
annotation <- annotation[!duplicated(annotation[,c('plate','drug')]),]
rownames(annotation) <- paste0(annotation$plate,'_', annotation$drug)
annotation <- annotation[,c('run','cellLine'),drop=F]
memoised_pheatmap(m.scaled,annotation=annotation,
                  scale='none',
                  cluster_rows = FALSE,
                  clustering_distance_cols = "euclidean",
                  border_color = NA)



#Drug viability across dosages
plot_drugViab_vs_dosage <- function(df){
  drugViab_dosages <- dcast(df, drug+plate+row+cellLine+Untreated+T0 ~ conc, value.var="viability",
                            fun.aggregate = function(x) mean(x))
  m <- drugViab_dosages[,-c(1:4)]
  m.scaled <- t(scale(t(m)))
  memoised_pheatmap(m.scaled,
                    scale='none',
                    cluster_cols = FALSE,
                    clustering_distance_cols = "euclidean",
                    border_color = NA)
}

png("plots/UCF_DrugScreens_drugViab_vs_dosage_heatmap.png",
    width=10, height=8, units="in", res=300)
plot_drugViab_vs_dosage(ds_norm_by_meanDMSO)
dev.off()


plot_drugViab_vs_dosage(ds)
plot_drugViab_vs_dosage(ds_norm_by_meanUntreated)
