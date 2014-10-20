options(stringsAsFactors = F)

library(plyr)
library(dplyr)
library(ggvis)
library(ggplot2)
library(grid)
library(reshape2)
library(gridExtra)
library("synapseClient")

synapseLogin()
setwd("~/dev/Synodos_NF2/R_Scripts/DrugScreening")
source("common_functions.R")
source("~/dev/apRs/plotting_helpers.R")

#get the cleaned MGH Drug Screen Data
MGH_cleaned_data <- 'syn2773611'
MGH_DS <- synGet(MGH_cleaned_data)
MGH_DS <- read.delim(MGH_DS@filePath,sep="\t", check.names=F)

#convert conc to molar from microMolar
MGH_DS$conc = as.numeric(MGH_DS$conc) * 1e-6
MGH_DS$cellLine = gsub('syn', 'Syn', MGH_DS$cellLine)

#normalize the drug names to be in syn with UCF AND MGH data
MGH_DS[grepl('Everol', MGH_DS$drug, ignore.case=T), 'drug'] = 'Everolimus'
MGH_DS[grepl('Selumetinib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Selumetinib'
MGH_DS[grepl('Bortezomib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Bortezomib'
MGH_DS[grepl('Lapatinib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Lapatinib'
MGH_DS[grepl('Panobinostat', MGH_DS$drug, ignore.case=T), 'drug'] = 'Panobinostat'
MGH_DS[grepl('Vorinostat', MGH_DS$drug, ignore.case=T), 'drug'] = 'Vorinostat'
MGH_DS[grepl('GDC-0941', MGH_DS$drug, ignore.case=T), 'drug'] = 'GDC0941'
MGH_DS[grepl('Vismodegib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Vismodegib'
MGH_DS[grepl('OSU-03012', MGH_DS$drug, ignore.case=T), 'drug'] = 'OSU03012'
MGH_DS[grepl('Ganetespib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Ganetespib'
MGH_DS[grepl('AR-42', MGH_DS$drug, ignore.case=T), 'drug'] = 'AR42'
MGH_DS[grepl('Trametinib', MGH_DS$drug, ignore.case=T), 'drug'] = 'Trametinib'
MGH_DS[grepl('GDC-0980', MGH_DS$drug, ignore.case=T), 'drug'] = 'GDC0980'
MGH_DS[grepl('CUDC-907', MGH_DS$drug, ignore.case=T), 'drug'] = 'CUDC907'
MGH_DS[grepl('Perifosine', MGH_DS$drug, ignore.case=T), 'drug'] = 'Perifosine'
MGH_DS[grepl('GSK2126458', MGH_DS$drug, ignore.case=T), 'drug'] = 'GSK2126458'


#nromalize by DMSO 
#IMP here to understand the plate structure, refer: syn2765805
#each cellLine across experiment was done on a single plate
MGH_meanDMSO <- MGH_DS %>%
                  group_by(cellLine, experiment ) %>%
                  summarize(meanDMSO = mean(viability[drug == 'DMSO']))

MGH_DS <- merge(MGH_DS, MGH_meanDMSO)
MGH_DS['normViability'] = MGH_DS$viability / MGH_DS$meanDMSO

#filter out DMSO and media values
MGH_DS <- filter(MGH_DS,  ! drug %in% c('DMSO', 'media'))


tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$conc, df$normViability, useLog=T)  
  },error=function(e){
      print(dim(df))
      print(df$conc)
      print(df$normViability)
      print(unique(df$cellLine))
      print(unique(df$drug))
      print(unique(df$experiment))
      stop('stopped')
  })
}

MGH_DS_dr_fit_by_meanDMSO <- ddply(.data=MGH_DS, .variables = c('experiment', 'cellLine', 'drug'), 
                                   .fun=tmp_iterator)
MGH_DS_dr_fit_by_meanDMSO_flt <- filter(MGH_DS_dr_fit_by_meanDMSO, goodNess_of_fit > .8)




#store in synapse
#DMSO norm data
write.table(MGH_DS, file="MGH_DrugScreen_DMSONorm_data.tsv", col.names=T, sep="\t", quote=F, row.names=F)
synStore(File("MGH_DrugScreen_DMSONorm_data.tsv", parentId="syn2740273"), 
         used = MGH_cleaned_data,
         executed = )


#ICvals
write.table(MGH_DS_dr_fit_by_meanDMSO_flt, file="MGH_DrugScreen_ICVals.tsv", col.names=T, sep="\t", quote=F, row.names=F)
synStore(File("MGH_DrugScreen_ICVals.csv", parentId="syn2740273"), used = )



















drug_levels <- ddply(.data=MGH_DS_dr_fit_by_meanDMSO_flt, .variables=c('drug'), .fun=function(x) mean(x$IC50,na.rm=T))
drug_levels <- arrange(drug_levels, desc(V1))
drug_levels <- drug_levels$drug

View(MGH_DS_dr_fit_by_meanDMSO)

#IC50 across three cell lines across two runs
p <- ggplot(data=MGH_DS_dr_fit_by_meanDMSO_flt, aes(x=factor(drug,levels=drug_levels), y=log10(IC50), group=experiment)) + geom_line(aes(color=experiment)) 
p <- p + facet_grid(cellLine ~ .) + theme_bw() + geom_point(color='grey50') 
p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug')



MGH_DS_dr_fit_by_meanDMSO_flt[which(MGH_DS_dr_fit_by_meanDMSO_flt$IC50 < 1e-25),]




doseResp <- filter(MGH_DS_dr_fit_by_meanDMSO_flt, drug %in% c('Lapatinib', 'Panobinostat'))

p <- ggplot(data=doseResp, aes(x=fittedX, y=fittedY, group=cellLine))
p <- p + geom_line(aes(color=cellLine)) + facet_grid(experiment ~ drug) + theme_bw()

p
p <- p + geom_hline(aes(yintercept=0.5), color='red3', linetype='dashed')
p <- p + xlab('molar conc (log10)') + ylab('cell viability %')
p




###########
#QC code
###########
exp1 <- filter(MGH_DS, experiment == 'exp_1')
exp2 <- filter(MGH_DS, experiment == 'exp_2')









