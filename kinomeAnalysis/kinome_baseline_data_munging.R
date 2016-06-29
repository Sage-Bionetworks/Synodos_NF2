library(synapseClient)
library("gdata")
library("tidyr")
library("dplyr")
library("reshape2")
require("parallel")
library("plyr")
library("doMC")
library("gdata")
registerDoMC(4)
synapseLogin()


#schwannoma baseline data
baseline_schwannoma_human_synid <- "syn4214458"
baseline_schwannoma_human <- synGet(baseline_schwannoma_human_synid)
baseline_schwannoma_human <- read.xls(xls=baseline_schwannoma_human@filePath,sheet="trimmed")
rownames(baseline_schwannoma_human) <- baseline_schwannoma_human$Gene
baseline_schwannoma_human <- baseline_schwannoma_human[,grepl('LFQ',colnames(baseline_schwannoma_human))]
baseline_schwannoma_human$Gene <- toupper(rownames(baseline_schwannoma_human))

baseline_schwannoma_mouse_synid <- "syn4214457"
baseline_schwannoma_mouse <- synGet(baseline_schwannoma_mouse_synid)
baseline_schwannoma_mouse <- read.xls(xls=baseline_schwannoma_mouse@filePath,sheet="trimmed")
rownames(baseline_schwannoma_mouse) <- baseline_schwannoma_mouse$Gene
baseline_schwannoma_mouse <- baseline_schwannoma_mouse[,grepl('LFQ',colnames(baseline_schwannoma_mouse))]
baseline_schwannoma_mouse$Gene <- toupper(rownames(baseline_schwannoma_mouse))


#vestibular schwannoma data
baseline_vest_schwannoma_synid <- "syn4942532"
baseline_vest_schwannoma <- synGet(baseline_vest_schwannoma_synid)
baseline_vest_schwannoma <- read.xls(xls=baseline_vest_schwannoma@filePath,sheet="Kinase LFQ")
rownames(baseline_vest_schwannoma) <- baseline_vest_schwannoma$Gene.names
baseline_vest_schwannoma <- baseline_vest_schwannoma[,grepl('LFQ',colnames(baseline_vest_schwannoma))]
baseline_vest_schwannoma$Gene <- toupper(rownames(baseline_vest_schwannoma))
View(baseline_vest_schwannoma)


#meningioma baseline data
baseline_meningioma_synid <- "syn3104723"
baseline_meningioma <- synGet(baseline_meningioma_synid)
baseline_meningioma <- read.xls(xls=baseline_meningioma@filePath,sheet="Kinases with LFQ in all")
rownames(baseline_meningioma) <- baseline_meningioma$Gene.names
baseline_meningioma <- baseline_meningioma[,grepl('LFQ',colnames(baseline_meningioma))]
baseline_meningioma$Gene <- toupper(rownames(baseline_meningioma))

kinome_baseline_LFQ_data  <- Reduce(function(x,y) merge(x,y,by="Gene", all=T), list(baseline_meningioma, baseline_schwannoma_human,baseline_schwannoma_mouse, baseline_vest_schwannoma))
colnames(kinome_baseline_LFQ_data) <- gsub('LFQ.intensity.', '',colnames(kinome_baseline_LFQ_data))

#upload to synapse
outfile <- "Synodos_kinome_baseline_LFQ_data.tsv"
write.table(kinome_baseline_LFQ_data, file=outfile, sep="\t", row.names=F, col.names=T, quote=F)
synStore(File(outfile, parentId = "syn4259360"),
         used = c(baseline_meningioma_synid, baseline_schwannoma_human_synid, baseline_schwannoma_mouse_synid,
                  baseline_vest_schwannoma_synid),
         executed = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/kinomeAnalysis/kinome_baseline_data_munging.R")
unlink(outfile)


###
#process new baseline data - iTRAQ based
###
source("kinomeAnalysis/kinome_data_functions.R")
#new baseline triplicate data - protein level
meningioma_baseline_triplicates_proteins_synid <- 'syn5575201'
schwannoma_baseline_triplicates_proteins_synid <- 'syn5575204'
meningioma_baseline_triplicates_proteins <- get_new_kinome_proteinData(meningioma_baseline_triplicates_proteins_synid)
schwannoma_baseline_triplicates_proteins <- get_new_kinome_proteinData(schwannoma_baseline_triplicates_proteins_synid)
Kinome_baseline_iTRAQ_data <- rbind(meningioma_baseline_triplicates_proteins,schwannoma_baseline_triplicates_proteins)
Kinome_baseline_iTRAQ_data <- Kinome_baseline_iTRAQ_data %>% mutate(condition = sample) %>%
  separate(sample, into=c('cellLine', 'referenceSample'), sep="/") %>%
  separate(cellLine, into=c('cellLine', 'replicate'), sep="_") %>%
  separate(referenceSample, into=c('referenceSample', 'refSampleReplicate'), sep="_") %>%
  mutate(log2ratio = log2(ratio))

View(Kinome_baseline_iTRAQ_data)

#upload to synapse
outfile <- "Synodos_kinome_baseline_iTRAQ_data.tsv"
write.table(Kinome_baseline_iTRAQ_data, file=outfile, sep="\t", row.names=F, col.names=T, quote=F)
synStore(File(outfile, parentId = "syn4259360"),
         used = c(meningioma_baseline_triplicates_proteins_synid,
                  schwannoma_baseline_triplicates_proteins_synid),
         executed = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/kinomeAnalysis/kinome_baseline_data_munging.R")
unlink(outfile)

