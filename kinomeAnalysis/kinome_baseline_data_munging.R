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
baseline_schwannoma_human$Gene <- rownames(baseline_schwannoma_human)


baseline_schwannoma_mouse_synid <- "syn4214457"
baseline_schwannoma_mouse <- synGet(baseline_schwannoma_mouse_synid)
baseline_schwannoma_mouse <- read.xls(xls=baseline_schwannoma_mouse@filePath,sheet="trimmed")
rownames(baseline_schwannoma_mouse) <- baseline_schwannoma_mouse$Gene
baseline_schwannoma_mouse <- baseline_schwannoma_mouse[,grepl('LFQ',colnames(baseline_schwannoma_mouse))]
baseline_schwannoma_mouse$Gene <- rownames(baseline_schwannoma_mouse)


#vestibular schwannoma data
baseline_vest_schwannoma_synid <- "syn4942532"
baseline_vest_schwannoma <- synGet(baseline_vest_schwannoma_synid)
baseline_vest_schwannoma <- read.xls(xls=baseline_vest_schwannoma@filePath,sheet="Kinase LFQ")
rownames(baseline_vest_schwannoma) <- baseline_vest_schwannoma$Gene.names
baseline_vest_schwannoma <- baseline_vest_schwannoma[,grepl('LFQ',colnames(baseline_vest_schwannoma))]
baseline_vest_schwannoma$Gene <- rownames(baseline_vest_schwannoma)

#meningioma baseline data
baseline_meningioma_synid <- "syn3104723"
baseline_meningioma <- synGet(baseline_meningioma_synid)
baseline_meningioma <- read.xls(xls=baseline_meningioma@filePath,sheet="Kinases with LFQ in all")
rownames(baseline_meningioma) <- baseline_meningioma$Gene.names
baseline_meningioma <- baseline_meningioma[,grepl('LFQ',colnames(baseline_meningioma))]
baseline_meningioma$Gene <- rownames(baseline_meningioma)

kinome_baseline_data  <- Reduce(function(x,y) merge(x,y,by="Gene", all=T), list(baseline_meningioma, baseline_schwannoma_human,baseline_schwannoma_mouse, baseline_vest_schwannoma))
colnames(kinome_baseline_data) <- gsub('LFQ.intensity.', '',colnames(kinome_baseline_data))

#upload to synapse
outfile <- "Synodos_kinome_baseline_data.tsv"
write.table(kinome_baseline_data, file=outfile, sep="\t", row.names=F, col.names=T, quote=F)
synStore(File(outfile, parentId = "syn4259360"),
         used = c(baseline_meningioma_synid, baseline_schwannoma_human_synid, baseline_schwannoma_mouse_synid,
                  baseline_vest_schwannoma_synid),
         executed = )
