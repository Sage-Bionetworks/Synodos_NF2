library(synapseClient)
library(tidyr)
library(dplyr)
library(tibble)
library(Biobase)

##### HUMAN BASELINE
kin.base <- read.table(synGet("syn6182638")@filePath, sep = "\t", header = TRUE, comment.char = "") %>% 
  filter(cellLine %in% c("HS01","Syn5"), referenceSample %in% c("HS11", "Syn1")) %>% 
  select(Gene, cellLine, replicate, referenceSample, refSampleReplicate, log2ratio) %>% 
  unite(comp, cellLine, referenceSample) %>% 
  group_by(Gene, comp) %>% 
  dplyr::summarise(avgLog2Norm = mean(log2ratio)) %>% 
  spread(comp, avgLog2Norm) %>% 
  remove_rownames() %>% 
  as.data.frame() %>% 
  column_to_rownames("Gene")

saveRDS(kin.base, "baseline_kinome.rds")

this.file = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/kinomeAnalysis/kinomeDataforSyDE.R"
synStore(File("baseline_kinome.rds", parentId = "syn10845573"), used = c("syn6181167"), executed = this.file)


##### HUMAN TREATED 
kinome <- read.table(synGet("syn6181167")@filePath, sep = "\t", header = TRUE) %>% 
  filter(time == "24h", cellLine %in% c("HS01","HS11","Syn1","Syn5","Syn6")) %>%
  select(cellLine, drug, replicate, Gene, log2NormRatio) %>% 
  group_by(cellLine, drug, Gene) %>% 
  dplyr::summarize(avgLog2NormRatio = mean(log2NormRatio)) %>% 
  unite(sample, cellLine, drug, sep = "_") %>% 
  spread(sample, avgLog2NormRatio) %>% 
  remove_rownames() %>% 
  as.data.frame() %>% 
  column_to_rownames("Gene")

kinome <- as.matrix(kinome)

meta <- as.data.frame(colnames(kinome))
colnames(meta) <- "sampleNames"
meta$Genotype <- NA
meta$Treatment <- NA
meta$`Cell_Line` <- NA
meta$`Cell_Type` <- NA

meta$Genotype[grep("HS01", meta$sampleNames)] <- "-/-"
meta$Genotype[grep("Syn5", meta$sampleNames)] <- "-/-"
meta$Genotype[grep("Syn6", meta$sampleNames)] <- "-/-"
meta$Genotype[grep("HS11", meta$sampleNames)] <- "+/+"
meta$Genotype[grep("Syn1", meta$sampleNames)] <- "+/+"

meta$`Cell_Line`[grep("HS01", meta$sampleNames)] <- "HS01"
meta$`Cell_Line`[grep("Syn5", meta$sampleNames)] <- "Syn5"
meta$`Cell_Line`[grep("Syn6", meta$sampleNames)] <- "Syn6"
meta$`Cell_Line`[grep("HS11", meta$sampleNames)] <- "HS11"
meta$`Cell_Line`[grep("Syn1", meta$sampleNames)] <- "Syn1"

meta$`Cell_Type`[grep("HS01", meta$sampleNames)] <- "schwann"
meta$`Cell_Type`[grep("Syn5", meta$sampleNames)] <- "arachnoid"
meta$`Cell_Type`[grep("Syn6", meta$sampleNames)] <- "meningioma"
meta$`Cell_Type`[grep("HS11", meta$sampleNames)] <- "schwann"
meta$`Cell_Type`[grep("Syn1", meta$sampleNames)] <- "arachnoid"

meta$Treatment[grep("CUDC", meta$sampleNames)] <- "CUDC907"
meta$Treatment[grep("GSK", meta$sampleNames)] <- "GSK2126458"
meta$Treatment[grep("Pano", meta$sampleNames)] <- "Panobinostat"

rownames(meta) <- meta$sampleNames
meta <- select(meta, -sampleNames)

df<-AnnotatedDataFrame(meta)
kinome<-Biobase::ExpressionSet(assayData = kinome, phenoData = df)

saveRDS(kinome, "treated_kinome.rds")

this.file = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/kinomeAnalysis/kinomeDataforSyDE.R"
synStore(File("treated_kinome.rds", parentId = "syn10845573"), used = c("syn6181167"), executed = this.file)

