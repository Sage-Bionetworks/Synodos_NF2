library(dplyr)
library(synapseClient)
library(tidyr)
synapseLogin()

all <- read.table(synGet("syn7437782")@filePath, sep = "\t", header = TRUE, na.strings = c("", "NA"))
idx <- grep("ERCC.+", all$ensemblId)
all <- all[-idx,]

homs <- read.table("HOM_MouseHumanSequence.rpt.txt", sep = "\t", header = TRUE) %>% 
  select(HomoloGene.ID, Common.Organism.Name, Symbol)

homs$Common.Organism.Name <- gsub("mouse, laboratory", "mouseGene", homs$Common.Organism.Name)
homs$Common.Organism.Name <- gsub("human", "humanGene", homs$Common.Organism.Name)

homs <- homs %>%
  group_by(Common.Organism.Name, HomoloGene.ID) %>%
  dplyr::mutate(ind = row_number()) %>%
  spread(Common.Organism.Name, Symbol) %>% 
  ungroup() %>% 
  dplyr::select(humanGene, mouseGene) %>% 
  filter(!is.na(mouseGene))

some <- all %>% filter(cellLine1 == "MS03", cellLine2 == "MS12")

dmso.1 <- some %>% 
  filter(treatment1 == "DMSO" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, BH)

colnames(dmso.1)[2:4] <- c("mouseGene","logFC.DMSO", "BH.DMSO")

#gsk.1 <- some %>% 
#  filter(treatment1 == "GSK458" & treatment2 == "GSK458") %>% 
#  select(ensemblId, geneName, logFC, logCPM, BH)
#
#colnames(gsk)[2:4] <- c("mouseGene","logFC.GSK", "logCPM.GSK", "BH.GSK")
# 
# pano.1 <- some %>% 
#   filter(treatment1 == "Pano" & treatment2 == "Pano") %>% 
#   select(ensemblId, geneName, logFC, logCPM, BH)
# 
# colnames(pano)[2:4] <- c("mouseGene","logFC.Pano", "logCPM.Pano", "BH.Pano")
# 
# cudc.1 <- some %>% 
#   filter(treatment1 == "CUDC" & treatment2 == "CUDC")%>% 
#   select(ensemblId, geneName, logFC, logCPM, BH)
# 
# colnames(cudc)[2:4] <- c("mouseGene","logFC.CUDC", "logCPM.CUDC",  "BH.CUDC")
# 
# some2 <- dmso1 %>% full_join(cudc.1) %>% full_join(gsk.1) %>% full_join(pano.1)
# 
# 
# write.table(all2, "MS03vsMS12-treated.txt", sep = "\t")
# write.table(dmso, "MS03vsMS12-DMSO.txt", sep = "\t")

#opened in excel and finished for collaborators there. Supplemental Table 7 for paper,

##within cell line comparisions 

###ms03

some <- all %>% filter(cellLine1 == "MS03", cellLine2 == "MS03")
gsk.2 <- some %>% 
  filter(treatment1 == "GSK458" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, BH)

colnames(gsk.2)[2:4] <- c("mouseGene","logFC.GSK.2", "BH.GSK.2")

pano.2 <- some %>% 
  filter(treatment1 == "Pano" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, BH)

colnames(pano.2)[2:4] <- c("mouseGene","logFC.Pano.2", "BH.Pano.2")

cudc.2 <- some %>% 
  filter(treatment1 == "CUDC" & treatment2 == "DMSO")%>% 
  select(ensemblId, geneName, logFC, BH)

colnames(cudc.2)[2:4] <- c("mouseGene","logFC.CUDC.2",  "BH.CUDC.2")

some <- all %>% filter(cellLine1 == "MS12", cellLine2 == "MS12")

gsk.3 <- some %>% 
  filter(treatment1 == "GSK458" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, BH)

colnames(gsk.3)[2:4] <- c("mouseGene","logFC.GSK.3", "BH.GSK.3")

pano.3 <- some %>% 
  filter(treatment1 == "Pano" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, BH)

colnames(pano.3)[2:4] <- c("mouseGene","logFC.Pano.3", "BH.Pano.3")

cudc.3 <- some %>% 
  filter(treatment1 == "CUDC" & treatment2 == "DMSO")%>% 
  select(ensemblId, geneName, logFC, BH)

colnames(cudc.3)[2:4] <- c("mouseGene","logFC.CUDC.3",  "BH.CUDC.3")

some2 <- dmso.1 %>% full_join(cudc.2) %>% full_join(pano.2) %>% full_join(gsk.2) %>% 
  full_join(cudc.3) %>% full_join(pano.3) %>% full_join(gsk.3)

all2<-left_join(some2, homs)

write.table(all2, "MouseDEJoined.txt", row.names = F, quote = F, sep = "\t")

