library(dplyr)
library(synapseClient)
library(tidyr)
synapseLogin()

all <- read.table(synGet("syn7437782")@filePath, sep = "\t", header = TRUE, na.strings = c("", "NA")) 

gsk <- all %>% 
  filter(treatment1 == "GSK458" & treatment2 == "GSK458") %>% 
  select(ensemblId, geneName, logFC, logCPM, BH)

colnames(gsk)[2:5] <- c("mouseGene","logFC.GSK", "logCPM.GSK", "BH.GSK")

pano <- all %>% 
  filter(treatment1 == "Pano" & treatment2 == "Pano") %>% 
  select(ensemblId, geneName, logFC, logCPM, BH)

colnames(pano)[2:5] <- c("mouseGene","logFC.Pano", "logCPM.Pano", "BH.Pano")

cudc <- all %>% 
  filter(treatment1 == "CUDC" & treatment2 == "CUDC")%>% 
  select(ensemblId, geneName, logFC, logCPM, BH)

colnames(cudc)[2:5] <- c("mouseGene","logFC.CUDC", "logCPM.CUDC",  "BH.CUDC")

all2 <- cudc %>% full_join(gsk) %>% full_join(pano)

homs <- read.table("HOM_MouseHumanSequence.rpt.txt", sep = "\t", header = TRUE) %>% 
  select(HomoloGene.ID, Common.Organism.Name, Symbol)

homs$Common.Organism.Name <- gsub("mouse, laboratory", "mouseGene", homs$Common.Organism.Name)
homs$Common.Organism.Name <- gsub("human", "humanGene", homs$Common.Organism.Name)

homs <- homs %>%
  group_by(Common.Organism.Name, HomoloGene.ID) %>%
  mutate(ind = row_number()) %>%
  spread(Common.Organism.Name, Symbol) %>% 
  ungroup() %>% 
  select(humanGene, mouseGene) %>% 
  filter(!is.na(mouseGene))

all2<-left_join(all2, homs)

dmso <- all %>% 
  filter(treatment1 == "DMSO" & treatment2 == "DMSO") %>% 
  select(ensemblId, geneName, logFC, logCPM, BH)

colnames(dmso)[2:5] <- c("mouseGene","logFC.DMSO", "logCPM.DMSO", "BH.DMSO")

dmso <- left_join(dmso, homs)

write.table(all2, "MS03vsMS12-treated.txt", sep = "\t")
write.table(dmso, "MS03vsMS12-DMSO.txt", sep = "\t")

#opened in excel and finished for collaborators there. Supplemental Table 7 for paper,
