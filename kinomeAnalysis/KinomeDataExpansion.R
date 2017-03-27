library(synapseClient)
library(dplyr)
synapseLogin()

dat<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "")
neighbors<-read.table(synGet("syn8517064")@filePath, sep = "\t", header = TRUE, comment.char = "")
colnames(neighbors) <- c("Gene", "Neighbor", "Confidence")

degenes<-read.table(synGet("syn6038243")@filePath, sep = "\t", header = TRUE, comment.char = "")

meningioma.kin.up<-dat %>% 
  filter(cellLine=="Syn5", referenceSample=="Syn1") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  filter(`mean(log2ratio)`>=0.1) %>% 
  arrange(desc(`mean(log2ratio)`))

men.kin.up.neigh.1 <- meningioma.kin.up %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

men.kin.up.neigh.2 <- men.kin.up.neigh.1 %>% 
  select(Neighbor, Confidence)
colnames(men.kin.up.neigh.2) <- c("Gene", "origConf")
men.kin.up.neigh.2 <- men.kin.up.neigh.2 %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

men.kin.up.2nd<-unique(c(as.character(meningioma.kin.up$Gene), as.character(men.kin.up.neigh.1$Neighbor), 
                     as.character(men.kin.up.neigh.2$Neighbor)))

men.upgenes<-degenes %>% filter(diffExptest=="Syn5.DMSO-Syn1.DMSO" & adj.P.Val<=0.1 & logFC>=0)
base::intersect(men.upgenes$geneName, men.kin.up.2nd)

#################################################################################################
meningioma.kin.down<-dat %>% 
  filter(cellLine=="Syn5", referenceSample=="Syn1") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  filter(`mean(log2ratio)`<=-0.1) %>% 
  arrange(desc(`mean(log2ratio)`))

men.kin.down.neigh.1 <- meningioma.kin.down %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

men.kin.down.neigh.2 <- men.kin.down.neigh.1 %>% 
  select(Neighbor, Confidence)
colnames(men.kin.down.neigh.2) <- c("Gene", "origConf")
men.kin.down.neigh.2 <- men.kin.down.neigh.2 %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

men.kin.down.2nd<-unique(c(as.character(meningioma.kin.down$Gene), as.character(men.kin.down.neigh.1$Neighbor), 
                         as.character(men.kin.down.neigh.2$Neighbor)))

men.downgenes<-degenes %>% filter(diffExptest=="Syn5.DMSO-Syn1.DMSO" & adj.P.Val<=0.1 & logFC<=0)
base::intersect(men.downgenes$geneName, men.kin.down.2nd)


##################################
base::intersect(men.downgenes$geneName, men.kin.up.2nd)
base::intersect(men.upgenes$geneName, men.kin.down.2nd)


##################################
schwannoma.kin.up<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  filter(`mean(log2ratio)`>=0.1) %>% 
  arrange(desc(`mean(log2ratio)`))

sch.kin.up.neigh.1 <- schwannoma.kin.up %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

sch.kin.up.neigh.2 <- sch.kin.up.neigh.1 %>% 
  select(Neighbor, Confidence)
colnames(sch.kin.up.neigh.2) <- c("Gene", "origConf")
sch.kin.up.neigh.2 <- sch.kin.up.neigh.2 %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

sch.kin.up.2nd<-unique(c(as.character(schwannoma.kin.up$Gene), as.character(sch.kin.up.neigh.1$Neighbor), 
                         as.character(sch.kin.up.neigh.2$Neighbor)))

sch.upgenes<-degenes %>% filter(diffExptest=="HS01.DMSO-HS11.DMSO" & adj.P.Val<=0.1 & logFC>=0)
base::intersect(sch.upgenes$geneName, sch.kin.up.2nd)

#################################################################################################
schwannoma.kin.down<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  filter(`mean(log2ratio)`<=-0.1) %>% 
  arrange(desc(`mean(log2ratio)`))

sch.kin.down.neigh.1 <- schwannoma.kin.down %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

sch.kin.down.neigh.2 <- sch.kin.down.neigh.1 %>% 
  select(Neighbor, Confidence)
colnames(sch.kin.down.neigh.2) <- c("Gene", "origConf")
sch.kin.down.neigh.2 <- sch.kin.down.neigh.2 %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.8)

sch.kin.down.2nd<-unique(c(as.character(schwannoma.kin.down$Gene), as.character(sch.kin.down.neigh.1$Neighbor), 
                           as.character(sch.kin.down.neigh.2$Neighbor)))

sch.downgenes<-degenes %>% filter(diffExptest=="HS01.DMSO-HS11.DMSO" & adj.P.Val<=0.1 & logFC<=0)
base::intersect(sch.downgenes$geneName, sch.kin.down.2nd)

##################################
base::intersect(sch.downgenes$geneName, sch.kin.up.2nd)
base::intersect(sch.upgenes$geneName, sch.kin.down.2nd)

base::intersect(sch.kin.up.2nd, sch.kin.down.2nd)
