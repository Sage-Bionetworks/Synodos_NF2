library(synapseClient)
library(dplyr)
library(biomaRt)
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
  filter(Confidence>=0.6)

men.kin.up<-unique(c(as.character(meningioma.kin.up$Gene), as.character(men.kin.up.neigh.1$Neighbor)))

men.upgenes<-degenes %>% filter(diffExptest=="Syn5.DMSO-Syn1.DMSO" & adj.P.Val<=0.1 & logFC>=0)
#base::intersect(men.upgenes$geneName, men.kin.up)

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
  filter(Confidence>=0.6)

men.kin.down<-unique(c(as.character(meningioma.kin.down$Gene), as.character(men.kin.down.neigh.1$Neighbor)))

#men.downgenes<-degenes %>% filter(diffExptest=="Syn5.DMSO-Syn1.DMSO" & adj.P.Val<=0.1 & logFC<=0)
#base::intersect(men.downgenes$geneName, men.kin.down.2nd)


##################################
#base::intersect(men.downgenes$geneName, men.kin.up.2nd)
#base::intersect(men.upgenes$geneName, men.kin.down.2nd)

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
  filter(Confidence>=0.6)

sch.kin.up<-as.data.frame(unique(c(as.character(schwannoma.kin.up$Gene), as.character(sch.kin.up.neigh.1$Neighbor))))

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
  filter(Confidence>=0.6)

sch.kin.down<-as.data.frame(unique(c(as.character(schwannoma.kin.down$Gene), as.character(sch.kin.down.neigh.1$Neighbor))))

sch.downgenes<-degenes %>% filter(diffExptest=="MS03.GSK212458-MS11.GSK212458" & adj.P.Val<=0.1 & logFC>=0)
base::intersect(sch.downgenes$geneName, sch.kin.down)

##################################
#base::intersect(sch.downgenes$geneName, sch.kin.up.2nd)
#base::intersect(sch.upgenes$geneName, sch.kin.down.2nd)

write.table(sch.kin.down, "Sch_Kin_up_expanded.txt",col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(sch.kin.up, "Sch_Kin_down_expanded.txt",col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(schwannoma.kin.up$Gene, "Sch_Kin_Up.txt",col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(schwannoma.kin.down$Gene, "Sch_Kin_Down.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")



#################################################################################################
dat2<-read.table(synGet("syn5840664")@filePath, sep = "\t", header = TRUE, comment.char = "") %>% 
  mutate(ratio=log2(MS03/MS12)) %>% select(Gene, ratio) %>%  
  filter(!is.na(ratio) & ratio!=Inf & ratio!=-Inf)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
mousemap<-getLDS(attributes = c("hgnc_symbol"), mart = human,
       attributesL = c("mgi_symbol"), martL = mouse)

mousemap$MGI.symbol<-toupper(mousemap$MGI.symbol)
colnames(mousemap) <- c("hugo", "Gene")
dat2<-left_join(dat2, mousemap)

mschwannoma.kin.up<-dat2 %>% 
  filter(ratio>=0.1) %>% 
  arrange(desc(ratio)) %>% 
  dplyr::select(hugo)

colnames(mschwannoma.kin.up)<-"Gene"

msch.kin.up.neigh.1 <- mschwannoma.kin.up %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.6)

msch.kin.up<-as.data.frame(unique(c(as.character(mschwannoma.kin.up$Gene), as.character(msch.kin.up.neigh.1$Neighbor))))

write.table(mschwannoma.kin.up, "ms_Kin_up.txt",col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(msch.kin.up, "ms_Kin_up_expanded.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#################################################################################################

mschwannoma.kin.down<-dat2 %>% 
  filter(ratio<=-0.1) %>% 
  arrange(desc(ratio)) %>% 
  dplyr::select(hugo)

colnames(mschwannoma.kin.down)<-"Gene"

msch.kin.down.neigh.1 <- mschwannoma.kin.down %>% 
  left_join(neighbors) %>% 
  filter(Confidence>=0.6)

msch.kin.down<-as.data.frame(unique(c(as.character(mschwannoma.kin.down$Gene), as.character(msch.kin.down.neigh.1$Neighbor))))

write.table(mschwannoma.kin.down, "ms_Kin_down.txt",col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(msch.kin.down, "ms_Kin_down_expanded.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

