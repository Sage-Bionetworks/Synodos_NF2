library(synapseClient)
library(pheatmap)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
synapseLogin()

kin<-read.table(synGet("syn4951080")@filePath, sep = "\t", header = TRUE) %>% 
  filter(time == "24h", cellLine %in% c("HS01", "HS11", "Syn1", "Syn5", "Syn6")) %>% 
  unite(treatment, cellLine, drug, sep = "_") %>% 
  select(treatment, Gene, log2NormRatio) %>% 
  group_by(treatment, Gene) %>% 
  dplyr::summarize(mean = mean(log2NormRatio, na.rm = TRUE)) %>% 
  spread(Gene, mean) %>% 
  ungroup()

kin <- as.data.frame(kin)
rownames(kin) <- kin$treatment
kin <- kin[,-1]
kin <- t(kin)

##only considering complete cases works 
kin2<-kin[complete.cases(kin),]

pdf("kinome-mini-complete.pdf")
pheatmap(kin2, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 3)
dev.off()

pdf("kinome-zoom-complete.pdf", height = 15)
pheatmap(kin2, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 5,
         fontsize_col = 5, cellwidth = 5, cellheight = 5, height = 11)
dev.off()

##more than 8 NAs/case and hclust spits out error 
kin2<-kin[rowSums(is.na(kin))<8,]
pdf("kinome-mini-half.pdf", height = 15)
pheatmap(kin2, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 3)
dev.off()

pdf("kinome-zoom-half.pdf", height = 20)
pheatmap(kin2, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 5,
         fontsize_col = 5, cellwidth = 5, cellheight = 5, height = 11)
dev.off()

##setting to 0 doesn't produce good results
kin3 <- kin
kin3[is.na(kin3)] <- 0
pheatmap(kin3, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 3)

##filtering out non variant samples or 0 sd doesnt work either
kin3<-kin[!is.na(apply(kin, 1, function(x) sd(x, na.rm = TRUE))),]
kin3<-kin3[(apply(kin3, 1, function(x) sd(x, na.rm = TRUE))) != 0,]
pheatmap(kin3, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 3)

kin3<-kin[!is.na(apply(kin, 1, function(x) var(x, na.rm = TRUE))),]
pheatmap(kin3, color = viridis(1000, direction = 1, option = "D"), border_color = NA, fontsize_row = 3)

hist(as.matrix(kin))



##Get intersecting kinases

dekin<-read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE) %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.1) 

HS01.HS11.PANO <- dekin %>% filter(drug == "Pano", cellLine1 =="HS01", cellLine2 =="HS11")
HS01.HS11.GSK <- dekin %>% filter(drug == "GSK458", cellLine1 =="HS01", cellLine2 =="HS11")
HS01.HS11.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="HS01", cellLine2 =="HS11")
prots.hs<-intersect(intersect(HS01.HS11.CUDC$protein, HS01.HS11.GSK$protein), HS01.HS11.PANO$protein)

Syn5.Syn1.PANO <- dekin %>% filter(drug == "Pano", cellLine1 =="Syn5", cellLine2 =="Syn1")
Syn5.Syn1.GSK <- dekin %>% filter(drug == "GSK458", cellLine1 =="Syn5", cellLine2 =="Syn1")
Syn5.Syn1.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="Syn5", cellLine2 =="Syn1")

prots.syn<-intersect(intersect(Syn5.Syn1.CUDC$protein, Syn5.Syn1.GSK$protein), Syn5.Syn1.PANO$protein)
prots.all<-intersect(prots.syn,prots.hs)



hs01 <- dekin %>% 
  filter(drug %in% c("Pano", "GSK458", "CUDC"), 
         cellLine1 =="HS01", cellLine2 =="HS11", protein %in% prots.all) %>% 
  select(drug, protein, cellLine1, avgRatio_proteinRep_cond1)

hs11 <- dekin %>% 
  filter(drug %in% c("Pano", "GSK458", "CUDC"), 
         cellLine1 =="HS01", cellLine2 =="HS11", protein %in% prots.all) %>% 
  select(drug, protein, cellLine2, avgRatio_proteinRep_cond2)
colnames(hs11) <- c("drug", "protein", "cellLine1", "avgRatio_proteinRep_cond1")

hs <- rbind(hs01, hs11)

ggplot(data = hs, aes(x=interaction(protein, cellLine1, drug) %>% factor(levels = interaction(protein, cellLine1, drug)[order(avgRatio_proteinRep_cond1)]), 
                      y=avgRatio_proteinRep_cond1, 
                      fill = cellLine1)) +
  geom_bar(stat = "identity") +
  coord_flip()

