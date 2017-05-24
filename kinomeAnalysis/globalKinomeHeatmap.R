library(synapseClient)
library(pheatmap)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(miscTools)
library(matrixStats)
synapseLogin()

kin<-read.table(synGet("syn4951080")@filePath, sep = "\t", header = TRUE) %>% 
  filter(time == "24h", cellLine %in% c("HS01", "HS11", "Syn1", "Syn5", "Syn6")) %>% 
  unite(treatment, cellLine, drug, sep = "_") %>% 
  select(treatment, Gene, log2ratio) %>% 
  group_by(treatment, Gene) %>% 
  dplyr::summarize(mean = mean(log2ratio, na.rm = TRUE)) %>% 
  spread(Gene, mean) %>% 
  ungroup()

kin <- as.data.frame(kin)
rownames(kin) <- kin$treatment
kin <- kin[,-1]
kin <- t(kin)

##only considering complete cases works 
kin2<-kin[complete.cases(kin),]

medians <- colMedians(kin2)
sd <- colSds(kin2)
kin.center <- scale(kin2, center = medians, scale = FALSE)
kin.scale <- scale(kin2, center = medians, scale = sd)

pdf("kinome-mini-complete.pdf")
pheatmap(kin.scale, border_color = NA, fontsize_row = 3)
dev.off()

pdf("kinome-zoom-complete.pdf", height = 15)
pheatmap(kin.scale, border_color = NA, fontsize_row = 5,
         fontsize_col = 5, cellwidth = 5, cellheight = 5, height = 11)
dev.off()

##more than 8 NAs/case and hclust spits out error 
kin2<-kin
medians <- colMedians(kin2, na.rm = TRUE)
sd <- colSds(kin2, na.rm = TRUE)
kin.center <- scale(kin2, center = medians, scale = FALSE)
kin.scale <- scale(kin2, center = medians, scale = sd)

kin2<-kin2[rowSums(is.na(kin))<8,]


pdf("kinome-mini-half.pdf", height = 15)
pheatmap(kin2, border_color = NA, fontsize_row = 3, cellheight = 1.5, na_col = "grey", fontsize_col = 15)
dev.off()

pdf("kinome-zoom-half.pdf", height = 20)
pheatmap(kin2, border_color = NA, fontsize_row = 5,
         fontsize_col = 5, cellwidth = 5, cellheight = 5, height = 11, na_col = "grey")
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

##Get intersecting kinases (Top Baseline Hits)

basekin<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "") 

HS01.HS11.base <- basekin %>% filter(cellLine=="HS01", referenceSample=="HS11") %>% 
  group_by(Gene) %>% 
  dplyr::summarize(mean.log2ratio = mean(log2ratio, na.rm = TRUE), sem = (sd(log2ratio, na.rm = TRUE)/(sqrt(length(log2ratio)))), comp = "HS01_HS11") %>% 
  filter(abs(sem)<abs(mean.log2ratio)) %>% 
  filter(abs(mean.log2ratio)>0.1) 

HS01.HS11.base$Gene <- reorder(HS01.HS11.base$Gene, HS01.HS11.base$mean.log2ratio)

ggplot(data = HS01.HS11.base, aes(x=Gene, 
                        y=mean.log2ratio, 
                        fill = comp, group = comp)) +
  geom_errorbar(aes(x=Gene, ymin=mean.log2ratio-sem, ymax = mean.log2ratio+sem, color = comp),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#00a9ff","#ff9000")) +
  scale_color_manual(values = c("#00a9ff","#ff9000")) +
  coord_flip() +
  labs(y = "log2ratio (mean)", x = "Gene", title = "HS01 vs HS11, baseline") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("HS01-HS11-baseline_waterfall.pdf")

Syn5.Syn1.base <- basekin %>% filter(cellLine=="Syn5", referenceSample=="Syn1") %>% 
  group_by(Gene) %>% 
  dplyr::summarize(mean.log2ratio = mean(log2ratio, na.rm = TRUE), sem = (sd(log2ratio, na.rm = TRUE)/(sqrt(length(log2ratio)))), comp = "Syn5_Syn1")  %>% 
  filter(abs(sem)<abs(mean.log2ratio)) %>% 
  filter(abs(mean.log2ratio)>0.1) 

Syn5.Syn1.base$Gene <- reorder(Syn5.Syn1.base$Gene, Syn5.Syn1.base$mean.log2ratio)

ggplot(data = Syn5.Syn1.base, aes(x=Gene, 
                                  y=mean.log2ratio, 
                                  fill = comp, group = comp)) +
  geom_errorbar(aes(x=Gene, ymin=mean.log2ratio-sem, ymax = mean.log2ratio+sem, color = comp),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#ff9000")) +
  scale_color_manual(values = c("#ff9000")) +
  coord_flip() +
  labs(y = "log2ratio (mean)", x = "Gene", title = "Syn5 vs Syn1, baseline") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("Syn5-Syn1-baseline_waterfall.pdf")

common<-intersect(HS01.HS11.base$Gene, Syn5.Syn1.base$Gene)

base <- rbind(HS01.HS11.base, Syn5.Syn1.base) %>% filter(Gene %in% common)
base$Gene <- reorder(base$Gene, base$mean.log2ratio)

ggplot(data = base, aes(x=Gene, 
                      y=mean.log2ratio, 
                      fill = comp, group = comp)) +
  geom_errorbar(aes(x=Gene, ymin=mean.log2ratio-sem, ymax = mean.log2ratio+sem, color = comp),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#00a9ff","#ff9000"), name = "Comparison") +
  scale_color_manual(values = c("#00a9ff","#ff9000"), name = "Comparison") +
  coord_flip() +
  labs(y = "log2ratio (mean)", x = "Gene", title = "Common Kinases - HS01 vs HS11 and Syn5 vs Syn1") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Men-Schwann-common-baseline_waterfall.pdf")

##Get intersecting kinases (DE after treatment)

dekin<-read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE) %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.1)

HS01.HS11.PANO <- dekin %>% filter(drug == "Pano", cellLine1 =="HS01", cellLine2 =="HS11")
HS01.HS11.GSK <- dekin %>% filter(drug == "GSK458", cellLine1 =="HS01", cellLine2 =="HS11")
HS01.HS11.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="HS01", cellLine2 =="HS11")

Syn5.Syn1.PANO <- dekin %>% filter(drug == "Pano", cellLine1 =="Syn5", cellLine2 =="Syn1")
Syn5.Syn1.GSK <- dekin %>% filter(drug == "GSK458", cellLine1 =="Syn5", cellLine2 =="Syn1")
Syn5.Syn1.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="Syn5", cellLine2 =="Syn1")

prots.pano <-intersect(HS01.HS11.PANO$protein, Syn5.Syn1.PANO$protein)
prots.gsk <-intersect(HS01.HS11.GSK$protein, Syn5.Syn1.GSK$protein)
prots.cudc <-intersect(HS01.HS11.CUDC$protein, Syn5.Syn1.CUDC$protein)

dekin.pano <- filter(dekin, protein %in% prots.pano & drug == "Pano", cellLine1 %in% c("HS01", "Syn5") ) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1/avgRatio_proteinRep_cond2, comp = paste(cellLine1, "_", cellLine2, sep = "")) 
dekin.gsk <- filter(dekin, protein %in% prots.gsk & drug == "CUDC", cellLine1 %in% c("HS01", "Syn5") ) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1/avgRatio_proteinRep_cond2, comp = paste(cellLine1, "_", cellLine2, sep = "")) 
dekin.cudc <- filter(dekin, protein %in% prots.cudc & drug == "GSK458", cellLine1 %in% c("HS01", "Syn5") ) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1/avgRatio_proteinRep_cond2, comp = paste(cellLine1, "_", cellLine2, sep = "")) 

dekin.pano$protein <- reorder(dekin.pano$protein, dekin.pano$avgRatio)

ggplot(data = dekin.pano, aes(x=protein, 
                        y=avgRatio, 
                        fill = comp, group = comp)) +
  #geom_errorbar(aes(x=protein, ymin=avgRatio_proteinRep_cond1-se, ymax = avgRatio_proteinRep_cond1+se, color = cellLine1),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#DD1C1A", "#086788")) +
  scale_color_manual(values = c("#DD1C1A","#086788")) +
  labs(y = "NF1- treated/NF1+ treated ratio", x = "Gene", title = "Panobinostat Common Kinases") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_grid(. ~ comp) +
  coord_flip()

dekin.gsk$protein <- reorder(dekin.gsk$protein, dekin.gsk$avgRatio)

ggplot(data = dekin.gsk, aes(x=protein, 
                              y=avgRatio, 
                              fill = comp, group = comp)) +
  #geom_errorbar(aes(x=protein, ymin=avgRatio_proteinRep_cond1-se, ymax = avgRatio_proteinRep_cond1+se, color = cellLine1),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#DD1C1A", "#086788")) +
  scale_color_manual(values = c("#DD1C1A","#086788")) +
  labs(y = "NF1- treated/NF1+ treated ratio", x = "Gene", title = "GSK2126458 Common Kinases") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

dekin.cudc$protein <- reorder(dekin.cudc$protein, dekin.cudc$avgRatio)

ggplot(data = dekin.cudc, aes(x=protein, 
                             y=avgRatio, 
                             fill = comp, group = comp)) +
  #geom_errorbar(aes(x=protein, ymin=avgRatio_proteinRep_cond1-se, ymax = avgRatio_proteinRep_cond1+se, color = cellLine1),stat = "identity", position = "dodge") +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#DD1C1A", "#086788")) +
  scale_color_manual(values = c("#DD1C1A","#086788")) +
  labs(y = "NF1- treated/NF1+ treated ratio", x = "Gene", title = "CUDC907 Common Kinases") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_flip()

