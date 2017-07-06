library(synapseClient)
library(viridis)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
synapseLogin()

degenes <- read.table(synGet("syn9884855")@filePath, header = TRUE, sep = "\t")

dat<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "")
kinome <- read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)
sch.kin.tx <- kinome %>% 
  filter(cellLine1 == "HS01" & cellLine2 == "HS11" & time1 == "24h" & time2 == "24h") %>% 
  mutate(FC = medRatio_peptides_cond1-medRatio_peptides_cond2)

#HS01 - HS11 baseline
sch.kin<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  dplyr::summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(degenes, comparison=="HS01DMSOvsHS11DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin)
bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio>0] <- 1
bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio<0] <- 2
bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio<0] <- 3
bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = Mean_Kinome_Ratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 5) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-baseline-integrated.png", height = 8, width = 8)

ggplot(bar, aes(y = logFC, x = Mean_Kinome_Ratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(Mean_Kinome_Ratio)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-baseline-integrated_logFC.png", height = 8, width = 8)

#HS01 - HS11 tx
sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, comparison=="HS01CUDC907vsHS11CUDC907") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin.cudc)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
  axis.title = element_text(size = 0))

ggsave("HS01-HS11-cudc-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-cudc-integrated_logFC.png", height = 6.5, width = 6.5)

sch.kin.gsk<-sch.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, comparison=="HS01GSK2126458vsHS11GSK2126458") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin.gsk)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"), 
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-gsk-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75),                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"), 
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-gsk-integrated_logFC.png", height = 6.5, width = 6.5)


sch.kin.pano<-sch.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj) #%>% 
  #group_by(protein) %>% 
  #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
  #arrange(desc(`mean(FC)`))

colnames(sch.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, comparison=="HS01panobinostatvsHS11panobinostat") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin.pano)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-panobinostat-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("HS01-HS11-panobinostat-integrated_logFC.png", height = 6.5, width = 6.5)


ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(data = bar %>% filter(BH<0 & pval_adj <0)) +
  geom_label_repel(data = bar %>% filter(BH<0 & pval_adj <0), 
                   aes(label = Hugo_Gene, fill = BH<0.1), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("blankplot-integrated.png", height = 6.5, width = 6.5)

##meningioma 
degenes<-read.table(synGet("syn6166525")@filePath, sep = "\t", header = TRUE, comment.char = "")
degenes$Hugo_Gene <- gsub("\\|.+$", "", degenes$geneName)

mn.kin.tx <- kinome %>% 
  filter(cellLine1 == "Syn5" & cellLine2 == "Syn1" & time1 == "24h" & time2 == "24h") %>% 
  mutate(FC = medRatio_peptides_cond1-medRatio_peptides_cond2)

#Syn5 - Syn1 baseline
mn.kin<-dat %>% 
  filter(cellLine=="Syn5", referenceSample=="Syn1") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  dplyr::summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(mn.kin) <- c("Hugo_Gene", "FC")

bar <- filter(degenes, diffExptest=="Syn5.DMSO-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(mn.kin)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

## removed axis labels to clean up quad 3 collision 
ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 5,
                   max.iter = 20000) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-baseline-integrated.png", height = 8, width = 8)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4,
                   max.iter = 20000) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-baseline-integrated_logFC.png", height = 8, width = 8)


#Syn5 - Syn1 treatment
mn.kin.cudc<-mn.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, diffExptest=="Syn5.CUDC907-Syn1.CUDC907") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(mn.kin.cudc)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-CUDC-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-CUDC-integrated_logFC.png", height = 6.5, width = 6.5)

mn.kin.gsk<-mn.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, diffExptest=="Syn5.GSK2126458-Syn1.GSK2126458") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(mn.kin.gsk)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-GSK-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-GSK-integrated_logFC.png", height = 6.5, width = 6.5)


mn.kin.pano<-mn.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- filter(degenes, diffExptest=="Syn5.Panobinostat-Syn1.Panobinostat") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(mn.kin.pano)
bar$quadrant[bar$logFC>0 & bar$FC>0] <- 1
bar$quadrant[bar$logFC>0 & bar$FC<0] <- 2
bar$quadrant[bar$logFC<0 & bar$FC<0] <- 3
bar$quadrant[bar$logFC<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-Pano-integrated.png", height = 6.5, width = 6.5)

ggplot(bar, aes(y = logFC, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.75 | abs(FC)>0.75), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 4) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0))

ggsave("Syn5-Syn1-Pano-integrated_logFC.png", height = 6.5, width = 6.5)







#HS01 - HS11 GO overlap plots for UCF
sch.kin<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  dplyr::summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(degenes, comparison=="HS01DMSOvsHS11DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin)
bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio>0] <- 1
bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio<0] <- 2
bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio<0] <- 3
bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

upgenes <- read.table("../RNASeq_analysis/HS01vsHS11_DMSO_GOup_genes.txt", header = TRUE, sep = "\t")
downgenes <- read.table("../RNASeq_analysis/HS01vsHS11_DMSO_GOdown_genes.txt", header = TRUE, sep = "\t")

for(i in colnames(upgenes)){
  print(i)
  int <- base::intersect(bar$Hugo_Gene, upgenes[,i])
  print(int)
  bar2 <- filter(bar, Hugo_Gene %in% int)

  ggplot(bar2, aes(y = logFC, x = Mean_Kinome_Ratio)) +
    theme_bw() +
    geom_hline(aes(yintercept = 0)) + 
    geom_vline(aes(xintercept = 0)) +
    geom_point() +
    geom_label_repel(data = bar2, 
                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.25, "lines"),
                    size = 5) +
    scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
    theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
    labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
  

ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 8, width = 8)
}

for(i in colnames(downgenes)){
  print(i)
  int <- base::intersect(bar$Hugo_Gene, downgenes[,i])
  print(int)
  bar2 <- filter(bar, Hugo_Gene %in% int)
  
  ggplot(bar2, aes(y = logFC, x = Mean_Kinome_Ratio)) +
    theme_bw() +
    geom_hline(aes(yintercept = 0)) + 
    geom_vline(aes(xintercept = 0)) +
    geom_point() +
    geom_label_repel(data = bar2, 
                     aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.25, "lines"),
                     size = 5) +
    scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
    theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
    labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
  
  ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 8, width = 8)
}
