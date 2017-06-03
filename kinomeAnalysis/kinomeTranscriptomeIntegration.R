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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Baseline") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HS01-HS11-baseline-integrated.png", height = 6, width = 6)

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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "CUDC907") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HS01-HS11-cudc-integrated.png", height = 6, width = 6)

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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "GSK2126458") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HS01-HS11-gsk-integrated.png", height = 6, width = 6)

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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Panobinostat") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HS01-HS11-panobinostat-integrated.png", height = 6, width = 6)

ggplot(bar, aes(y = logFC, x = FC)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(data = bar %>% filter(BH<0 & pval_adj <0)) +
  geom_label_repel(data = bar %>% filter(BH<0 & pval_adj <0), 
                   aes(label = Hugo_Gene, fill = BH<0.1), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c('TRUE' = '#3FA34D', 'FALSE' = '#A09F9D'), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Integrated Analysis") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("blankplot-integrated.png", height = 6, width = 6)


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

ggplot(bar, aes(y = logFC, x = FC)) +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Baseline") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Syn5-Syn1-baseline-integrated.png", height = 6, width = 6)


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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Baseline") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Syn5-Syn1-CUDC-integrated.png", height = 6, width = 6)

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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "GSK2126458") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Syn5-Syn1-GSK-integrated.png", height = 6, width = 6)


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
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(BH<0.05 & pval_adj <0.05), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines")) +
  scale_fill_manual(values = c("1" = "#F55D3E", "2" = "#F6AE2D", "3" = "#67AFEA", "4" = "#F6AE2D"), guide = "none") +
  labs(x = "FC Kinome", y = "log(2)FC Transcriptome", title = "Panobinostat") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Syn5-Syn1-GSK-integrated.png", height = 6, width = 6)