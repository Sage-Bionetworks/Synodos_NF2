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

# #HS01 - HS11 tx -- DEPRECATED!!
# sch.kin.cudc<-sch.kin.tx %>% 
#   filter(drug == "CUDC") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# bar <- select(schwann.cpm, HS01) %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(sch.kin.cudc)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# #ggplot(bar, aes(y = ratio, x = FC)) +
# #  theme_bw() +
# #  geom_hline(aes(yintercept = 0)) + 
# #  geom_vline(aes(xintercept = 0)) +
# #  geom_point() +
# #  geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
# #                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
# #                    box.padding = unit(0.35, "lines"),
# #                    point.padding = unit(0.25, "lines"),
# #                    size = 6) +
# #  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
# #  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
# #  axis.title = element_text(size = 0))
# 
# #ggsave("HS01-HS11-cudc-integrated.png", height = 6.5, width = 6.5)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(pval_adj <0.001), 
#                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"),
#                    size = 6) +
#   xlab("KINOME ratio") +
#   ylab("HS01CUDC907vsHS01DMSO/HS11CUDC907vsHS11DMSO") +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 8))
# 
# ggsave("test1.png", height = 6.5, width = 6.5)
# 
# barHS11 <- filter(degenes, comparison %in% c("HS11CUDC907vsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01CUDC907vsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barOLD <- filter(degenes, comparison %in% c("HS01CUDC907vsHS11CUDC907")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("HS01-HS11-cudc-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)
# 
# 
# sch.kin.gsk<-sch.kin.tx %>% 
#   filter(drug == "GSK458") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# colnames(sch.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# barHS11 <- filter(degenes, comparison %in% c("HS11GSK2126458vsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01GSK2126458vsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(sch.kin.gsk)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                     aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                     box.padding = unit(0.35, "lines"),
#                     point.padding = unit(0.25, "lines"), 
#                     size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("HS01-HS11-gsk-integrated.png", height = 6.5, width = 6.5)
# 
# barHS11 <- filter(degenes, comparison %in% c("HS11GSK2126458vsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01GSK2126458vsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barOLD <- filter(degenes, comparison %in% c("HS01GSK2126458vsHS11GSK2126458")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("HS01-HS11-gsk-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)
# 
# sch.kin.pano<-sch.kin.tx %>% 
#   filter(drug == "Pano") %>%
#   select(protein, FC, pval_adj) #%>% 
#   #group_by(protein) %>% 
#   #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#   #arrange(desc(`mean(FC)`))
# 
# colnames(sch.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# barHS11 <- filter(degenes, comparison %in% c("HS11panobinostatvsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01panobinostatvsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(sch.kin.pano)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"), 
#                    size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("HS01-HS11-panobinostat-integrated.png", height = 6.5, width = 6.5)
# 
# barHS11 <- filter(degenes, comparison %in% c("HS11panobinostatvsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01panobinostatvsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barOLD <- filter(degenes, comparison %in% c("HS01panobinostatvsHS11panobinostat")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("HS01-HS11-panobinostat-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)

# 
# ggplot(bar, aes(y = logFC, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point(data = bar %>% filter(BH<0 & pval_adj <0)) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("blankplot-integrated.png", height = 6.5, width = 6.5)

##meningioma 
degenes<-read.table(synGet("syn6166525")@filePath, sep = "\t", header = TRUE, comment.char = "")
degenes$Hugo_Gene <- gsub("\\|.+$", "", degenes$geneName)

mn.kin.tx <- kinome %>% 
  filter(cellLine1 == "Syn5" & cellLine2 == "Syn1" & time1 == "24h" & time2 == "24h") %>% 
  mutate(FC = medRatio_peptides_cond1-medRatio_peptides_cond2)


##some p-vals cannot be calculated due to lack of counts, set those to 1 
mn.kin.tx$pval_adj[is.na(mn.kin.tx$pval_adj)] <- 1

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


# #Syn5 - Syn1 treatment
# mn.kin.cudc<-mn.kin.tx %>% 
#   filter(drug == "CUDC") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(mn.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.CUDC907-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.CUDC907-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(mn.kin.cudc)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"), 
#                    size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-CUDC-integrated.png", height = 6.5, width = 6.5)
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.CUDC907-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.CUDC907-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# barOLD <- filter(degenes, diffExptest %in% c("Syn5.CUDC907-Syn1.CUDC907")) %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-cudc-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)
# 
# mn.kin.gsk<-mn.kin.tx %>% 
#   filter(drug == "GSK458") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(mn.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.GSK2126458-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.GSK2126458-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(mn.kin.gsk)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"), 
#                    size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-GSK-integrated.png", height = 6.5, width = 6.5)
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.GSK2126458-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.GSK2126458-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# barOLD <- filter(degenes, diffExptest %in% c("Syn5.GSK2126458-Syn1.GSK2126458")) %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-gsk-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)
# 
# 
# mn.kin.pano<-mn.kin.tx %>% 
#   filter(drug == "Pano") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(mn.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.Panobinostat-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.Panobinostat-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(mn.kin.pano)
# 
# bar$quadrant[bar$ratio>0 & bar$FC>0] <- 1
# bar$quadrant[bar$ratio>0 & bar$FC<0] <- 2
# bar$quadrant[bar$ratio<0 & bar$FC<0] <- 3
# bar$quadrant[bar$ratio<0 & bar$FC>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"), 
#                    size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-Pano-integrated.png", height = 6.5, width = 6.5)
# 
# barSyn5 <- filter(degenes, diffExptest=="Syn5.Panobinostat-Syn5.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest)  
# barSyn1 <- filter(degenes, diffExptest=="Syn1.Panobinostat-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# barOLD <- filter(degenes, diffExptest %in% c("Syn5.Panobinostat-Syn1.Panobinostat")) %>% select(Hugo_Gene, logFC, BH, diffExptest) 
# bar <- inner_join(barSyn5, barSyn1, by = "Hugo_Gene") %>% mutate(ratio = logFC.x-logFC.y) %>% inner_join(barOLD, by ="Hugo_Gene")
# 
# ggplot(bar, aes(y = ratio, x = logFC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 0))
# 
# ggsave("Syn5-Syn1-pano-oldvsnewTranscriptome.png", height = 6.5, width = 6.5)
# 




#### NEW ANALYSIS using CPM() edgeR output - use this. 
schwann.cpm <- read.table(synGet("syn10595624")@filePath, sep = "\t", header = T)

schwann.cpm$Gene <- rownames(schwann.cpm)
schwann.cpm <- mutate(schwann.cpm, HS01_DMSO_average = (HS01_DMSO_Run1_S1+HS01_DMSO_Run2_S1+HS01_DMSO_Run4_S9)/3)
schwann.cpm <- mutate(schwann.cpm, HS11_DMSO_average = (HS11_DMSO_Run1_S5+HS11_DMSO_Run2_S5+HS11_DMSO_Run3_S5+HS11_DMSO_Run4_S9)/4)

kinome <- read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)
sch.kin.tx <- kinome %>% 
  filter(cellLine1 == "HS01" & cellLine2 == "HS11" & time1 == "24h" & time2 == "24h") %>% 
  mutate(FC = medRatio_peptides_cond1-medRatio_peptides_cond2)

##some p-vals cannot be calculated due to lack of counts, set those to 1 
sch.kin.tx$pval_adj[is.na(sch.kin.tx$pval_adj)] <- 1

schwann.cpm <- filter(schwann.cpm, Gene %in% sch.kin.tx$protein)

schwann.cpm <- mutate(schwann.cpm, HS01_CUDC907_Run1_S2_dmsocorr=(HS01_CUDC907_Run1_S2/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_CUDC907_Run2_S2_dmsocorr=(HS01_CUDC907_Run2_S2/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_CUDC907_Run3_S2_dmsocorr=(HS01_CUDC907_Run3_S2/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_GSK2126458_Run1_S3_dmsocorr=(HS01_GSK2126458_Run1_S3/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_GSK2126458_Run2_S3_dmsocorr=(HS01_GSK2126458_Run2_S3/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_GSK2126458_Run3_S3_dmsocorr=(HS01_GSK2126458_Run3_S3/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_panobinostat_Run1_S4_dmsocorr=(HS01_panobinostat_Run1_S4/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_panobinostat_Run2_S4_dmsocorr=(HS01_panobinostat_Run2_S4/HS01_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS01_panobinostat_Run3_S4_dmsocorr=(HS01_panobinostat_Run3_S4/HS01_DMSO_average))

schwann.cpm <- mutate(schwann.cpm, HS11_CUDC907_Run1_S6_dmsocorr=(HS11_CUDC907_Run1_S6/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_CUDC907_Run2_S6_dmsocorr=(HS11_CUDC907_Run2_S6/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_CUDC907_Run3_S6_dmsocorr=(HS11_CUDC907_Run3_S6/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_GSK2126458_Run1_S7_dmsocorr=(HS11_GSK2126458_Run1_S7/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_GSK2126458_Run2_S7_dmsocorr=(HS11_GSK2126458_Run2_S7/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_GSK2126458_Run3_S7_dmsocorr=(HS11_GSK2126458_Run3_S7/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_panobinostat_Run1_S8_dmsocorr=(HS11_panobinostat_Run1_S8/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_panobinostat_Run2_S8_dmsocorr=(HS11_panobinostat_Run2_S8/HS11_DMSO_average))
schwann.cpm <- mutate(schwann.cpm, HS11_panobinostat_Run3_S8_dmsocorr=(HS11_panobinostat_Run3_S8/HS11_DMSO_average))

row.names(schwann.cpm) <- schwann.cpm$Gene

CUDC <- as.data.frame(t(sapply(schwann.cpm$Gene, function(i){
  test <- t.test(schwann.cpm[i, 30:32], schwann.cpm[i, 39:41])
  ratio <- sum(schwann.cpm[i, 30:32])/sum(schwann.cpm[i, 39:41])
  c(test$p.value, ratio)
})))

colnames(CUDC) <- c("p", "ratio")
CUDC$BH <- p.adjust(CUDC$p, method = "BH")
CUDC$Hugo_Gene <- rownames(CUDC)

GSK <- as.data.frame(t(sapply(schwann.cpm$Gene, function(i){
  test <- t.test(schwann.cpm[i, 33:35], schwann.cpm[i, 42:44])
  ratio <- sum(schwann.cpm[i, 33:35])/sum(schwann.cpm[i, 42:44])
  c(test$p.value, ratio)
})))

colnames(GSK) <- c("p", "ratio")
GSK$BH <- p.adjust(GSK$p, method = "BH")
GSK$Hugo_Gene <- rownames(GSK)

PANO <- as.data.frame(t(sapply(schwann.cpm$Gene, function(i){
  test <- t.test(schwann.cpm[i, 36:38], schwann.cpm[i, 45:47])
  ratio <- sum(schwann.cpm[i, 36:38])/sum(schwann.cpm[i, 45:47])
  c(test$p.value, ratio)
})))
colnames(PANO) <- c("p", "ratio")
PANO$BH <- p.adjust(PANO$p, method = "BH")
PANO$Hugo_Gene <- rownames(PANO)

sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(CUDC, sch.kin.cudc)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Schwannoma_CUDC_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

sch.kin.gsk<-sch.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(GSK, sch.kin.gsk)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Schwannoma_GSK_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)


sch.kin.pano<-sch.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(PANO, sch.kin.pano)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Schwannoma_PANO_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)


mening.cpm <- read.table(synGet("syn10595790")@filePath, sep = "\t", header = T)

mening.cpm$Hugo_Gene <- sapply(rownames(mening.cpm), function(i) {
  split<-unlist(strsplit(i, split = "\\|"))
  split[1]
  })

kinome <- read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)

mn.kin.tx <- kinome %>% 
  filter(cellLine1 == "Syn5" & cellLine2 == "Syn1" & time1 == "24h" & time2 == "24h") %>% 
  mutate(FC = medRatio_peptides_cond1-medRatio_peptides_cond2)

mening.cpm<-filter(mening.cpm,Hugo_Gene %in% sch.kin.tx$protein) %>% 
  select(Syn1.1.DMSO,Syn1.2.DMSO,Syn1.3.DMSO,Syn1.1.CUDC907,Syn1.2.CUDC907,
         Syn1.3.CUDC907,Syn1.1.Panobinostat,Syn1.2.Panobinostat,Syn1.3.Panobinostat,Syn1.1.GSK2126458,
         Syn1.2.GSK2126459,Syn1.3.GSK2126460,Syn5.1.DMSO,Syn5.2.DMSO,Syn5.3.DMSO,
         Syn5.1.CUDC907,Syn5.2.CUDC907,Syn5.3.CUDC907,Syn5.1.Panobinostat,Syn5.2.Panobinostat,
         Syn5.3.Panobinostat,Syn5.1.GSK2126458,Syn5.2.GSK2126459,Syn5.3.GSK2126460,Hugo_Gene)

mening.cpm <- mutate(mening.cpm, Syn1_DMSO_average = (Syn1.1.DMSO+Syn1.2.DMSO+Syn1.3.DMSO)/3)
mening.cpm <- mutate(mening.cpm, Syn5_DMSO_average = (Syn5.1.DMSO+Syn5.2.DMSO+Syn5.3.DMSO)/3)

mening.cpm <- mutate(mening.cpm, Syn1.1.CUDC907_dmsocorr=(Syn1.1.CUDC907/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.2.CUDC907_dmsocorr=(Syn1.2.CUDC907/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.3.CUDC907_dmsocorr=(Syn1.3.CUDC907/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.1.GSK2126458_dmsocorr=(Syn1.1.GSK2126458/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.2.GSK2126459_dmsocorr=(Syn1.2.GSK2126459/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.3.GSK2126460_dmsocorr=(Syn1.3.GSK2126460/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.1.Panobinostat_dmsocorr=(Syn1.1.Panobinostat/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.2.Panobinostat_dmsocorr=(Syn1.2.Panobinostat/Syn1_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn1.3.Panobinostat_dmsocorr=(Syn1.3.Panobinostat/Syn1_DMSO_average))

mening.cpm <- mutate(mening.cpm, Syn5.1.CUDC907_dmsocorr=(Syn5.1.CUDC907/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.2.CUDC907_dmsocorr=(Syn5.2.CUDC907/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.3.CUDC907_dmsocorr=(Syn5.3.CUDC907/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.1.GSK2126458_dmsocorr=(Syn5.1.GSK2126458/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.2.GSK2126459_dmsocorr=(Syn5.2.GSK2126459/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.3.GSK2126460_dmsocorr=(Syn5.3.GSK2126460/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.1.Panobinostat_dmsocorr=(Syn5.1.Panobinostat/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.2.Panobinostat_dmsocorr=(Syn5.2.Panobinostat/Syn5_DMSO_average))
mening.cpm <- mutate(mening.cpm, Syn5.3.Panobinostat_dmsocorr=(Syn5.3.Panobinostat/Syn5_DMSO_average))

mening.cpm <- filter(mening.cpm) %>% group_by(Hugo_Gene) %>% summarize_all(funs(sum)) %>% ungroup()
row.names(mening.cpm) <- mening.cpm$Hugo_Gene
mening.cpm <- as.data.frame(mening.cpm)

CUDC <- as.data.frame(t(sapply(mening.cpm$Hugo_Gene, function(i){
  #test <- t.test(mening.cpm[i, 28:30], mening.cpm[i, 37:39], paired = FALSE)
  ##not enough variance to calculate p-values
  ratio <- sum(mening.cpm[i, 37:39])/sum(mening.cpm[i, 28:30])
  c(NA, ratio)
})))

colnames(CUDC) <- c("p", "ratio")
CUDC$BH <- p.adjust(CUDC$p, method = "BH")
CUDC$Hugo_Gene <- rownames(CUDC)

GSK <- as.data.frame(t(sapply(mening.cpm$Hugo_Gene, function(i){
  #test <- t.test(mening.cpm[i, 31:33], mening.cpm[i, 40:42])
  ##not enough variance to calculate p-values
  ratio <- sum(mening.cpm[i, 40:42])/sum(mening.cpm[i, 31:33])
  c(NA, ratio)
})))

colnames(GSK) <- c("p", "ratio")
GSK$BH <- p.adjust(GSK$p, method = "BH")
GSK$Hugo_Gene <- rownames(GSK)

PANO <- as.data.frame(t(sapply(mening.cpm$Hugo_Gene, function(i){
  #test <- t.test(mening.cpm[i, 34:36], mening.cpm[i, 43:45])
  ##not enough variance to calculate p-values
  ratio <- sum(mening.cpm[i, 43:45])/sum(mening.cpm[i, 34:36])
  c(NA, ratio)
})))
colnames(PANO) <- c("p", "ratio")
PANO$BH <- p.adjust(PANO$p, method = "BH")
PANO$Hugo_Gene <- rownames(PANO)

mn.kin.cudc<-mn.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(CUDC, mn.kin.cudc)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Meningioma_CUDC_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

mn.kin.gsk<-mn.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(GSK, mn.kin.gsk)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Meningioma_GSK_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)


mn.kin.pano<-mn.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(PANO, mn.kin.pano)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point(aes(color = pval_adj <0.05)) +
  geom_label_repel(data = bar %>% filter(abs(FC)>1 | abs(logratio) > 1), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFF791", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-3.5,3.5)) +
  theme(legend.position="none")

ggsave("Meningioma_PANO_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)








#### Comparison of new and old approaches

sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")
CUDC$Hugo_Gene <- rownames(CUDC)
barHS11 <- filter(degenes, comparison %in% c("HS11CUDC907vsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
barHS01 <- filter(degenes, comparison %in% c("HS01CUDC907vsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio_manual = logFC.x-logFC.y) %>% inner_join(CUDC)

cor(bar$ratio_manual, log2(bar$ratio))

ggplot(bar, aes(x = ratio_manual, y = log2(ratio))) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  xlab("HS01CUDC907vsHS01DMSO/HS11CUDC907vsHS11DMSO") +
  ylab("(HS01CUDC_CPM/HS01DMSO_CPM)/HS11CUDC_CPM/HS11DMSO_CPM") +
  ggtitle("cor: 0.39" ) +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8))

ggsave("test2.png", height = 6.5, width = 6.5)


#HS01 - HS11 tx
sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, FC, pval_adj) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")

bar <- inner_join(CUDC, sch.kin.cudc)
bar$logratio <- log2(bar$ratio)
bar$quadrant[bar$logratio>0 & bar$FC>0] <- 1
bar$quadrant[bar$logratio>0 & bar$FC<0] <- 2
bar$quadrant[bar$logratio<0 & bar$FC<0] <- 3
bar$quadrant[bar$logratio<0 & bar$FC>0] <- 4
bar$quadrant <- as.factor(bar$quadrant)

ggplot(bar, aes(y = logratio, x = FC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar %>% filter(pval_adj <0.001), 
                   aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
  xlab("KINOME ratio") +
  ylab("(HS01CUDC_CPM/HS01DMSO_CPM)/HS11CUDC_CPM/HS11DMSO_CPM") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 8))

ggsave("test3.png", height = 6.5, width = 6.5)




# #HS01 - HS11 GO overlap plots for UCF --- old stuff!
# sch.kin<-dat %>% 
#   filter(cellLine=="HS01", referenceSample=="HS11") %>%
#   select(Gene, log2ratio) %>% 
#   group_by(Gene) %>% 
#   dplyr::summarise(mean(log2ratio)) %>%
#   arrange(desc(`mean(log2ratio)`))
# 
# colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")
# 
# bar <- filter(degenes, comparison=="HS01DMSOvsHS11DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin)
# bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio>0] <- 1
# bar$quadrant[bar$logFC>0 & bar$Mean_Kinome_Ratio<0] <- 2
# bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio<0] <- 3
# bar$quadrant[bar$logFC<0 & bar$Mean_Kinome_Ratio>0] <- 4
# bar$quadrant <- as.factor(bar$quadrant)
# 
# upgenes <- read.table("../RNASeq_analysis/HS01vsHS11_DMSO_GOup_genes.txt", header = TRUE, sep = "\t")
# downgenes <- read.table("../RNASeq_analysis/HS01vsHS11_DMSO_GOdown_genes.txt", header = TRUE, sep = "\t")

# for(i in colnames(upgenes)){
#   print(i)
#   int <- base::intersect(bar$Hugo_Gene, upgenes[,i])
#   print(int)
#   bar2 <- filter(bar, Hugo_Gene %in% int)
#   
#   ggplot(bar2, aes(y = logFC, x = Mean_Kinome_Ratio)) +
#     theme_bw() +
#     geom_hline(aes(yintercept = 0)) + 
#     geom_vline(aes(xintercept = 0)) +
#     geom_point() +
#     geom_label_repel(data = bar2, 
#                      aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#                      box.padding = unit(0.35, "lines"),
#                      point.padding = unit(0.25, "lines"),
#                      size = 5) +
#     scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#     theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
#     labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
#   
#   
#   ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 8, width = 8)
# }

# for(i in colnames(downgenes)){
#   print(i)
#   int <- base::intersect(bar$Hugo_Gene, downgenes[,i])
#   print(int)
#   bar2 <- filter(bar, Hugo_Gene %in% int)
#   
#   ggplot(bar2, aes(y = logFC, x = Mean_Kinome_Ratio)) +
#     theme_bw() +
#     geom_hline(aes(yintercept = 0)) + 
#     geom_vline(aes(xintercept = 0)) +
#     geom_point() +
#     #geom_label_repel(data = bar2, 
#     #                 aes(label = Hugo_Gene, fill = quadrant), min.segment.length = unit(0, "lines"),
#     #                 box.padding = unit(0.35, "lines"),
#     #                 point.padding = unit(0.25, "lines"),
#     #                 size = 5) +
#     scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#     theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
#     labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
#   
#   ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 8, width = 8)
# }
