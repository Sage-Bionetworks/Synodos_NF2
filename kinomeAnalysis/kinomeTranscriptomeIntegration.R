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
##2022 note: syn5840701 refers to a private file, the same file is available publicly at syn6182638
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
bar$labelcolor <- "0"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio>0.5] <- "1"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio< -0.5] <- "2"
bar$labelcolor[bar$logFC>-1 &  bar$logFC<1 & bar$Mean_Kinome_Ratio>0.5] <- "3"
bar$labelcolor[bar$logFC>-1 & bar$logFC<1 & bar$Mean_Kinome_Ratio< -0.5] <- "4"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio>-0.5 & bar$Mean_Kinome_Ratio<0.5] <- "5"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio>-0.5 & bar$Mean_Kinome_Ratio<0.5] <- "6"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio< -0.5] <- "7"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio>0.5] <- "8"
bar$label <- ""
bar$label[abs(bar$logFC) > 1 | abs(bar$Mean_Kinome_Ratio) > 0.5] <- bar$Hugo_Gene[abs(bar$logFC) > 1 | abs(bar$Mean_Kinome_Ratio) > 0.5]

ggplot(bar, aes(x = Mean_Kinome_Ratio, y = logFC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 30) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("HS01-HS11-baseline-integrated.png", height = 6.5, width = 6.5)

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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# #ggplot(bar, aes(y = ratio, x = FC)) +
# #  theme_bw() +
# #  geom_hline(aes(yintercept = 0)) + 
# #  geom_vline(aes(xintercept = 0)) +
# #  geom_point() +
# #  geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
# #                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                     aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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

colnames(mn.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(degenes, diffExptest=="Syn5.DMSO-Syn1.DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(mn.kin)
bar$labelcolor <- "0"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio>0.5] <- "1"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio< -0.5] <- "2"
bar$labelcolor[bar$logFC>-1 &  bar$logFC<1 & bar$Mean_Kinome_Ratio>0.5] <- "3"
bar$labelcolor[bar$logFC>-1 & bar$logFC<1 & bar$Mean_Kinome_Ratio< -0.5] <- "4"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio>-0.5 & bar$Mean_Kinome_Ratio<0.5] <- "5"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio>-0.5 & bar$Mean_Kinome_Ratio<0.5] <- "6"
bar$labelcolor[bar$logFC>1 & bar$Mean_Kinome_Ratio< -0.5] <- "7"
bar$labelcolor[bar$logFC< -1 & bar$Mean_Kinome_Ratio>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)
bar$label <- ""
bar$label[abs(bar$logFC) > 1 | abs(bar$Mean_Kinome_Ratio) > 0.5] <- bar$Hugo_Gene[abs(bar$logFC) > 1 | abs(bar$Mean_Kinome_Ratio) > 0.5]

ggplot(bar, aes(x = Mean_Kinome_Ratio, y = logFC)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 30) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Syn5-Syn1-baseline-integrated.png", height = 6.5, width = 6.5)


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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
# bar$labelcolor[bar$ratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$ratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$ratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$ratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = ratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(abs(FC) > 1 | abs(ratio) > 1), 
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
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
schwann.cpm <- read.table(synGet("syn10912550")@filePath, sep = "\t", header = T)

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
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj", "HS01_ratio", "HS01_ratio_se")

bar <- inner_join(CUDC, sch.kin.cudc)
bar$logratio <- log2(bar$ratio)

##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Schwannoma_CUDC_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_rect(xmin = 0.5, xmax = Inf, ymin = 1, ymax = Inf, fill = "#FF3A45", alpha = 0.02) +
  geom_rect(xmin = -Inf, xmax = -0.5, ymin = -Inf, ymax = -1, fill = "#4EC1F3", alpha = 0.005) +
  geom_rect(xmin = 0.5, xmax = Inf, ymin = -1, ymax = 1, fill = "#FC9790", alpha = 0.01) +
  geom_rect(xmin = -Inf, xmax = -0.5, ymin = -1, ymax = 1, fill = "#FFFF9F", alpha = 0.01) +
  geom_rect(xmin = -0.5, xmax = 0.5, ymin = 1, ymax = Inf, fill = "#E7FF7D", alpha = 0.03) +
  geom_rect(xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = -1, fill = "#BDBDFF", alpha = 0.01) +
  geom_rect(xmin = -Inf, xmax = -0.5, ymin = 1, ymax = Inf, fill = "#FFFF54", alpha = 0.01) +
  geom_rect(xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = -1, fill = "#C152E7", alpha = 0.01) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = HS01_ratio), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_gradient(low = "grey40", high = "white") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")
ggsave("Schwannoma_CUDC_Integrated_DMSOnorm_alternativecolors.png", height = 6.5, width = 6.5)

bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$HS01_ratio)])
bar2 <- bar %>% top_n(10,-HS01_ratio)
bar3 <- bar %>% top_n(10, HS01_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = HS01_ratio, fill = HS01_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio>0), aes(x = Hugo_Gene, ymax=HS01_ratio+HS01_ratio_se, ymin = HS01_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio<0), aes(x = Hugo_Gene, ymin=HS01_ratio-HS01_ratio_se, ymax = HS01_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")

  
ggsave("Schwannoma_CUDC_HS01_Top20Kinases.png", height = 6.5, width = 5.5)


sch.kin.gsk<-sch.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj", "HS01_ratio", "HS01_ratio_se")

bar <- inner_join(GSK, sch.kin.gsk)
bar$logratio <- log2(bar$ratio)
##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")


ggsave("Schwannoma_GSK_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$HS01_ratio)])
bar2 <- bar %>% top_n(10,-HS01_ratio)
bar3 <- bar %>% top_n(10, HS01_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = HS01_ratio, fill = HS01_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio>0), aes(x = Hugo_Gene, ymax=HS01_ratio+HS01_ratio_se, ymin = HS01_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio<0), aes(x = Hugo_Gene, ymin=HS01_ratio-HS01_ratio_se, ymax = HS01_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")


ggsave("Schwannoma_GSK_HS01_Top20Kinases.png", height = 6.5, width = 5.5)


sch.kin.pano<-sch.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(sch.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj", "HS01_ratio", "HS01_ratio_se")

bar <- inner_join(PANO, sch.kin.pano)
bar$logratio <- log2(bar$ratio)
##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Schwannoma_PANO_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$HS01_ratio)])
bar2 <- bar %>% top_n(10,-HS01_ratio)
bar3 <- bar %>% top_n(10, HS01_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = HS01_ratio, fill = HS01_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio>0), aes(x = Hugo_Gene, ymax=HS01_ratio+HS01_ratio_se, ymin = HS01_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(HS01_ratio<0), aes(x = Hugo_Gene, ymin=HS01_ratio-HS01_ratio_se, ymax = HS01_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")

ggsave("Schwannoma_Pano_HS01_Top20Kinases.png", height = 6.5, width = 5.5)


##meningioma
mening.cpm <- read.table(synGet("syn10912525")@filePath, sep = "\t", header = T)

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
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj", "Syn5_ratio", "Syn5_ratio_se")

bar <- inner_join(CUDC, mn.kin.cudc)
bar$logratio <- log2(bar$ratio)
##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Meningioma_CUDC_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$Syn5_ratio)])
bar2 <- bar %>% top_n(10,-Syn5_ratio)
bar3 <- bar %>% top_n(10, Syn5_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = Syn5_ratio, fill = Syn5_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio>0), aes(x = Hugo_Gene, ymax=Syn5_ratio+Syn5_ratio_se, ymin = Syn5_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio<0), aes(x = Hugo_Gene, ymin=Syn5_ratio-Syn5_ratio_se, ymax = Syn5_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")


ggsave("Meningioma_CUDC_Syn5_Top20Kinases.png", height = 6.5, width = 5.5)


mn.kin.gsk<-mn.kin.tx %>% 
  filter(drug == "GSK458") %>%
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.gsk) <- c("Hugo_Gene", "FC", "pval_adj", "Syn5_ratio", "Syn5_ratio_se")

bar <- inner_join(GSK, mn.kin.gsk)
bar$logratio <- log2(bar$ratio)
##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Meningioma_GSK_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)


bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$Syn5_ratio)])
bar2 <- bar %>% top_n(10,-Syn5_ratio)
bar3 <- bar %>% top_n(10, Syn5_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = Syn5_ratio, fill = Syn5_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio>0), aes(x = Hugo_Gene, ymax=Syn5_ratio+Syn5_ratio_se, ymin = Syn5_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio<0), aes(x = Hugo_Gene, ymin=Syn5_ratio-Syn5_ratio_se, ymax = Syn5_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")


ggsave("Meningioma_GSK_Syn5_Top20Kinases.png", height = 6.5, width = 5.5)

mn.kin.pano<-mn.kin.tx %>% 
  filter(drug == "Pano") %>%
  select(protein, FC, pval_adj, medRatio_peptides_cond1, se_peptides_cond1) #%>% 
#group_by(protein) %>% 
#dplyr::summarise(mean(FC), mean(pval_adj)) %>%
#arrange(desc(`mean(FC)`))

colnames(mn.kin.pano) <- c("Hugo_Gene", "FC", "pval_adj", "Syn5_ratio", "Syn5_ratio_se")

bar <- inner_join(PANO, mn.kin.pano)
bar$logratio <- log2(bar$ratio)
##logratio - logFCtranscriptome, FC - logFCkinome
bar$labelcolor <- "0"
bar$labelcolor[bar$logratio>1 & bar$FC>0.5] <- "1"
bar$labelcolor[bar$logratio< -1 & bar$FC< -0.5] <- "2"
bar$labelcolor[bar$logratio>-1 &  bar$logratio<1 & bar$FC>0.5] <- "3"
bar$labelcolor[bar$logratio>-1 & bar$logratio<1 & bar$FC< -0.5] <- "4"
bar$labelcolor[bar$logratio>1 & bar$FC>-0.5 & bar$FC<0.5] <- "5"
bar$labelcolor[bar$logratio< -1 & bar$FC>-0.5 & bar$FC<0.5] <- "6"
bar$labelcolor[bar$logratio>1 & bar$FC< -0.5] <- "7"
bar$labelcolor[bar$logratio< -1 & bar$FC>0.5] <- "8"
bar$labelcolor <- as.factor(bar$labelcolor)

bar$label <- ""
bar$label[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)] <- 
  bar$Hugo_Gene[abs(bar$logratio) > 1 | abs(bar$FC) > 0.5 & !is.na(bar$FC) & !is.na(bar$logratio)]

ggplot(bar, aes(x = FC, y = logratio)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0)) +
  geom_point() +
  geom_label_repel(data = bar,
                   aes(label = label, fill = labelcolor), min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6,
                   force = 1) +
  xlab("Kinome") +
  ylab("Transcriptome") +
  scale_fill_manual(values = c("1" = "#FF3A45", "2" = "#4EC1F3", "3" = "#FC9790", "4" = "#FFFF9F",
                               "5" = "#E7FF7D", "6" = "#BDBDFF", "7" = "#FFFF54", "8" = "#C152E7"), guide = "none") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(-4,4)) +
  theme(legend.position="none")

ggsave("Meningioma_PANO_Integrated_DMSOnorm.png", height = 6.5, width = 6.5)

bar$Hugo_Gene <- factor(bar$Hugo_Gene, 
                        levels = bar$Hugo_Gene[order(bar$Syn5_ratio)])
bar2 <- bar %>% top_n(10,-Syn5_ratio)
bar3 <- bar %>% top_n(10, Syn5_ratio) %>% bind_rows(bar2)

ggplot(bar3, aes(x = logratio)) +
  theme_bw() +
  geom_bar(aes(x= Hugo_Gene, y = Syn5_ratio, fill = Syn5_ratio), stat = "identity", color = "black") +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio>0), aes(x = Hugo_Gene, ymax=Syn5_ratio+Syn5_ratio_se, ymin = Syn5_ratio),
                width=0.2) +
  geom_errorbar(data = bar3 %>% filter(Syn5_ratio<0), aes(x = Hugo_Gene, ymin=Syn5_ratio-Syn5_ratio_se, ymax = Syn5_ratio),
                width=0.2) +
  scale_fill_viridis_c(option = "C") +
  theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title = element_text(size = 0)) +
  ylab("Ratio to DMSO") +
  coord_flip(ylim = c(-2,2))+
  theme(legend.position="none")

ggsave("Meningioma_PANO_Syn5_Top20Kinases.png", height = 6.5, width = 5.5)





# 
# #### Comparison of new and old approaches
# 
# sch.kin.cudc<-sch.kin.tx %>% 
#   filter(drug == "CUDC") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")
# CUDC$Hugo_Gene <- rownames(CUDC)
# barHS11 <- filter(degenes, comparison %in% c("HS11CUDC907vsHS11DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# barHS01 <- filter(degenes, comparison %in% c("HS01CUDC907vsHS01DMSO")) %>% select(Hugo_Gene, logFC, BH, comparison) 
# bar <- inner_join(barHS01, barHS11, by = "Hugo_Gene") %>% mutate(ratio_manual = logFC.x-logFC.y) %>% inner_join(CUDC)
# 
# cor(bar$ratio_manual, log2(bar$ratio))
# 
# ggplot(bar, aes(x = ratio_manual, y = log2(ratio))) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   xlab("HS01CUDC907vsHS01DMSO/HS11CUDC907vsHS11DMSO") +
#   ylab("(HS01CUDC_CPM/HS01DMSO_CPM)/HS11CUDC_CPM/HS11DMSO_CPM") +
#   ggtitle("cor: 0.39" ) +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 8))
# 
# ggsave("test2.png", height = 6.5, width = 6.5)
# 
# 
# #HS01 - HS11 tx
# sch.kin.cudc<-sch.kin.tx %>% 
#   filter(drug == "CUDC") %>%
#   select(protein, FC, pval_adj) #%>% 
# #group_by(protein) %>% 
# #dplyr::summarise(mean(FC), mean(pval_adj)) %>%
# #arrange(desc(`mean(FC)`))
# 
# colnames(sch.kin.cudc) <- c("Hugo_Gene", "FC", "pval_adj")
# 
# bar <- inner_join(CUDC, sch.kin.cudc)
# bar$logratio <- log2(bar$ratio)
# bar$labelcolor[bar$logratio>0 & bar$FC>0] <- 1
# bar$labelcolor[bar$logratio>0 & bar$FC<0] <- 2
# bar$labelcolor[bar$logratio<0 & bar$FC<0] <- 3
# bar$labelcolor[bar$logratio<0 & bar$FC>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
# 
# ggplot(bar, aes(y = logratio, x = FC)) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 0)) + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_point() +
#   geom_label_repel(data = bar %>% filter(pval_adj <0.001), 
#                    aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
#                    box.padding = unit(0.35, "lines"),
#                    point.padding = unit(0.25, "lines"),
#                    size = 6) +
#   scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#   xlab("KINOME ratio") +
#   ylab("(HS01CUDC_CPM/HS01DMSO_CPM)/HS11CUDC_CPM/HS11DMSO_CPM") +
#   theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm"),
#         axis.title = element_text(size = 8))
# 
# ggsave("test3.png", height = 6.5, width = 6.5)
# 
# 


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
# bar$labelcolor[bar$logFC>0 & bar$Mean_Kinome_Ratio>0] <- 1
# bar$labelcolor[bar$logFC>0 & bar$Mean_Kinome_Ratio<0] <- 2
# bar$labelcolor[bar$logFC<0 & bar$Mean_Kinome_Ratio<0] <- 3
# bar$labelcolor[bar$logFC<0 & bar$Mean_Kinome_Ratio>0] <- 4
# bar$labelcolor <- as.factor(bar$labelcolor)
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
#                      aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
#                      box.padding = unit(0.35, "lines"),
#                      point.padding = unit(0.25, "lines"),
#                      size = 5) +
#     scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#     theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
#     labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
#   
#   
#   ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 6.5, width = 6.5)
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
#     #                 aes(label = Hugo_Gene, fill = labelcolor), min.segment.length = unit(0, "lines"),
#     #                 box.padding = unit(0.35, "lines"),
#     #                 point.padding = unit(0.25, "lines"),
#     #                 size = 5) +
#     scale_fill_manual(values = c("1" = "#FF9C99", "2" = "#FFC251", "3" = "#91CDFF", "4" = "#FFC251"), guide = "none") +
#     theme(axis.text = element_text(size = 18), plot.margin = unit(c(1,1,1,1), "cm")) +
#     labs(x= "Kinome Ratio", y = "log2FC Transcriptome", title = i)
#   
#   ggsave(paste("HS01-HS11-baseline-GO-",i,".png"), height = 6.5, width = 6.5)
# }
