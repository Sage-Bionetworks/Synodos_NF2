library(synapseClient)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggrepel)
files<-dir(path = "../UNC_edgeR/degenes")

degenes<-lapply(files, function(x){
  foo<-read.table(paste("../UNC_edgeR/degenes/",x, sep = ""), header = T)
  comp <- gsub("_edgeR_quasilikelihoodFtest.txt", "", x)
  foo$comparison <- c(rep(comp, nrow(foo)))
  foo$Hugo_Gene <- rownames(foo)
  foo
})

degenes <- ldply(degenes)

for(x in unique(degenes$comparison)){
  print(x)
  bar <- filter(degenes, comparison==x)
  bar2 <- filter(bar, BH < 0.1) %>% group_by(Hugo_Gene) %>% top_n(1, -logFC) %>% ungroup()
  
  p<-ggplot(bar2, aes(x=Hugo_Gene %>% reorder(logFC) , y=logFC, fill = logFC)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 8)) +
    labs(x = paste("Gene (n = ", nrow(bar2),")", sep =""), y = "log2(FC), BH<0.1", title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots/",x,".png", sep = ""), height = 8, width = 5)
  
  bar3 <- distinct(rbind((bar2 %>% top_n(15, logFC)), (bar2 %>% top_n(15, -logFC))))
  
  p<-ggplot(bar3, aes(x=Hugo_Gene %>% reorder(logFC), y=logFC, fill = logFC)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots/",x,"top30.png", sep = ""), height = 8, width = 5)
  
  #synStore(File(paste("UNCDEgenePlots/",x,".png", sep = ""), parentId = "syn9779338"), used = syns, executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCSchwannRNASeq.R")
  #synStore(File(paste("UNCDEgenePlots/",x,"top30.png", sep = ""), parentId = "syn9779338"), used = syns, executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCSchwannRNASeq.R")
}

dat<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "")
kinome <- read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)
sch.kin.tx <- kinome %>% filter(cellLine1 == "HS01" & cellLine2 == "HS11" & pval_adj < 0.05 & time1 == "24h" & time2 == "24h")

#HS01 - HS11 baseline
sch.kin<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(samples2, comparison=="HS01.DMSO-HS11.DMSO") %>% select(Hugo_Gene, log2FoldChange, padj) %>% inner_join(sch.kin)

ggplot(bar, aes(y = log2FoldChange, x = Mean_Kinome_Ratio)) +
  geom_point() +
  geom_point(data = bar %>% filter(padj<=0.1)) +
  geom_label_repel(data = bar %>% filter(abs(log2FoldChange)>0.5 | abs(Mean_Kinome_Ratio)>0.25), 
                   aes(label = Hugo_Gene, fill = padj<0.1)) +
  scale_fill_manual(values = c('TRUE' = '#3FA34D', 'FALSE' = '#A09F9D'), guide = "none") +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0))

#HS01 - HS11 tx
sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(samples2, comparison=="HS01.DMSO-HS11.DMSO") %>% select(Hugo_Gene, log2FoldChange, padj) %>% right_join(sch.kin)

ggplot(bar, aes(y = log2FoldChange, x = Mean_Kinome_Ratio)) +
  geom_point() +
  geom_point(data = bar %>% filter(padj<=0.1)) +
  geom_label_repel(data = bar %>% filter(abs(log2FoldChange)>0.5 | abs(Mean_Kinome_Ratio)>0.25), 
                   aes(label = Hugo_Gene, fill = padj<0.1)) +
  scale_fill_manual(values = c('TRUE' = '#3FA34D', 'FALSE' = '#A09F9D'), guide = "none") +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0))

