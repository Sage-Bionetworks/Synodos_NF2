library(synapseClient)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(ggrepel)
synapseLogin()

syns <- c("syn9773229", "syn9773231", "syn9773234", "syn9773236", 
          "syn9773238", "syn9773240", "syn9773348", "syn9773392",
          "syn9773393", "syn9773394")


samples <- lapply(syns, function(x){
  foo <- synGet(x)
  y <- read.table(foo@filePath, sep = "\t", header = TRUE)
  y$sample <- rep(foo@fileHandle$fileName, nrow(y))
  y
})

samples2<-rbindlist(samples, use.names = TRUE)

newname<-c(
  "CUDC907-HS01out.txt",
  "CUDC907-HS11out.txt",     
  "DMSO-HS11-out.txt",        
  "GSK-HS01out.txt",         
  "GSK-HS11out.txt",          
  "panobinostat-HS01out.txt",
  "panobinostat-HS11out.txt", 
  "panobinostatout.txt",     
  "GSKout.txt",               
  "CUDC907out.txt")

newname<-as.data.frame(cbind(newname,c(
  "HS01.CUDC907-HS01.DMSO",
  "HS11.CUDC907-HS11.DMSO",
  "HS01.DMSO-HS11.DMSO",
  "HS01.GSK2126458-HS01.DMSO",
  "HS11.GSK2126458-HS11.DMSO",
  "HS01.Panobinostat-HS01.DMSO",
  "HS11.Panobinostat-HS11.DMSO",
  "HS01.Panobinostat-HS11.Panobinostat",
  "HS01.GSK2126458-HS11.GSK2126458",
  "HS01.CUDC907-HS11.CUDC907"
)))
colnames(newname) <- c("sample", "comparison")

samples2 <- full_join(samples2, newname)
colnames(samples2)[colnames(samples2)=="X"] <- "Hugo_Gene"

for(x in unique(samples2$comparison)){
  print(x)
  bar <- filter(samples2, comparison==x)
  bar2 <- filter(bar, padj < 0.1) %>% group_by(Hugo_Gene) %>% top_n(1, -log2FoldChange) %>% ungroup()
  
  p<-ggplot(bar2, aes(x=Hugo_Gene %>% reorder(log2FoldChange) , y=log2FoldChange, fill = log2FoldChange)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 8)) +
    labs(x = paste("Gene (n = ", nrow(bar2),")", sep =""), y = "log2(FC), BH<0.1", title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots/",x,".png", sep = ""), height = 8, width = 5)
  
  bar3 <- distinct(rbind((bar2 %>% top_n(15, log2FoldChange)), (bar2 %>% top_n(15, -log2FoldChange))))
  
  p<-ggplot(bar3, aes(x=Hugo_Gene %>% reorder(log2FoldChange), y=log2FoldChange, fill = log2FoldChange)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots/",x,"top30.png", sep = ""), height = 8, width = 5)
  
  synStore(File(paste("UNCDEgenePlots/",x,".png", sep = ""), parentId = "syn9779338"), used = syns, executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCSchwannRNASeq.R")
  synStore(File(paste("UNCDEgenePlots/",x,"top30.png", sep = ""), parentId = "syn9779338"), used = syns, executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCSchwannRNASeq.R")
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
