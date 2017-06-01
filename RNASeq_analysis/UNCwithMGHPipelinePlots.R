library(synapseClient)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggrepel)
synapseLogin()

files <- c("syn9884535", "syn9884502", "syn9884491", "syn9884548", "syn9884513",
           "syn9884562", "syn9884524", "syn9884581", "syn9884581", "syn9884592",
           "syn9884617")

degenes<-lapply(files, function(x){
  bar<-synGet(x)
  foo<-read.table(bar@filePath, header = T)
  comp <- bar@fileHandle$fileName
  comp <- gsub("_edgeR_quasilikelihoodFtest.txt", "", comp)
  foo$comparison <- c(rep(comp, nrow(foo)))
  foo$Hugo_Gene <- rownames(foo)
  foo
})

degenes <- ldply(degenes)

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCwithMGHPipelinePlots.R"
write.table(degenes, "schwannoma_degenes_reseq_edgeR.txt", sep = "\t")
synStore(File("schwannoma_degenes_reseq_edgeR.txt", parentId="syn9884455"), 
         used = files,
         executed = this.file)

degenes <- read.table(synGet("syn9884855")@filePath, header = TRUE, sep = "\t")

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
  ggsave(paste("UNCDEgenePlots_",x,".png", sep = ""), height = 8, width = 5)
  synStore(File(paste("UNCDEgenePlots_",x,".png", sep = ""), parentId="syn9884859"), used = "syn9884855", executed = this.file)  

  bar3 <- distinct(rbind((bar2 %>% top_n(15, logFC)), (bar2 %>% top_n(15, -logFC))))
  
  p<-ggplot(bar3, aes(x=Hugo_Gene %>% reorder(logFC), y=logFC, fill = logFC)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots_",x,"top30.png", sep = ""), height = 8, width = 5)
  synStore(File(paste("UNCDEgenePlots_",x,"top30.png", sep = ""), parentId="syn9884859"), used = "syn9884855", executed = this.file)  
  
}

dat<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "")
kinome <- read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)
sch.kin.tx <- kinome %>% filter(cellLine1 == "HS01" & cellLine2 == "HS11" & pval_adj < 0.05 & time1 == "24h" & time2 == "24h")

#HS01 - HS11 baseline
sch.kin<-dat %>% 
  filter(cellLine=="HS01", referenceSample=="HS11") %>%
  select(Gene, log2ratio) %>% 
  group_by(Gene) %>% 
  dplyr::summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(degenes, comparison=="HS01DMSOvsHS11DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin)

ggplot(bar, aes(y = logFC, x = Mean_Kinome_Ratio)) +
  geom_point() +
  geom_point(data = bar %>% filter(BH<=0.1)) +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.5 | abs(Mean_Kinome_Ratio)>0.25), 
                   aes(label = Hugo_Gene, fill = BH<0.1)) +
  scale_fill_manual(values = c('TRUE' = '#3FA34D', 'FALSE' = '#A09F9D'), guide = "none") +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0))

#HS01 - HS11 tx
sch.kin.cudc<-sch.kin.tx %>% 
  filter(drug == "CUDC") %>%
  select(protein, log2ratio) %>% 
  group_by(protein) %>% 
  dplyr::summarise(mean(log2ratio)) %>%
  arrange(desc(`mean(log2ratio)`))

colnames(sch.kin.cudc) <- c("Hugo_Gene", "Mean_Kinome_Ratio")

bar <- filter(degenes, comparison=="HS01DMSOvsHS11DMSO") %>% select(Hugo_Gene, logFC, BH) %>% inner_join(sch.kin)

ggplot(bar, aes(y = logFC, x = Mean_Kinome_Ratio)) +
  geom_point() +
  geom_point(data = bar %>% filter(BH<=0.1)) +
  geom_label_repel(data = bar %>% filter(abs(logFC)>0.5 | abs(Mean_Kinome_Ratio)>0.25), 
                   aes(label = Hugo_Gene, fill = BH<0.1)) +
  scale_fill_manual(values = c('TRUE' = '#3FA34D', 'FALSE' = '#A09F9D'), guide = "none") +
  geom_hline(aes(yintercept = 0)) + 
  geom_vline(aes(xintercept = 0))

