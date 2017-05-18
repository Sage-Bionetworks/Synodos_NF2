library(synapseClient)
library(dplyr)
library(data.table)
library(ggplot2)
library(viridis)
library(ggrepel)
synapseLogin()

samples<-read.table(synGet("syn6166525")@filePath, sep = "\t", header = TRUE)
samples$Hugo_Gene <- gsub("\\|.+","",samples$geneName)

for(x in unique(samples$diffExptest)){
  print(x)
  bar <- filter(samples, diffExptest==x)
  bar2 <- filter(bar, BH < 0.1) %>% group_by(Hugo_Gene) %>% top_n(1, -logFC) %>% ungroup()
  
  p<-ggplot(bar2, aes(x=Hugo_Gene %>% reorder(logFC) , y=logFC, fill = logFC)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 8)) +
    labs(x = paste("Gene (n =",nrow(bar2),")"), y = "log2(FC), adj p<0.1", title = x) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_flip()
  ggsave(paste("MGHDEgenePlots/",x,".png", sep = ""), height = 8, width = 6)
  
  bar3 <- distinct(rbind((bar2 %>% top_n(15, logFC)), (bar2 %>% top_n(15, -logFC))))
  
  p<-ggplot(bar3 , aes(x=Hugo_Gene %>% reorder(logFC), y=logFC, fill = logFC)) +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("MGHDEgenePlots/",x,"top30.png", sep = ""), height = 6, width = 6)
  
  synStore(File(paste("MGHDEgenePlots/",x,".png", sep = ""), parentId = "syn9779337"), used = "syn6166525", executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/MGHMeninSchwannRNASeq.R")
  synStore(File(paste("MGHDEgenePlots/",x,"top30.png", sep = ""), parentId = "syn9779337"), used = "syn6166525", executed = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/MGHMeninSchwannRNASeq.R")
}
