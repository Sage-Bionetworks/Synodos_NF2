library(synapseClient)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
synapseLogin()

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/MGHwithMGHPipelinePlots.R"

degenes <- read.table(synGet("syn6166525")@filePath, header = TRUE, sep = "\t")
degenes$Hugo_Gene <- gsub("\\|.+$", "", degenes$geneName)
write.table(degenes,"MGHDEgenes.txt", sep = "\t")

for(x in unique(degenes$diffExptest)){
  print(x)
  bar <- filter(degenes, diffExptest==x)
  bar2 <- filter(bar, BH < 0.05) %>% group_by(Hugo_Gene) %>% top_n(1, -logFC) %>% ungroup()
  
  p<-ggplot(bar2, aes(x=Hugo_Gene %>% reorder(logFC) , y=logFC, fill = logFC)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 8)) +
    labs(x = paste("Gene (n = ", nrow(bar2),")", sep =""), y = "log2(FC), BH<0.05", title = x) +
    coord_flip()
  ggsave(paste("MGHDEgenePlots_",x,".png", sep = ""), height = 8, width = 5)
  synStore(File(paste("MGHDEgenePlots_",x,".png", sep = ""), parentId="syn9884859"), used = "syn6166525", executed = this.file)  
  
  bar3 <- distinct(rbind((bar2 %>% top_n(15, logFC)), (bar2 %>% top_n(15, -logFC))))
  
  p<-ggplot(bar3, aes(x=Hugo_Gene %>% reorder(logFC), y=logFC, fill = logFC)) +
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("MGHDEgenePlots_",x,"top30.png", sep = ""), height = 8, width = 5)
  synStore(File(paste("MGHDEgenePlots_",x,"top30.png", sep = ""), parentId="syn9884859"), used = "syn6166525", executed = this.file)  
  
}



