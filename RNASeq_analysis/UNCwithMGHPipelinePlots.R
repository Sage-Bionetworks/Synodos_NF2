library(synapseClient)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(pheatmap)
library(VennDiagram)
synapseLogin()

##for merged table 

files <- c("syn9884535", "syn9884502", "syn9884491", "syn9884548", "syn9884513",
           "syn9884562", "syn9884524", "syn9884581", "syn9884592", "syn9884617")

degenes<-lapply(files, function(x){
 bar<-synGet(x)
 foo<-read.table(bar@filePath, header = T)
  comp <- bar@fileHandle$fileName
  comp <- gsub("_edgeR_quasilikelihoodFtest.txt", "", comp)
  foo$comparison <- c(rep(comp, nrow(foo)))
  foo$Hugo_Gene <- foo$gene
  foo <- select(foo, -gene)
  foo
})

degenes <- ldply(degenes)

this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/RNASeq_analysis/UNCwithMGHPipelinePlots.R"

write.table(degenes, "schwannoma_degenes_reseq_edgeR.txt", sep = "\t")
synStore(File("schwannoma_degenes_reseq_edgeR.txt", parentId="syn9884455"), 
         used = files,
         executed = this.file)

##merged table is on synapse
degenes <- read.table(synGet("syn9884855")@filePath, header = TRUE, sep = "\t")

for(x in unique(degenes$comparison)){
  print(x)
  bar <- filter(degenes, comparison==x)
  bar2 <- filter(bar, BH < 0.1) %>% group_by(Hugo_Gene) %>% top_n(1, -logFC) %>% ungroup()
  
  p<-ggplot(bar2, aes(x=Hugo_Gene %>% reorder(logFC) , y=logFC, fill = logFC)) +
    theme_bw() +
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
    theme_bw() +
    geom_bar(stat = "identity") +
    scale_fill_viridis(option = "C") +
    theme(legend.position="none", axis.text.x=element_text(angle=60,hjust=1)) +
    theme(plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 13)) +
    labs(x = paste("Gene (n = ", nrow(bar3),")", sep =""), y = paste("log2(FC)", nrow(bar3), "most significant"), title = x) +
    coord_flip()
  ggsave(paste("UNCDEgenePlots_",x,"top30.png", sep = ""), height = 8, width = 5)
  synStore(File(paste("UNCDEgenePlots_",x,"top30.png", sep = ""), parentId="syn9884859"), used = "syn9884855", executed = this.file)  
  
}



##top30 heatmap

top20 <- degenes %>%   
  filter(comparison %in% c("HS01DMSOvsHS11DMSO", "HS01CUDC907vsHS11CUDC907",
                           "HS01GSK2126458vsHS11GSK2126458", "HS01panobinostatvsHS11panobinostat")) %>% 
  group_by(comparison) %>% 
  top_n(20, abs(logFC))

top.genes <- degenes %>% 
  filter(comparison %in% c("HS01DMSOvsHS11DMSO", "HS01CUDC907vsHS11CUDC907",
                           "HS01GSK2126458vsHS11GSK2126458", "HS01panobinostatvsHS11panobinostat")) %>% 
  filter(Hugo_Gene %in% c(unique(top20$Hugo_Gene))) %>% 
  select(logFC, comparison, Hugo_Gene) %>% 
  spread(comparison, logFC)

rownames(top.genes) <- top.genes$Hugo_Gene
top.genes <- top.genes[,-1]
pdf("top-genes.pdf", height = 20)
pheatmap(top.genes, cellwidth = 10, cellheight = 10, cluster_rows = FALSE)
dev.off()
##venn diagram

sig.genes <- degenes %>% filter(BH<0.1 & logFC>1)

HS01.HS11.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS11CUDC907") %>% 
  select(Hugo_Gene)
HS01.HS11.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS11GSK2126458") %>% 
  select(Hugo_Gene)
HS01.HS11.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS11panobinostat") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS11.CUDC, HS01.HS11.GSK, HS01.HS11.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_hs11_venn_up_noDMSO.tiff", col = "red", fill = "red", label.col = "darkred", alpha = 0.25)

HS01.HS11.DMSO <- sig.genes %>% filter(comparison == "HS01DMSOvsHS11DMSO") %>% 
  select(Hugo_Gene)
HS01.HS11.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS11CUDC907") %>% 
  select(Hugo_Gene)
HS01.HS11.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS11GSK2126458") %>% 
  select(Hugo_Gene)
HS01.HS11.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS11panobinostat") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS11.DMSO, HS01.HS11.CUDC, HS01.HS11.GSK, HS01.HS11.PANO))
names(list) <- c("DMSO", "CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_hs11_venn_up.tiff", col = "red", fill = "red", label.col = "darkred", alpha = 0.25)

HS01.HS01.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS01DMSO") %>% 
  select(Hugo_Gene)
HS01.HS01.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS01DMSO") %>% 
  select(Hugo_Gene)
HS01.HS01.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS01DMSO") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS01.CUDC, HS01.HS01.GSK, HS01.HS01.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_bydrug_venn_up.tiff", col = "red", fill = "red", label.col = "darkred", alpha = 0.25)

HS11.HS11.CUDC <- sig.genes %>% filter(comparison == "HS11CUDC907vsHS11DMSO") %>% 
  select(Hugo_Gene)
HS11.HS11.GSK <- sig.genes %>% filter(comparison == "HS11GSK2126458vsHS11DMSO") %>% 
  select(Hugo_Gene)
HS11.HS11.PANO <- sig.genes %>% filter(comparison == "HS11panobinostatvsHS11DMSO") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS11.HS11.CUDC, HS11.HS11.GSK, HS11.HS11.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs11_bydrug_venn_up.tiff", col = "red", fill = "red", label.col = "darkred", alpha = 0.25)


sig.genes <- degenes %>% filter(BH<0.1 & logFC<=-1)

HS01.HS11.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS11CUDC907") %>% 
  select(Hugo_Gene)
HS01.HS11.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS11GSK2126458") %>% 
  select(Hugo_Gene)
HS01.HS11.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS11panobinostat") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS11.CUDC, HS01.HS11.GSK, HS01.HS11.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_hs11_venn_down_noDMSO.tiff", col = "blue", fill = "blue", label.col = "darkblue", alpha = 0.25)


HS01.HS11.DMSO <- sig.genes %>% filter(comparison == "HS01DMSOvsHS11DMSO") %>% 
  select(Hugo_Gene)
HS01.HS11.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS11CUDC907") %>% 
  select(Hugo_Gene)
HS01.HS11.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS11GSK2126458") %>% 
  select(Hugo_Gene)
HS01.HS11.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS11panobinostat") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS11.DMSO, HS01.HS11.CUDC, HS01.HS11.GSK, HS01.HS11.PANO))
names(list) <- c("DMSO", "CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_hs11_venn_down.tiff", col = "blue", fill = "blue", label.col = "darkblue", alpha = 0.25)


HS01.HS01.CUDC <- sig.genes %>% filter(comparison == "HS01CUDC907vsHS01DMSO") %>% 
  select(Hugo_Gene)
HS01.HS01.GSK <- sig.genes %>% filter(comparison == "HS01GSK2126458vsHS01DMSO") %>% 
  select(Hugo_Gene)
HS01.HS01.PANO <- sig.genes %>% filter(comparison == "HS01panobinostatvsHS01DMSO") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS01.HS01.CUDC, HS01.HS01.GSK, HS01.HS01.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs01_bydrug_venn_down.tiff", col = "blue", fill = "blue", label.col = "darkblue", alpha = 0.25)


HS11.HS11.CUDC <- sig.genes %>% filter(comparison == "HS11CUDC907vsHS11DMSO") %>% 
  select(Hugo_Gene)
HS11.HS11.GSK <- sig.genes %>% filter(comparison == "HS11GSK2126458vsHS11DMSO") %>% 
  select(Hugo_Gene)
HS11.HS11.PANO <- sig.genes %>% filter(comparison == "HS11panobinostatvsHS11DMSO") %>% 
  select(Hugo_Gene)

list <- as.list(c(HS11.HS11.CUDC, HS11.HS11.GSK, HS11.HS11.PANO))
names(list) <- c("CUDC907", "GSK2126458", "Panobinostat")
venn.diagram(list, "hs11_bydrug_venn_down.tiff", col = "blue", fill = "blue", label.col = "darkblue", alpha = 0.25)

for(x in unique(degenes$comparison)){
  p <- filter(degenes, comparison == x & BH<0.1 & abs(logFC)>1)
  write.table(p, paste(x,"_degenes.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}



