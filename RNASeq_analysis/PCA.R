library(synapseClient)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)
synapseLogin()

##meningioma raw counts mgh
syn<-read.table(synGet("syn6045994")@filePath,sep="\t",header=T)

syn<-separate(syn,gene,sep="\\|",c("gene","ensembl"))%>%dplyr::select(gene,Syn1.2_CUDC.907,Syn1.3_CUDC.907,Syn5.1_CUDC.907,Syn5.2_CUDC.907,Syn5.3_CUDC.907,Syn1.1_Panobinostat,
Syn1.2_Panobinostat,Syn1.3_Panobinostat,Syn5.1_Panobinostat,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,Syn5.1_DMSO,
Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn1.1_CUDC.907,Syn5.1_GSK2126458,Syn5.2_GSK2126458,Syn5.3_GSK2126458,
Syn6.1_GSK2126458,Syn6.2_GSK2126458,Syn6.3_DMSO,Syn6.3_CUDC.907,Syn6.3_Panobinostat,Syn6.3_GSK2126458,Syn6.1_CUDC.907,
Syn6.2_CUDC.907,Syn5.2_Panobinostat,Syn5.3_Panobinostat,Syn6.1_Panobinostat,Syn6.2_Panobinostat,Syn1.1_GSK2126458,
Syn1.2_GSK2126458,Syn1.3_GSK2126458, Syn10.1_DMSO, Syn10.2_DMSO)

##schwannoma raw counts (unc)
hs<-read.table(synGet("syn9925491")@filePath,sep="\t",header=T)
colnames(hs)[1]<-"gene"

  
all<-inner_join(hs,syn)
all2<-all
all<-all[,-1]

all.rank<-apply(t(all),1,rank,ties.method='min')

baseline <- as.data.frame(all.rank) %>% dplyr::select(HS11_DMSO_Run1_S5,HS11_DMSO_Run2_S5,HS11_DMSO_Run3_S5,HS11_DMSO_Run4_S9,HS01_DMSO_Run1_S1,
                            HS01_DMSO_Run2_S1,HS01_DMSO_Run3_S1,HS01_DMSO_Run4_S9,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,
                            Syn5.1_DMSO,Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn6.3_DMSO, Syn10.1_DMSO, Syn10.2_DMSO)
pca.baseline <- prcomp(t(baseline), center = T)

pcaplot <- as.data.frame(pca.baseline$x)
pcaplot$samp <- rownames(pcaplot)

ggplot(pcaplot, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = samp))


##mouse raw counts mgh
ms<-read.table(synGet("syn7467412")@filePath,sep=",",header=T, comment.char = "")
ms<-dplyr::select(ms, gene, MS03_2, MS03_3, MS03_4, MS12_2,MS12_3, MS12_4)
ms$gene <- gsub("^.+\\|","",ms$gene)
ms$gene <- toupper(ms$gene)
all2<-inner_join(ms, all2)
all2 <- aggregate(. ~ gene, data = all2, FUN = sum)
rownames(all2) <- all2$gene
all2 <- all2[,-1]
all2.rank<-apply(t(all2),1,rank,ties.method='min')

baseline <- as.data.frame(all2.rank) %>% 
  dplyr::select(HS11_DMSO_Run1_S5,HS11_DMSO_Run2_S5,HS11_DMSO_Run3_S5,HS11_DMSO_Run4_S9,HS01_DMSO_Run1_S1,
                HS01_DMSO_Run2_S1,HS01_DMSO_Run3_S1,HS01_DMSO_Run4_S9,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,
                Syn5.1_DMSO,Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn6.3_DMSO, MS03_2, MS03_3, 
                MS03_4, MS12_2,MS12_3, MS12_4) %>% 
  dplyr::transmute(HS11=(HS11_DMSO_Run1_S5+HS11_DMSO_Run2_S5+HS11_DMSO_Run3_S5+HS11_DMSO_Run4_S9)/4,
                HS01=(HS01_DMSO_Run1_S1+HS01_DMSO_Run2_S1+HS01_DMSO_Run3_S1+HS01_DMSO_Run4_S9)/4,
                Syn1=(Syn1.1_DMSO+Syn1.2_DMSO+Syn1.3_DMSO)/3,
                Syn5=(Syn5.1_DMSO+Syn5.2_DMSO+Syn5.3_DMSO)/3,
                Syn6=(Syn6.1_DMSO+Syn6.2_DMSO+Syn6.3_DMSO)/3,
                MS03=(MS03_2+MS03_3+MS03_4)/3,
                MS12=(MS12_2+MS12_3+MS12_4)/3)



pca.baseline <- prcomp(t(baseline), center = T)

pcaplot <- as.data.frame(pca.baseline$x)
pcaplot$samp <- rownames(pcaplot)
vals <- summary(pca.baseline)

ggplot(pcaplot, aes(x = PC1, y = PC2)) +
  theme_bw() + 
  geom_point() +
  geom_label_repel(aes(label = samp, fill = samp), 
                   max.iter = 200000,
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.25, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("HS11" = "#FF9C99", "HS01" = "#FF9C99", "MS12" = "#91CDFF", "MS03" = "#91CDFF", "Syn6" = "#B9A6DD", "Syn5"  = "#B9A6DD","Syn1"  = "#B9A6DD"), guide = "none") +
  scale_x_continuous(name=paste0("PC1 (",signif(vals$importance[2,1]*100,3),"% variance explained)"),
                     labels = scales::comma) +
  scale_y_continuous(name=paste0("PC2 (",signif(vals$importance[2,2]*100,3),"% variance explained)"),
                     labels = scales::comma)


ggsave("RNASeqPCAPlot.png", height = 7, width = 7)



baseline <- as.data.frame(all2.rank) %>% 
  dplyr::select(HS11_DMSO_Run1_S5,HS11_DMSO_Run2_S5,HS11_DMSO_Run3_S5,HS11_DMSO_Run4_S9,HS01_DMSO_Run1_S1,
                HS01_DMSO_Run2_S1,HS01_DMSO_Run3_S1,HS01_DMSO_Run4_S9,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,
                Syn5.1_DMSO,Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn6.3_DMSO, MS03_2, MS03_3, 
                MS03_4, MS12_2,MS12_3, MS12_4, Syn10.1_DMSO, Syn10.2_DMSO) %>% 
  dplyr::transmute(HS11=(HS11_DMSO_Run1_S5+HS11_DMSO_Run2_S5+HS11_DMSO_Run3_S5+HS11_DMSO_Run4_S9)/4,
                   HS01=(HS01_DMSO_Run1_S1+HS01_DMSO_Run2_S1+HS01_DMSO_Run3_S1+HS01_DMSO_Run4_S9)/4,
                   Syn1=(Syn1.1_DMSO+Syn1.2_DMSO+Syn1.3_DMSO)/3,
                   Syn5=(Syn5.1_DMSO+Syn5.2_DMSO+Syn5.3_DMSO)/3,
                   Syn6=(Syn6.1_DMSO+Syn6.2_DMSO+Syn6.3_DMSO)/3,
                   MS03=(MS03_2+MS03_3+MS03_4)/3,
                   MS12=(MS12_2+MS12_3+MS12_4)/3,
                   Syn10=(Syn10.1_DMSO+Syn10.2_DMSO)/2)



pca.baseline <- prcomp(t(baseline), center = T)

pcaplot <- as.data.frame(pca.baseline$x)
pcaplot$samp <- rownames(pcaplot)
vals <- summary(pca.baseline)
ggplot(pcaplot, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_label_repel(aes(label = samp, fill = samp), 
                   max.iter = 20000,
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.3, "lines"),
                   size = 6) +
  scale_fill_manual(values = c("HS11" = "#FF9C99", "HS01" = "#FF9C99", "MS12" = "#FFC251", "MS03" = "#FFC251", "Syn6" = "#91CDFF", "Syn5"  = "#91CDFF","Syn1"  = "#91CDFF", "Syn10"  = "#91CDFF"), guide = "none") +
  labs(x = paste("PC1 (",signif(vals$importance[2,1]*100,3),"% variance explained)", sep = ""),
       y = paste("PC2 (",signif(vals$importance[2,2]*100,3),"% variance explained)", sep = ""))
ggsave("RNASeqPCAPlotwithSyn10.png", height = 6, width = 6)

