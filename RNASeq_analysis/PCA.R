library(synapseClient)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)
synapseLogin()

syn<-read.table(synGet("syn6045994")@filePath,sep="\t",header=T)

syn<-separate(syn,gene,sep="\\|",c("gene","ensembl"))%>%select(gene,Syn1.2_CUDC.907,Syn1.3_CUDC.907,Syn5.1_CUDC.907,Syn5.2_CUDC.907,Syn5.3_CUDC.907,Syn1.1_Panobinostat,
Syn1.2_Panobinostat,Syn1.3_Panobinostat,Syn5.1_Panobinostat,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,Syn5.1_DMSO,
Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn1.1_CUDC.907,Syn5.1_GSK2126458,Syn5.2_GSK2126458,Syn5.3_GSK2126458,
Syn6.1_GSK2126458,Syn6.2_GSK2126458,Syn6.3_DMSO,Syn6.3_CUDC.907,Syn6.3_Panobinostat,Syn6.3_GSK2126458,Syn6.1_CUDC.907,
Syn6.2_CUDC.907,Syn5.2_Panobinostat,Syn5.3_Panobinostat,Syn6.1_Panobinostat,Syn6.2_Panobinostat,Syn1.1_GSK2126458,
Syn1.2_GSK2126458,Syn1.3_GSK2126458)

hs<-read.table(synGet("syn9779230")@filePath,sep="\t",header=T)
colnames(hs)[1]<-"gene"


ms<-read.table(synGet("syn7467412")@filePath,sep="\t",header=T)

  
all<-inner_join(hs,syn)
all<-aggregate(all[,-1],list(all$gene),mean)
rownames(all)<-all$Group.1
all<-all[,-1]

all.rank<-apply(t(all),1,rank,ties.method='min')

baseline <- as.data.frame(all.rank) %>% select(HS11_DMSO_Run1_S5,HS11_DMSO_Run2_S5,HS11_DMSO_Run3_S5,HS11_DMSO_Run4_S9,HS01_DMSO_Run1_S1,
                            HS01_DMSO_Run2_S1,HS01_DMSO_Run3_S1,HS01_DMSO_Run4_S9,Syn1.1_DMSO,Syn1.2_DMSO,Syn1.3_DMSO,
                            Syn5.1_DMSO,Syn5.2_DMSO,Syn5.3_DMSO,Syn6.1_DMSO,Syn6.2_DMSO,Syn6.3_DMSO)
pca.baseline <- prcomp(t(baseline), center = T)

pcaplot <- as.data.frame(pca.baseline$x)
pcaplot$samp <- rownames(pcaplot)

ggplot(pcaplot, aes(x = PC1, y = PC2)) +
  geom_point() +
  geom_text_repel(aes(label = samp))

