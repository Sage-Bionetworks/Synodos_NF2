library(refGenome)
library(dplyr)
library(tidyr)
library(synapseClient)
synapseLogin()

men <- read.table(synGet("syn10881885")@filePath, sep = "\t", header = T) %>% 
  separate(geneName, c("hugo", "ens"), sep = "\\|") %>% 
  dplyr::select(hugo, ens, diffExptest, logFC, BH)

sch <- read.table(synGet("syn9884855")@filePath, sep = "\t", header = T) %>% 
  dplyr::rename("hugo" = Hugo_Gene, "ens" = ensembl) %>% 
  dplyr::select(hugo, ens, comparison, logFC, BH)


Syn5DMSOSyn1DMSO <- men %>% filter(diffExptest=="Syn5.DMSO-Syn1.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("1.FC" = logFC, "1.BH" = BH)
HS01DMSOHS11DMSO <- sch %>% filter(comparison=="HS01DMSOvsHS11DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("2.FC" = logFC, "2.BH" = BH)
Syn1CUDCvsSyn1DMSO <- men %>% filter(diffExptest=="Syn1.CUDC907-Syn1.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("3.FC" = logFC, "3.BH" = BH)
Syn1PanovsSyn1DMSO <- men %>% filter(diffExptest=="Syn1.Panobinostat-Syn1.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("4.FC" = logFC, "4.BH" = BH)
Syn1GSKvsSyn1CUDC <- men %>% filter(diffExptest=="Syn1.GSK2126458-Syn1.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("5.FC" = logFC, "5.BH" = BH)
Syn5CUDCvsSyn5DMSO <- men %>% filter(diffExptest=="Syn5.CUDC907-Syn5.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("6.FC" = logFC, "6.BH" = BH)
Syn5PanovsSyn5DMSO <- men %>% filter(diffExptest=="Syn5.Panobinostat-Syn5.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("7.FC" = logFC, "7.BH" = BH)
Syn5GSKvsSyn5DMSO <- men %>% filter(diffExptest=="Syn5.GSK2126458-Syn5.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("8.FC" = logFC, "8.BH" = BH)
Syn6CUDCvsSyn6DMSO <- men %>% filter(diffExptest=="Syn6.CUDC907-Syn6.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("9.FC" = logFC, "9.BH" = BH)
Syn6PanovsSyn6DMSO <- men %>% filter(diffExptest=="Syn6.Panobinostat-Syn6.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("10.FC" = logFC, "10.BH" = BH)
Syn6GSKvsSyn6DMSO <- men %>% filter(diffExptest=="Syn6.GSK2126458-Syn6.DMSO") %>% 
  select(-diffExptest) %>% 
  dplyr::rename("11.FC" = logFC, "11.BH" = BH)
HS01CUDCvsHS01DMSO <- sch %>% filter(comparison=="HS01CUDC907vsHS01DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("12.FC" = logFC, "12.BH" = BH)
HS01PanovsHS01DMSO <- sch %>% filter(comparison=="HS01panobinostatvsHS01DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("13.FC" = logFC, "13.BH" = BH)
HS01GSKvsHS01DMSO<- sch %>% filter(comparison=="HS01GSK2126458vsHS01DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("14.FC" = logFC, "14.BH" = BH)
HS11CUDCvsHS11DMSO <- sch %>% filter(comparison=="HS11CUDC907vsHS11DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("15.FC" = logFC, "15.BH" = BH)
HS11PanovsHS11DMSO <- sch %>% filter(comparison=="HS11panobinostatvsHS11DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("16.FC" = logFC, "16.BH" = BH)
HS11GSKvsHS11DMSO <- sch %>% filter(comparison=="HS11GSK2126458vsHS11DMSO") %>% 
  select(-comparison) %>% 
  dplyr::rename("17.FC" = logFC, "17.BH" = BH)

all <- full_join(Syn5DMSOSyn1DMSO, HS01DMSOHS11DMSO) %>% 
  full_join(Syn1CUDCvsSyn1DMSO) %>% 
  full_join(Syn1PanovsSyn1DMSO) %>% 
  full_join(Syn1GSKvsSyn1CUDC) %>% 
  full_join(Syn5CUDCvsSyn5DMSO) %>% 
  full_join(Syn5PanovsSyn5DMSO) %>% 
  full_join(Syn5GSKvsSyn5DMSO) %>% 
  full_join(Syn6CUDCvsSyn6DMSO) %>% 
  full_join(Syn6PanovsSyn6DMSO) %>% 
  full_join(Syn6GSKvsSyn6DMSO) %>% 
  full_join(HS01CUDCvsHS01DMSO) %>% 
  full_join(HS01PanovsHS01DMSO) %>% 
  full_join(HS01GSKvsHS01DMSO) %>% 
  full_join(HS11CUDCvsHS11DMSO) %>% 
  full_join(HS11PanovsHS11DMSO) %>% 
  full_join(HS11GSKvsHS11DMSO)

write.table(all, "all.txt", sep ="\t")
