library(synapseClient)
library(dplyr)
synapseLogin()

degenes.sage<-read.table(synGet("syn6038243")@filePath, sep = "\t", header = TRUE, comment.char = "")
degenes.mgh<-read.table(synGet("syn6166525")@filePath, sep = "\t", header = TRUE, comment.char = "")
degenes.ms<-read.table(synGet("syn7437782")@filePath, sep = "\t", header = TRUE, comment.char = "")

###Sage Baseline
hs01hs11.s.base <- nrow(filter(degenes.sage, diffExptest == 'HS01.DMSO-HS11.DMSO' & adj.P.Val <= 0.1))
syn5syn1.s.base <- nrow(filter(degenes.sage, diffExptest == 'Syn5.DMSO-Syn1.DMSO' & adj.P.Val <= 0.1))
syn6syn1.s.base <- nrow(filter(degenes.sage, diffExptest == 'Syn6.DMSO-Syn1.DMSO' & adj.P.Val <= 0.1))

###MGH Baseline
hs01hs11.m.base <- nrow(filter(degenes.mgh, diffExptest == 'HS01.DMSO-HS11.DMSO' & bonferroni <= 0.1))
syn5syn1.m.base <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.DMSO-Syn1.DMSO' & bonferroni <= 0.1))
syn6syn1.m.base <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.DMSO-Syn1.DMSO' & bonferroni <= 0.1))

###Sage CUDC
hs01hs11.s.cudc <- nrow(filter(degenes.sage, diffExptest == 'HS01.CUDC907-HS11.CUDC907' & adj.P.Val <= 0.1))
syn5syn1.s.cudc <- nrow(filter(degenes.sage, diffExptest == 'Syn5.CUDC907-Syn1.CUDC907' & adj.P.Val <= 0.1))
syn6syn1.s.cudc <- nrow(filter(degenes.sage, diffExptest == 'Syn6.CUDC907-Syn1.CUDC907' & adj.P.Val <= 0.1))

###MGH CUDC
hs01hs11.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'HS01.CUDC907-HS11.CUDC907' & bonferroni <= 0.1))
syn5syn1.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.CUDC907-Syn1.CUDC907' & bonferroni <= 0.1))
syn6syn1.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.CUDC907-Syn1.CUDC907' & bonferroni <= 0.1))

###Sage GSK
hs01hs11.s.gsk <- nrow(filter(degenes.sage, diffExptest == 'HS01.GSK2126458-HS11.GSK2126458' & adj.P.Val <= 0.1))
syn5syn1.s.gsk <- nrow(filter(degenes.sage, diffExptest == 'Syn5.GSK2126458-Syn1.GSK2126458' & adj.P.Val <= 0.1))
syn6syn1.s.gsk <- nrow(filter(degenes.sage, diffExptest == 'Syn6.GSK2126458-Syn1.GSK2126458' & adj.P.Val <= 0.1))

###MGH GSK
hs01hs11.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'HS01.GSK2126458-HS11.GSK2126458' & bonferroni <= 0.1))
syn5syn1.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.GSK2126458-Syn1.GSK2126458' & bonferroni <= 0.1))
syn6syn1.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.GSK2126458-Syn1.GSK2126458' & bonferroni <= 0.1))

###Sage Pano
hs01hs11.s.pano <- nrow(filter(degenes.sage, diffExptest == 'HS01.Panobinostat-HS11.Panobinostat' & adj.P.Val <= 0.1))
syn5syn1.s.pano <- nrow(filter(degenes.sage, diffExptest == 'Syn5.Panobinostat-Syn1.Panobinostat' & adj.P.Val <= 0.1))
syn6syn1.s.pano <- nrow(filter(degenes.sage, diffExptest == 'Syn6.Panobinostat-Syn1.Panobinostat' & adj.P.Val <= 0.1))

###MGH Pano
hs01hs11.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'HS01.Panobinostat-HS11.Panobinostat' & bonferroni <= 0.1))
syn5syn1.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.Panobinostat-Syn1.Panobinostat' & bonferroni <= 0.1))
syn6syn1.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.Panobinostat-Syn1.Panobinostat' & bonferroni <= 0.1))

###Mouse
ms03ms12.base <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "DMSO" & treatment2 == "DMSO" & bonferroni <= 0.1))
ms03ms12.cudc <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "CUDC" & treatment2 == "CUDC" & bonferroni <= 0.1))
ms03ms12.gsk <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "GSK458" & treatment2 == "GSK458" & bonferroni <= 0.1))
ms03ms12.pano <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "Pano" & treatment2 == "Pano" & bonferroni <= 0.1))


####tables
sage<-matrix(data=c(hs01hs11.s.base,hs01hs11.s.cudc,hs01hs11.s.gsk,hs01hs11.s.pano,
               syn5syn1.s.base,syn5syn1.s.cudc,syn5syn1.s.gsk, syn5syn1.s.pano,
               syn6syn1.s.base,syn6syn1.s.cudc,syn6syn1.s.gsk, syn6syn1.s.pano), nrow = 3, ncol = 4, byrow=TRUE)
rownames(sage) <- c("HS01-HS11", "Syn5-Syn1", "Syn6-Syn1")
colnames(sage) <- c("baseline", "CUDC907", "GSK2126458", "panobinostat")

mgh<-matrix(c(hs01hs11.m.base,hs01hs11.m.cudc,hs01hs11.m.gsk,hs01hs11.m.pano,
              ms03ms12.base, ms03ms12.cudc, ms03ms12.gsk, ms03ms12.pano,
              syn5syn1.m.base,syn5syn1.m.cudc,syn5syn1.m.gsk, syn5syn1.m.pano,
              syn6syn1.m.base,syn6syn1.m.cudc,syn6syn1.m.gsk, syn6syn1.m.pano), nrow = 4, ncol = 4, byrow=TRUE) 
rownames(mgh) <- c("HS01-HS11", "MS01-MS11", "Syn5-Syn1", "Syn6-Syn1")
colnames(mgh) <- c("baseline", "CUDC907", "GSK2126458", "panobinostat")


## BH instead of bon
###MGH Baseline
hs01hs11.m.base <- nrow(filter(degenes.mgh, diffExptest == 'HS01.DMSO-HS11.DMSO' & BH <= 0.1))
syn5syn1.m.base <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.DMSO-Syn1.DMSO' & BH <= 0.1))
syn6syn1.m.base <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.DMSO-Syn1.DMSO' & BH <= 0.1))

###MGH CUDC
hs01hs11.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'HS01.CUDC907-HS11.CUDC907' & BH <= 0.1))
syn5syn1.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.CUDC907-Syn1.CUDC907' & BH <= 0.1))
syn6syn1.m.cudc <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.CUDC907-Syn1.CUDC907' & BH <= 0.1))

###MGH GSK
hs01hs11.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'HS01.GSK2126458-HS11.GSK2126458' & BH <= 0.1))
syn5syn1.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.GSK2126458-Syn1.GSK2126458' & BH <= 0.1))
syn6syn1.m.gsk <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.GSK2126458-Syn1.GSK2126458' & BH <= 0.1))

###MGH Pano
hs01hs11.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'HS01.Panobinostat-HS11.Panobinostat' & BH <= 0.1))
syn5syn1.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'Syn5.Panobinostat-Syn1.Panobinostat' & BH <= 0.1))
syn6syn1.m.pano <- nrow(filter(degenes.mgh, diffExptest == 'Syn6.Panobinostat-Syn1.Panobinostat' & BH <= 0.1))

###Mouse
ms03ms12.base <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "DMSO" & treatment2 == "DMSO" & BH <= 0.1))
ms03ms12.cudc <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "CUDC" & treatment2 == "CUDC" & BH <= 0.1))
ms03ms12.gsk <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "GSK458" & treatment2 == "GSK458" & BH <= 0.1))
ms03ms12.pano <- nrow(filter(degenes.ms, cellLine1 == "MS03" & cellLine2 == "MS12" & treatment1 == "Pano" & treatment2 == "Pano" & BH <= 0.1))

mgh.bh<-matrix(c(hs01hs11.m.base,hs01hs11.m.cudc,hs01hs11.m.gsk,hs01hs11.m.pano,
              ms03ms12.base, ms03ms12.cudc, ms03ms12.gsk, ms03ms12.pano,
              syn5syn1.m.base,syn5syn1.m.cudc,syn5syn1.m.gsk, syn5syn1.m.pano,
              syn6syn1.m.base,syn6syn1.m.cudc,syn6syn1.m.gsk, syn6syn1.m.pano), nrow = 4, ncol = 4, byrow=TRUE) 
rownames(mgh.bh) <- c("HS01-HS11", "MS01-MS11", "Syn5-Syn1", "Syn6-Syn1")
colnames(mgh.bh) <- c("baseline", "CUDC907", "GSK2126458", "panobinostat")


