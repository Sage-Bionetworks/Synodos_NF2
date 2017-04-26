library(synapseClient)
library(dplyr)
library(openxlsx)
library(tidyr)
library(ggplot2)
synapseLogin()

mgh<-read.table(synGet("syn6166525")@filePath, sep = "\t", header = TRUE, comment.char = "")
mgh <- separate(mgh,geneName, into=c("gene"), sep ="\\|")
m.dmso <- filter(mgh, diffExptest == "HS01.DMSO-HS11.DMSO")
m.cudc <- filter(mgh, diffExptest == "HS01.CUDC907-HS11.CUDC907")
m.gsk <- filter(mgh, diffExptest == "HS01.GSK2126458-HS11.GSK2126458")
m.pano <- filter(mgh, diffExptest == "HS01.Panobinostat-HS11.Panobinostat")

g.hs01.dmso<-read.xlsx(synGet("syn5033489")@filePath, 1) 
g.hs01.cudc<-read.xlsx(synGet("syn5033361")@filePath, 1)
g.hs01.gsk<-read.xlsx(synGet("syn5033612")@filePath, 1)
g.hs01.pano<-read.xlsx(synGet("syn5033717")@filePath, 1)
g.hs11.dmso<-read.xlsx(synGet("syn5034722")@filePath, 1)
g.hs11.cudc<-read.xlsx(synGet("syn5033847")@filePath, 1)
g.hs11.gsk<-read.xlsx(synGet("syn5034724")@filePath, 1)
g.hs11.pano<-read.xlsx(synGet("syn5034726")@filePath, 1)

g.hs01.dmso<-g.hs01.dmso %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs01.cudc<-g.hs01.cudc %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs01.gsk<-g.hs01.gsk %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs01.pano<-g.hs01.pano %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs11.dmso<-g.hs11.dmso %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs11.cudc<-g.hs11.cudc %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs11.gsk<-g.hs11.gsk %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")
g.hs11.pano<-g.hs11.pano %>% tidyr::separate(gene_id, into=c("gene"), sep ="\\|") %>% filter(gene != "?" & gene != "SLC35E2")

g.dmso <- as.data.frame(log(g.hs01.dmso$TPM/g.hs11.dmso$TPM))
colnames(g.dmso)<- "logFC.unc"
g.dmso$gene <- g.hs01.dmso$gene

g.gsk <- as.data.frame(log(g.hs01.gsk$TPM/g.hs11.gsk$TPM))
colnames(g.gsk)<- "logFC.unc"
g.gsk$gene <- g.hs01.gsk$gene

g.cudc <- as.data.frame(log(g.hs01.cudc$TPM/g.hs11.cudc$TPM))
colnames(g.cudc)<- "logFC.unc"
g.cudc$gene <- g.hs01.cudc$gene

g.pano <- as.data.frame(log(g.hs01.pano$TPM/g.hs11.pano$TPM))
colnames(g.pano)<- "logFC.unc"
g.pano$gene <- g.hs01.pano$gene

g.dmso <- filter(g.dmso, logFC.unc != "NaN" & logFC.unc != Inf &  logFC.unc != -Inf)
g.gsk <- filter(g.gsk, logFC.unc != "NaN" & logFC.unc != Inf &  logFC.unc != -Inf)
g.cudc <- filter(g.cudc, logFC.unc != "NaN" & logFC.unc != Inf &  logFC.unc != -Inf)
g.pano <- filter(g.pano, logFC.unc != "NaN" & logFC.unc != Inf &  logFC.unc != -Inf)

dmso <- full_join(g.dmso, select(m.dmso, gene, logFC)) %>%  filter(!is.na(logFC) & !is.na(logFC.unc))
coef(lm(logFC.unc ~ logFC, data = dmso))
summary(lm(logFC.unc ~ logFC, data = dmso))

cudc <- full_join(g.cudc, select(m.cudc, gene, logFC)) %>%  filter(!is.na(logFC) & !is.na(logFC.unc))
coef(lm(logFC.unc ~ logFC, data = cudc))
summary(lm(logFC.unc ~ logFC, data = cudc))

gsk <- full_join(g.gsk, select(m.gsk, gene, logFC)) %>%  filter(!is.na(logFC) & !is.na(logFC.unc))
coef(lm(logFC.unc ~ logFC, data = gsk))
summary(lm(logFC.unc ~ logFC, data = gsk))

pano <- full_join(g.pano, select(m.pano, gene, logFC)) %>%  filter(!is.na(logFC) & !is.na(logFC.unc))
coef(lm(logFC.unc ~ logFC, data = pano))
summary(lm(logFC.unc ~ logFC, data = pano))

theme_update(plot.title = element_text(hjust = 0.5))

ggplot(dmso, aes(x=logFC, y=logFC.unc)) +
  geom_point() + 
  geom_abline(intercept = 0.1057, slope = 0.1548) +
  labs(title = "DMSO, R^2 = 0.0829", x = "logFC(MGH)", y = "logFC(UNC)")

ggplot(g.hs01.dmso, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.dmso, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5) + 
  labs(title = "HS01 DMSO", x = "log(CPM)", y = "density")

ggplot(g.hs11.pano, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.pano, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5) + 
  labs(title = "HS11 DMSO", x = "log(CPM)", y = "density")

ggplot(cudc, aes(x=logFC, y=logFC.unc)) +
  geom_point() + 
  geom_abline(intercept = 0.0479, slope = 0.0201) +
  labs(title = "CUDC, R^2 = 0.0013", x = "logFC(MGH)", y = "logFC(UNC)")
  
ggplot(g.hs01.cudc, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.cudc, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5)  + 
  labs(title = "HS01 CUDC", x = "log(CPM)", y = "density")

ggplot(g.hs11.pano, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.pano, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5) + 
  labs(title = "HS11 CUDC", x = "log(CPM)", y = "density")

ggplot(gsk, aes(x=logFC, y=logFC.unc)) +
  geom_point() + 
  geom_abline(intercept = -0.0405, slope = 0.1510) +
  labs(title = "GSK, R^2 = 0.0536", x = "logFC(MGH)", y = "logFC(UNC)")
  
ggplot(g.hs01.gsk, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.gsk, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5)  + 
  labs(title = "HS01 GSK", x = "log(CPM)", y = "density")

ggplot(g.hs11.gsk, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.pano, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5)  + 
  labs(title = "HS11 GSK", x = "log(CPM)", y = "density")

ggplot(pano, aes(x=logFC, y=logFC.unc)) +
  geom_point() + 
  geom_abline(intercept = -0.0556, slope = -0.0201) +
  labs(title = "PANO, R^2 = 0.0011", x = "logFC(MGH)", y = "logFC(UNC)")

ggplot(g.hs01.pano, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.pano, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5)  + 
  labs(title = "HS01 Pano", x = "log(CPM)", y = "density")

ggplot(g.hs11.pano, aes(x = log(TPM))) +
  geom_density(kernel = "gaussian", fill = "blue", alpha = 0.5) +
  geom_density(data = m.pano, aes(x = logCPM), kernel = "gaussian", fill = "red", alpha = 0.5)  + 
  labs(title = "HS11 Pano", x = "log(CPM)", y = "density")


