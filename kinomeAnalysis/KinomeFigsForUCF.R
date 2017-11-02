library(tidyr)
library(dplyr)
library(ggplot2)
library(synapseClient)
synapseLogin()

##Get intersecting kinases (DE after treatment)

dekin<-read.table(synGet("syn4975368")@filePath, sep = "\t", header = TRUE)

HS01.HS11.CUDC.up <- dekin %>% filter(drug == "CUDC", cellLine1 =="HS01", cellLine2 =="HS11") %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.05) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1-avgRatio_proteinRep_cond2) %>% 
  top_n(10, avgRatio) %>% 
  filter(avgRatio > 0)

HS01.HS11.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="HS01", cellLine2 =="HS11") %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.05) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1-avgRatio_proteinRep_cond2) %>% 
  top_n(10,-avgRatio) %>% 
  filter(avgRatio < 0) %>% 
  bind_rows(HS01.HS11.CUDC.up)

HS01.HS11.CUDC$protein <- reorder(HS01.HS11.CUDC$protein, HS01.HS11.CUDC$avgRatio)

ggplot(data = HS01.HS11.CUDC, aes(x=protein, y=avgRatio)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = avgRatio)) +
  scale_fill_viridis_c(option = "C") +
  coord_flip() +
  labs(y = "log2ratio (mean)", x = "Gene", title = "HS01 vs HS11, CUDC") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("HS01_HS11_kinome_24h_cudc.png", height = 5, width = 3)

# HS01.HS11.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="HS01", cellLine2 =="HS11") %>% 
#   filter(time1 == "2h", time2 == "2h", pval_adj <= 0.05) %>% 
#   mutate(avgRatio = avgRatio_proteinRep_cond1-avgRatio_proteinRep_cond2)
# 
# HS01.HS11.CUDC$protein <- reorder(HS01.HS11.CUDC$protein, HS01.HS11.CUDC$avgRatio)
# 
# ggplot(data = HS01.HS11.CUDC, aes(x=protein, y=avgRatio)) +
#   geom_bar(stat = "identity", position = "dodge", aes(fill = avgRatio)) +
#   scale_fill_viridis_c(option = "C") +
#   coord_flip() +
#   labs(y = "log2ratio (mean)", x = "Gene", title = "HS01 vs HS11, CUDC") + 
#   theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
# 
# ggsave("HS01_HS11_kinome_2h_cudc.png", height = 5, width = 3)


MS03.MS12.CUDC.up <- dekin %>% filter(drug == "CUDC", cellLine1 =="MS03", cellLine2 =="MS12") %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.05) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1-avgRatio_proteinRep_cond2) %>% 
  top_n(10, avgRatio)

MS03.MS12.CUDC <- dekin %>% filter(drug == "CUDC", cellLine1 =="MS03", cellLine2 =="MS12") %>% 
  filter(time1 == "24h", time2 == "24h", pval_adj <= 0.05) %>% 
  mutate(avgRatio = avgRatio_proteinRep_cond1-avgRatio_proteinRep_cond2) %>% 
  top_n(10, -avgRatio) %>% 
  bind_rows(MS03.MS12.CUDC.up)


MS03.MS12.CUDC$protein <- reorder(MS03.MS12.CUDC$protein, MS03.MS12.CUDC$avgRatio)

ggplot(data = MS03.MS12.CUDC, aes(x=protein, y=avgRatio)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = avgRatio)) +
  scale_fill_viridis_c(option = "C") +
  coord_flip() +
  labs(y = "log2ratio (mean)", x = "Gene", title = "MS03 vs MS12, CUDC") + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave("MS03_MS12_kinome_24h_cudc.png", height = 5, width = 3)
