library(tidyverse)
library(synapseClient)
synapseLogin()

####schwannoma 

ncats <-read.table(synGet("syn8314523")@filePath, header = T, sep = "\t")

ncats.null <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS01") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgNull = mean(AUC))

ncats.wt <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS12") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgWT = mean(AUC))

ncats.min <- full_join(ncats.null, ncats.wt) %>% 
  dplyr::mutate(diffr = avgNull - avgWT) %>% 
  magrittr::set_colnames(c("ncgc", "avgNull", "avgWT", "diffr")) %>% 
  dplyr::left_join(map2) %>%

targs <- read.table(synGet("syn11609839")@filePath, header = T, sep = "\t", comment.char = "", quote = "\"")
sch.deg <- read.table(synGet("syn10881886")@filePath, header = T, sep = "\t") %>% 
  separate(geneName, c("ens","Hugo_Gene"), sep = "\\|") %>%
  filter(cellLine1 == "HS01", cellLine2 == "HS11", treatment1 == "DMSO", treatment2 == "DMSO")

ncats.mush <- ncats.min %>% left_join(targs) %>% 
  group_by(Hugo_Gene) %>%
  summarize(mean = mean(diffr), sd = sd(diffr)) %>% 
  ungroup() %>% 
  inner_join(sch.deg)

library(ggrepel)
library(ggplot2)
ggplot(ncats.mush, aes(y = mean, x = logFC)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), color = "black") +
  geom_label_repel(data = ncats.mush %>% filter(abs(mean)>20), aes(label = Hugo_Gene))

ggplot(ncats.mush, aes(y = mean, x = logFC)) +
  geom_point() +
  labs(x = "gene logFC", y = "mean differential response")
