library(synapseClient)
source("NCATS_helpers.R")
library(tidyverse)
library(magrittr)
synapseLogin()

map <- read.table(synGet("syn11559906")@filePath, sep = "\t", header = T, comment.char = "", quote = "")
ncats <-read.table(synGet("syn8314523")@filePath, header = T, sep = "\t")

fp <- parseInputFingerprint(as.character(map$smiles2))

sim <- fp.sim.matrix(fp)

colnames(sim) <- map$ncgc
rownames(sim) <- map$ncgc

ncats.null <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS01" | Cell.line == "Syn5") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgNull = mean(AUC))

ncats.wt <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS12" | Cell.line == "Syn1") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgWT = mean(AUC))

ncats.min <- full_join(ncats.null, ncats.wt) %>% 
  dplyr::mutate(diffr = avgNull - avgWT) %>% 
  set_colnames(c("ncgc", "avgNull", "avgWT", "diffr")) %>% 
  left_join(map)

sim2 <- sim %>% as.data.frame() %>% rownames_to_column(".id") %>% 
  gather("ncgc", "sim", -.id) %>% filter(sim > 0.5, .id != ncgc)

for (i in 1:nrow(sim2)){
  sim2[i, 1:2] = sort(sim2[i, 1:2])  ##eliminate duplicates from pairwise analysis
}
sim2 <- sim2 %>% distinct()
                        
write.table(ncats.min, "ncats_cytoscape_nodes.txt", sep = "\t", quote = F, row.names = F)
write.table(sim2, "ncats_cytoscape_edges.txt", sep = "\t", quote = F, row.names = F)

##use cytoscape to render these files

sim2 <- sim %>% as.data.frame() %>% rownames_to_column(".id") %>% 
  gather("ncgc", "sim", -.id) %>% filter(sim > 0.3, .id != ncgc)

for (i in 1:nrow(sim2)){
  sim2[i, 1:2] = sort(sim2[i, 1:2])  ##eliminate duplicates from pairwise analysis
}

sim2 <- sim2 %>% distinct()

write.table(sim2, "ncats_cytoscape_edges_sim_0pt3.txt", sep = "\t", quote = F, row.names = F)

## 

sim2 <- sim2 %>% filter(sim > 0.7, .id != ncgc)
write.table(sim2, "ncats_cytoscape_edges_sim_0pt7.txt", sep = "\t", quote = F, row.names = F)


library(PCSF)
library(igraph)

sim.filt <- sim
sim.filt[sim.filt<0.5] <- 0
sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = T, mode = "upper", diag = F)

ncat.up <- filter(ncats.min, diffr < 0)
ncat.up$diffr <- abs(ncat.up$diffr)
terminals <- ncat.up$diffr
names(terminals) <- ncat.up$ncgc

subnet.ncats <- PCSF_rand(sim.net, terminals, w = 1, b = 10, mu = 0.005) ##bit of hairball
plot(subnet.ncats)

ncats.pcsf <- as_data_frame(subnet.ncats)
write.table(ncats.pcsf, "PCSF_ncats_net.txt", sep = "\t", quote = F, row.names = F)

library(plyr)
grps<-cluster_edge_betweenness(subnet.ncats)
grps<-ldply(membership(grps))
colnames(grps) <- c("ncgc", "community")

library(randomcoloR)
pal <- as.data.frame(distinctColorPalette(108)) %>% 
  rownames_to_column("community") %>% 
  set_colnames(c("community", "community_hex"))

pal$community <- as.numeric(pal$community)

ncats.min.pcsf <- left_join(ncats.min, grps) %>% left_join(pal)

write.table(ncats.min.pcsf, "ncats_cytoscape_nodes_pcsf.txt", sep = "\t", quote = F, row.names = F)

library(Rtsne.multicore)
sim.dupes.rm <- sim[!duplicated.matrix(sim),!duplicated.matrix(sim)]
# tsne <- Rtsne.multicore(sim.dupes.rm,
#                         num_threads = 8, 
#                         perplexity = 100, 
#                         theta = 0, 
#                         max_iter = 1000,
#                         verbose = T) ##somewhat spread out
# 
# tsne <- Rtsne.multicore(sim.dupes.rm,
#                         num_threads = 8, 
#                         perplexity = 10, 
#                         theta = 0, 
#                         max_iter = 1000,
#                         verbose = T) ## more dense clusters but still "maplike"
# 
# tsne <- Rtsne.multicore(sim.dupes.rm,
#                         num_threads = 8, 
#                         perplexity = 1, 
#                         theta = 0, 
#                         max_iter = 1000,
#                         verbose = T) ## very dense clusters but probably too low of a perplexity - MANY clusters, looks almost like very high perplexity
# 
# tsne <- Rtsne.multicore(sim.dupes.rm,
#                         num_threads = 8, 
#                         perplexity = 50, 
#                         theta = 0, 
#                         max_iter = 1000,
#                         verbose = T) ## a large aggregate at low V1/V2
# 
# tsne <- Rtsne.multicore(sim.dupes.rm,
#                         num_threads = 8, 
#                         perplexity = 70, 
#                         theta = 0, 
#                         max_iter = 1000,
#                         verbose = T) ## one large aggregate at low V1/V2, might be non ideal perplexity?

tsne <- Rtsne.multicore(sim.dupes.rm,
                        num_threads = 8, 
                        perplexity = 25, 
                        theta = 0, 
                        max_iter = 1000,
                        verbose = T) ## seems fine!

res <- as.data.frame(tsne$Y)
res$ncgc <- colnames(sim.dupes.rm)

ncatannot <- ncats %>% dplyr::select(Sample.ID, Gene.Symbol) %>% distinct() %>% set_colnames(c("ncgc", "ncats.target"))

res.single <- left_join(res, ncatannot) %>% 
  add_column("true" = rep(TRUE, nrow(res))) %>% 
  filter(ncats.target != "") %>% 
  spread(ncats.target, true, fill = FALSE) ##add true/false columns for all annotated targets

targs4 <- read.table(synGet("syn11609840")@filePath,header = T)
res.poly <- inner_join(res, targs4) 

ggplot(data= res.single) +
  geom_point(aes(V1,V2)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.single) +
  geom_point(aes(V1,V2, color = CDK1)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.poly) +
  geom_point(aes(V1,V2, color = CDK1)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.single) +
  geom_point(aes(V1,V2, color = AURKA)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.poly) +
  geom_point(aes(V1,V2, color = AURKA)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.single) +
  geom_point(aes(V1,V2, color = CHEK1)) +
  labs(x = 'dimension 1', y = 'dimension 2')

ggplot(data= res.poly) +
  geom_point(aes(V1,V2, color = CHEK1)) +
  labs(x = 'dimension 1', y = 'dimension 2')


ggplot(data= res.single) +
  geom_hex(aes(V1,V2, fill = TUBB)) +
  guides(color = F) +
  scale_color_manual(values = c("TRUE" = "#179BED", "FALSE" = "#000000")) +
  theme_bw()
 
ggplot(data= res.poly) +
  geom_hex(aes(V1,V2, fill = TUBB)) +
  guides(color = F) +
  scale_color_manual(values = c("TRUE" = "#179BED", "FALSE" = "#000000")) +
  theme_bw()


ncats <-read.table(synGet("syn8314523")@filePath, header = T, sep = "\t")

ncats.null <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS01" | Cell.line == "Syn5") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgNull = mean(AUC))

ncats.wt <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line == "HS12" | Cell.line == "Syn1") %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgWT = mean(AUC))

ncats.min <- full_join(ncats.null, ncats.wt) %>% 
  dplyr::mutate(diffr = avgNull - avgWT) %>% 
  magrittr::set_colnames(c("ncgc", "avgNull", "avgWT", "diffr")) %>% 
  dplyr::left_join(map2)

res.single.2 <- left_join(res.single, ncats.min)
res.poly.2 <- left_join(res.poly, ncats.min)

ggplot(data= res.single.2) +
  geom_point(aes(V1,V2, color = -diffr)) +
  theme_bw()

ggplot(data= res.poly.2) +
  geom_point(aes(V1,V2, color = -diffr)) +
  theme_bw()
