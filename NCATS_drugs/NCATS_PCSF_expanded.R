library(plyr)
library(tidyverse)
library(magrittr)
library(synapseClient)
synapseLogin()
source("NCATS_helpers.R")

map<- read.table("ncats_drugs_curated.txt", sep = "\t", quote = "", comment.char = "", header = T) %>% 
  select(ncgc, smiles2)

map2<- read.table("ncats_drugs_curated.txt", sep = "\t", quote = "", comment.char = "", header = T)

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

##create basis network adding in molecules from evotec/dgidb shiny database

evo <- readRDS("evotec_dgidb.rds") %>% 
  select(Common_Name, Original_molecule_SMILES) %>%
  distinct() %>% 
  set_colnames(c("ncgc", "smiles2")) %>% 
  bind_rows(map)

fp <- parseInputFingerprint(as.character(evo$smiles2))
sim <- fp.sim.matrix(fp)

colnames(sim) <- evo$ncgc
rownames(sim) <- evo$ncgc

library(PCSF)
library(igraph)

sim.filt <- sim
sim.net <- graph_from_adjacency_matrix(sim, weighted = "tanimoto", mode = "upper", diag = F)
sim.deg <- as.data.frame(degree(sim.net))

sim.filt[sim.filt<0.5] <- 0
sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = "tanimoto", mode = "upper", diag = F)
sim.05 <- as.data.frame(degree(sim.net))
sim.05.filt <- filter(sim.05, `degree(sim.net)` > 0)

sim.filt[sim.filt<0.7] <- 0
sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = "tanimoto", mode = "upper", diag = F)
sim.07 <- as.data.frame(degree(sim.net))
sim.07.filt <- filter(sim.07, `degree(sim.net)` > 0)

sim.filt[sim.filt<0.9] <- 0
sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = "tanimoto", mode = "upper", diag = F)
sim.09 <- as.data.frame(degree(sim.net))
sim.09.filt <- filter(sim.09, `degree(sim.net)` > 0)


ggplot() +
  geom_density(data= sim.05.filt, aes(x = `degree(sim.net)`, y = ..count..), fill = "#B5E07D",alpha = 0.5) +
  geom_density(data= sim.07.filt, aes(x = `degree(sim.net)`, y = ..count..), fill = "#377896",alpha = 0.1) +
  geom_density(data= sim.09.filt, aes(x = `degree(sim.net)`, y = ..count..), fill = "#08415C",alpha = 0.5) +
    scale_y_continuous(trans = "log1p")

sim.filt <- sim
sim.filt[sim.filt<0.5] <- 0
sim.filt <- sim.filt[,!is.na(colnames(sim.filt))]
sim.filt <- sim.filt[!is.na(rownames(sim.filt)),]

sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = T, mode = "upper", diag = F)
sim.05 <- as.data.frame(degree(sim.net))


#sim.test <- abs(sim.filt)
#sim.net <- graph_from_adjacency_matrix(sim.test, weighted = "tanimoto", mode = "upper", diag = F)

ggplot() +
  geom_density(data= sim.05, aes(x = `degree(sim.net)`, y = ..count..), fill = "#B5E07D",alpha = 0.5) +
  scale_y_continuous(trans = "log1p")

ncat.up <- filter(ncats.min, diffr < 0)
ncat.up$diffr <- abs(ncat.up$diffr)
terminals <- ncat.up$diffr
names(terminals) <- ncat.up$ncgc

subnet.ncats <- PCSF(sim.net, terminals, w = 2, b = 1, mu = 5e-04) ##all points superimposed
subnet.ncats <- PCSF(sim.net, terminals, w = 20, b = 1, mu = 5e-04) ##all points superimposed
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 1, mu = 5e-04) ##nice looking network with fair amount of Steiners
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 10, mu = 5e-04) ##all points superimposed
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 100, mu = 5e-04) ##all points superimposed
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 0.1, mu = 5e-04) ##nice looking network with fair amount of Steiners, less dense than previous (more stringy)
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 0.1, mu = 5e-03) ##nice looking network with fair amount of Steiners, more dense than previous?
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 0.1, mu = 5e-02) ##nice looking network with fair amount of Steiners, much less dense than previous 
subnet.ncats <- PCSF(sim.net, terminals, w = 200, b = 0.1, mu = 5e-01) ##very stringy, few steiners
subnet.ncats <- PCSF_rand(sim.net, terminals, w = 200, b = 0.1, mu = 5e-02) ##pretty dense network, try to increase nubmer of trees
subnet.ncats <- PCSF_rand(sim.net, terminals, w = 300, b = 0.1, mu = 5e-02) ##pretty dense network, try to increase nubmer of trees
subnet.ncats <- PCSF_rand(sim.net, terminals, w = 1000, b = 0.1, mu = 5e-02) ##pretty dense network, try to increase nubmer of trees
subnet.ncats <- PCSF_rand(sim.net, terminals, w = 10000, b = 0.1, mu = 5e-02) ##no subnetwork identified
plot(subnet.ncats)

subnet.ncats <- PCSF_rand(sim.net, terminals, w = 200, b = 0.1, mu = 5e-02) ##nice looking network with fair amount of Steiners, much less dense than previous 
plot(subnet.ncats)

output <- as_data_frame(subnet.ncats) 

evo <- readRDS("evotec_dgidb.rds") %>% 
  select(Common_Name, Original_molecule_SMILES) %>%
  distinct()

names <- evo$Common_Name

evo <- evo %>% add_column("drug_name" = names) %>% 
  set_colnames(c("ncgc","smiles2","drug_name")) %>% 
  bind_rows(map2) %>% 
  full_join(ncat.up)

write.table(evo, 'ncats_nodes_PCSF.txt', sep = "\t", quote = F, row.names = F)
write.table(output, "ncats_edges_PCSF.txt", sep = "\t", quote = F, row.names = F)

##### selective for NF2 wt over NF2 null
sim.filt <- sim
sim.filt[sim.filt<0.5] <- 0
sim.filt <- sim.filt[,!is.na(colnames(sim.filt))]
sim.filt <- sim.filt[!is.na(rownames(sim.filt)),]

sim.net <- graph_from_adjacency_matrix(sim.filt, weighted = T, mode = "upper", diag = F)

ncat.down <- filter(ncats.min, diffr > 0)
ncat.down$diffr <- abs(ncat.down$diffr)
terminals <- ncat.down$diffr
names(terminals) <- ncat.down$ncgc

subnet.ncats <- PCSF_rand(sim.net, terminals, w = 200, b = 0.1, mu = 5e-02) ##nice looking network with fair amount of Steiners, much less dense than previous 
plot(subnet.ncats)
output <- as_data_frame(subnet.ncats) 
write.table(output, "ncats_edges_PCSF_nf2WT.txt", sep = "\t", quote = F, row.names = F)
