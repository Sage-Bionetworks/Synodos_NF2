library(synapseClient)
library(tidyverse)
library(magrittr)
synapseLogin()

ncats <-read.table(synGet("syn8314523")@filePath, header = T, sep = "\t")

this.file 
map <- ncats %>% 
  select(Sample.ID, Sample.Name) %>% 
  distinct() %>% 
  set_colnames(c("ncgc", "drug_name"))

write.table(map, "ncatsdrugs.txt", sep = "\t", quote = F, row.names = F) ##input both lists to pubchem CIR webtool 

names <- read.table("drugnamestosmiles.txt", sep = "\t", comment.char = "", quote = "")
names <- names %>% 
  set_colnames(c("drug_name", "smiles")) %>% 
  add_column(count = c(1:nrow(names))) %>% 
  group_by(drug_name) %>% 
  top_n(1, count) %>% 
  ungroup() %>% 
  select(drug_name, smiles) 

ncgc <- read.table("ncgctosmiles.txt", sep = "\t", comment.char = "")
ncgc <- ncgc %>% 
  set_colnames(c("ncgc", "smiles2")) %>% 
  add_column(count = c(1:nrow(ncgc))) %>% 
  group_by(ncgc) %>% 
  top_n(1, count) %>% 
  ungroup() %>% 
  select(ncgc, smiles2)

map <- full_join(map, names) 
map <- full_join(map, ncgc)

##sanity check to compare smiles from two sources
source("NCATS_helpers.R")

temp.map <- map %>% filter(smiles != " " & smiles != "" & smiles2 != "" & smiles2 != " ")
temp.map$sim <- 0

for(i in 1:nrow(temp.map)){

  x<-temp.map[i,3]
  y<-temp.map[i,4]

  x<-parseInputFingerprint(as.character(x))
  y<-parseInputFingerprint(as.character(y))

  d <- distance(x[[1]],y[[1]])

  temp.map[i,5] <- d
}

ggplot(data = temp.map) + 
  geom_histogram(aes(x=sim))

## in cases where smiles disagree, NCGC lookup seems to be the accurate one (cursory search on pubchem)
## use this for mapping 
map <- ncats %>% 
  select(Sample.ID, Sample.Name) %>% 
  distinct() %>% 
  set_colnames(c("ncgc", "drug_name"))

map <- full_join(map, ncgc) 

write.table(map, "ncats_drugs_curated.txt", sep = "\t", quote = F, row.names = F)

map <- read.table("ncats_drugs_curated.txt", sep = "\t", quote = "", comment.char = "", header = T)
synStore(File("ncats_drugs_curated.txt", parentId = "syn8682571"), executed = this.file)