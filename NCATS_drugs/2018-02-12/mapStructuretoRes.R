library(synapser)
library(plyr)
library(tidyverse)
library(rcdk)
library(fingerprint)
library(pbapply)
library(pheatmap)
synLogin()

is.smiles <- function(x, verbose = TRUE) { ##corrected version from webchem
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # x <- 'Clc(c(Cl)c(Cl)c1C(=O)O)c(Cl)c1Cl'
  if (length(x) > 1) {
    stop('Cannot handle multiple input strings.')
  }
  out <- try(rcdk::parse.smiles(x), silent = TRUE)
  if (inherits(out[[1]], "try-error") | is.null(out[[1]])) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

parseInputFingerprint <- function(input) {
  test_smiles <- is.smiles(input)
  if(is.smiles(input==TRUE)){
    input.mol <- parse.smiles(as.character(input))
    lapply(input.mol, do.typing)
    lapply(input.mol, do.aromaticity)
    lapply(input.mol, do.isotopes)
    fp.inp <- lapply(input.mol, get.fingerprint, type = "circular")
  }else{
    print('Please input a valid SMILES string.')
  }
}

ncats.structures <- read.table(synGet("syn11559906")$path, sep = "\t", comment.char = "", header = T, quote = "") %>% 
  filter(smiles2 != "")

fp.ncats <- sapply(ncats.structures$smiles2, parseInputFingerprint)
names(fp.ncats) <- ncats.structures$ncgc

ncats <- read.table(synGet("syn8314523")$path, sep = "\t", comment.char = "", header = T, quote = "")  

ncats.null <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line %in% c("HS01","Syn5", "Syn6")) %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgNull = mean(AUC))

ncats.wt <- ncats %>% 
  select(Sample.ID, Cell.line, AUC) %>% 
  filter(Cell.line %in% c("HS11","Syn1")) %>% 
  select(-Cell.line) %>% 
  dplyr::group_by(Sample.ID) %>% 
  dplyr::summarise(avgWT = mean(AUC))

ncats.min <- full_join(ncats.null, ncats.wt) %>% 
  dplyr::mutate(diffr = avgNull - avgWT) %>% 
  magrittr::set_colnames(c("ncgc", "avgNull", "avgWT", "diffr")) %>% 
  dplyr::inner_join(ncats.structures)

library(tibble)
ncats.mat <- fp.sim.matrix(fp.ncats) %>% as.data.frame() %>% set_names(names(fp.ncats))
ncats.mat$input <- names(fp.ncats)
ncats.tidy <- ncats.mat %>% gather("out", "connectivity", -input) %>% distinct()

ncats.targets <- read.table(synGet("syn11808774")$path, sep = "\t", header=T) %>% 
  select(ncgc, hugo_gene) %>% 
  mutate(connectivity = 1) %>% 
  set_names(c('input', "out", "connectivity")) %>% 
  filter(!is.na(out))

ncats.targets.reverse <- ncats.targets %>% select(out, input, connectivity) %>% 
  set_names(c('input', "out", "connectivity"))

ncats.graph <- bind_rows(ncats.tidy, ncats.targets) %>% 
  bind_rows(ncats.targets.reverse) %>% 
  distinct() %>% 
  spread(out, connectivity, fill = 0) %>% 
  remove_rownames() %>% 
  column_to_rownames("input") %>% 
  as.matrix()

ncats.min <- ncats.min %>% select(ncgc, diffr) %>% filter(ncgc %in% ncats.targets$input)
ncats.min$diffr[ncats.min$diffr >0] <- 0 
ncats.min$diffr <- abs(ncats.min$diffr)

names <- as.data.frame(rownames(ncats.graph)) %>% set_names("ncgc")
ncats.min <- full_join(names,ncats.min)
ncats.min$diffr[is.na(ncats.min$diffr)] <- 0

p0 <- ncats.min$diffr
p.norm <- (p0/sum(p0))

walk <- random.walk(p.norm, ncats.graph) %>% as.data.frame() %>%  set_names(c("walk"))
walk$names <- names$ncgc

heat <- heat.diffusion(p.norm, ncats.graph) %>% as.data.frame() %>% set_names(c("heat"))
heat$names <- names$ncgc

perms <- full_join(heat, walk)

freq <- bind_rows(ncats.tidy, ncats.targets) %>% 
  bind_rows(ncats.targets.reverse) %>% 
  distinct() %>% 
  group_by(input) %>% 
  tally() %>% 
  ungroup() %>% 
  select(input, n) %>% 
  set_names(c("names", "n"))

perms <- perms[grep("NCGC.+", invert = T, perms$names),]
perms <- left_join(perms,freq)
  
library(ggplot2)
library(ggrepel)
ggplot(perms, aes(x = walk, y = heat)) +
  geom_point() +
  geom_label_repel(data = perms %>% top_n(20, walk), aes(label = names, fill = log(n)))

cor(perms$heat, perms$n)
cor(perms$walk, perms$n)

perms.norm <- mutate(perms, walk = walk/n, heat = heat/n)

ggplot(perms.norm, aes(x = walk, y = heat)) +
  geom_point() +
  geom_label_repel(data = perms.norm %>% top_n(20, heat), aes(label = names))

