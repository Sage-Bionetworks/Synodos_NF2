library(tidyverse)
library(synapser)
library(rcdk)
library(fingerprint)
library(diffusr)
library(plyr)
library(matrixStats)
library(pbapply)
library(enrichR)
synLogin()
### drug data

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

is.whole <- function(a) { 
  (is.numeric(a) && floor(a)==a) ||
    (is.complex(a) && floor(Re(a)) == Re(a) && floor(Im(a)) == Im(a))
}


targets <- read.table(synGet("syn12063782")$path, comment.char = "", header = T, sep = "\t") %>% filter(sim >0.95)

target.network <- select(targets, ncgc, hugo_gene, confidence)

fp.ncats <- sapply(unique(targets$smiles), parseInputFingerprint)
names(fp.ncats) <- unique(targets$ncgc)

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
  dplyr::inner_join(target.network)

ncats.mat <- fp.sim.matrix(fp.ncats) %>% as.data.frame() %>% set_names(names(fp.ncats))
ncats.mat$input <- names(fp.ncats)
ncats.tidy <- ncats.mat %>% gather("out", "connectivity", -input) %>% distinct()

ncats.targets <- read.table(synGet("syn12063782")$path, sep = "\t", header=T) %>% 
  filter(ncgc %in% names(fp.ncats)) %>% 
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

ncats.min <- ncats.min %>% select(ncgc, diffr) %>% 
  filter(ncgc %in% ncats.targets$input) %>% 
  distinct()
ncats.min$diffr[ncats.min$diffr >0] <- 0 
ncats.min$diffr <- abs(ncats.min$diffr)

names <- as.data.frame(rownames(ncats.graph)) %>% set_names("ncgc")
ncats.min <- full_join(names,ncats.min)
ncats.min$diffr[is.na(ncats.min$diffr)] <- 0


###Random walk  analysis

 p0 <- ncats.min$diffr
 walk <- random.walk(p0, ncats.graph) %>% as.data.frame() %>%  set_names(c("walk"))
 walk$name <- names$ncgc

foo <- ncats.min

res.walk <- pbsapply(1:1000, function(x){
  foo$diffr[grep("NCGC.+",foo$ncgc)] <- sample(foo$diffr[grep("NCGC.+",foo$ncgc)])
  p0 <- foo$diffr
  walk <- random.walk(p0, ncats.graph, r = 0.5) %>% as.data.frame() %>% set_names(c("heat"))
})

res2<-ldply(res.walk) %>% select(-.id) %>% t() %>% as.data.frame()

res2$name <- c(foo$ncgc)
res3 <- res2
res3$mean <- rowMeans(res3[,1:1000]) 
res3$sd <- rowSds(as.matrix(res3[,1:1000]))
res3 <- res3 %>% dplyr::select(1001:1003) %>% 
  inner_join(walk) %>% 
  mutate(diff = walk-mean) %>% 
  mutate(test = diff>sd) 

db <- listEnrichrDbs()[,1]
db <- db[!db %in% "NCI-Nature_2015"]
res4 <- res3[!grepl("NCGC.+",res3$name),]
res4 <- res4 %>% filter(diff>0 & test==TRUE) %>% 
  top_n(50, diff)

enrichment <- enrichr(res4$name, db)
enrichment1 <- ldply(enrichment)


igraph <- igraph::graph_from_adjacency_matrix(ncats.graph)

library(dnet)

foo2 <- column_to_rownames(ncats.min, "ncgc")
rw <- dRWR(igraph, setSeeds = foo2, restart = 0.5)
rw <- as.matrix(rw) %>% as.data.frame()
rw$name <- rownames(foo2)


join <- full_join(res3, rw)

library(ggplot2)
library(plotly)

p <- ggplot(join %>% filter(!grepl("NCGC", name))) +
  geom_point(aes(x = walk, y = V1, color = V1+walk, label = name))
ggplotly(p)




