options(java.parameters = "-Xmx250g" ) 
library(synapser)
library(plyr)
library(tidyverse)
library(rcdk)
library(fingerprint)
synLogin("allawayr","Warbler44")
this.file = "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/NCATS_drugs/2018-02-12/mapNCATStoDTE.R"

parser <- get.smiles.parser()

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

db <- readRDS(synGet("syn11712148")$path) %>% 
  filter(!is.na(hugo_gene)) %>% 
  select(internal_id, common_name, hugo_gene, mean_pchembl, n_quantitative, n_qualitative, known_selectivity_index, confidence)

fp.db <- readRDS(synGet("syn11808628")$path)
fp.db <- fp.db[names(fp.db) %in% unique(db$internal_id)]

library(parallel)
library(pbmcapply)

sims.temp <- pbmclapply(fp.ncats, function(i) {
  print(names(i))
  sim <- lapply(fp.db, function(j) {
    distance(i, j)
  })
  bar <- ldply(sim)
  colnames(bar) <- c("match", "similarity")
  bar <- filter(bar, similarity >= .95)
  bar
}, mc.cores = 64)


sims <- ldply(sims.temp)
colnames(sims) <- c("smiles2","internal_id", "sim")
sims <- sims %>% left_join(db) %>% left_join(ncats.structures)
sims.network <- sims %>% select(ncgc, hugo_gene) %>% distinct() %>%  mutate(one = 1) %>% spread(hugo_gene, one , fill = "0")

write.table(sims, "NCATS_to_DTE_tanimoto_cutoff_0_95.txt", sep = "\t", row.names = F)
synStore(File("NCATS_to_DTE_tanimoto_cutoff_0_95.txt", parentId = "syn11808773"), 
         used = c("syn11808628","syn11712148", "syn11559906"), 
         executed = this.file)

