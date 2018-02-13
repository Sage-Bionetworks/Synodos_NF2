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

library(superheat)
ncats.mat <- fp.sim.matrix(fp.ncats) %>% as.data.frame()
ncats.mat$smiles2 <- names(fp.ncats)
ncats.mat <- left_join(ncats.mat, ncats.min)

superheat(ncats.mat %>% select(1:1911), 
          pretty.order.rows = TRUE,
          pretty.order.cols = TRUE,
          yr = ncats.mat$diffr)
