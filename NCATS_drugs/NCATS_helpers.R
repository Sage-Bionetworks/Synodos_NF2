library(rJava)
library(rcdk)
library(fingerprint)
library(webchem)
library(plyr)

##converts SMILES string to fingerprint
parseInputFingerprint <- function(input) {
  input.mol <- parse.smiles(input)
  fp.inp <- lapply(input.mol, get.fingerprint, type = "extended")
  
}

convertDrugToSmiles <- function(input) {
  filt <- filter(syns, Common_Name == input) %>% dplyr::select(Original_molecule_SMILES)
  filt
}

getTargetList <- function(selectdrugs) {
  targets <- filter(evo, Common_Name %in% selectdrugs) %>% dplyr::select(Common_Name, Hugo_Gene, MedianActivity_nM, N_quantitative, N_qualitative, 
                                                                         N_inactive, N_DGIDB, Confidence_Score) %>% arrange(-N_quantitative)
  if (nrow(targets) > 1) {
    targets
  } else {
    print("none found")
  }
}

getSimMols <- function(input, sim.thres) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.evo, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  sims <- ldply(sims)
  sims2 <- sims %>% filter(sim >= sim.thres) %>% arrange(-sim)
  sims2$Original_molecule_SMILES <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  targets <- left_join(sims2, evo) %>% dplyr::select(Common_Name, `Tanimoto Similarity`) %>% distinct()
}
