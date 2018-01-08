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
  targets <- left_join(sims2, evo) %>% dplyr::select(Common_Name, Original_molecule_SMILES, Structure_ID, `Tanimoto Similarity`) %>% distinct()
}



funct<-function(input, sim.thres) {
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
}


drug.resp <- readRDS("drugresp.rds")
ctrp.structures <- readRDS("ctrpstructures.rds")
fp.ctrp <- readRDS("fpctrp.rds")

plotSimCTRPDrugs <- function(input) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.ctrp, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  
  sims <- ldply(sims)
  sims2 <- sims %>% arrange(-sim)
  sims2$cpd_smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, ctrp.structures) %>% dplyr::select(makenames, cpd_name, `Tanimoto Similarity`) %>% distinct()
  
  top_drug <- top_n(drugs, 1, `Tanimoto Similarity`)
  
  drug.resp.single <- drug.resp[[top_drug$makenames]]

    cors<-sapply(colnames(drug.resp), function(x){
    test <- data.frame(drug.resp.single, drug.resp[[x]])
    if(nrow(test[complete.cases(test),])>1){
      cor<-cor.test(drug.resp.single, drug.resp[[x]], method = "spearman", use = "complete.obs")
      res <- c("p.val" = cor$p.value, cor$estimate)
    }else{
      res <- c("p.val" = -1, "rho" = 0)
    }
  })

  cors <- cors %>% 
    t() %>%
    as.data.frame() %>% 
    rownames_to_column("makenames") %>%
    inner_join(drugs) %>% 
    filter(p.val != -1)
  
  cors$Correlation <- cors$rho

  cors$`BH adj p.val` <- p.adjust(cors$p.val, method = "BH")
  cors
}

##same as above but for prepared sanger data

drug.resp.sang<-readRDS("drugresp_sang.rds")
sang.structures<-readRDS("sangstructures.rds")
fp.sang<-readRDS("fpsang.rds")


plotSimSangDrugs <- function(input) {
  input <- input
  fp.inp <- parseInputFingerprint(input)
  
  sims <- lapply(fp.inp, function(i) {
    sim <- sapply(fp.sang, function(j) {
      distance(i, j)
    })
    bar <- as.data.frame(sim)
    bar$match <- rownames(bar)
    bar
  })
  
  sims <- ldply(sims)
  sims2 <- sims %>% arrange(-sim)
  sims2$smiles <- as.character(sims2$match)
  sims2$`Tanimoto Similarity` <- signif(sims2$sim, 3)
  drugs <- left_join(sims2, sang.structures) %>% dplyr::select(makenames, sanger_names, `Tanimoto Similarity`) %>% distinct()
  
  top_drug <- top_n(drugs, 1, `Tanimoto Similarity`)
  
  drug.resp.single <- drug.resp.sang[[top_drug$makenames]]
    
  cors<-sapply(colnames(drug.resp.sang), function(x){
    test <- data.frame(drug.resp.single, drug.resp.sang[[x]])
    if(nrow(test[complete.cases(test),])>1){
      cor<-cor.test(drug.resp.single, drug.resp.sang[[x]], method = "spearman", use = "complete.obs")
      res <- c("p.val" = cor$p.value, cor$estimate)
    }else{
      res <- c("p.val" = -1, "rho" = 0)
    }
  })
  
  cors <- cors %>% 
    t() %>%
    as.data.frame() %>% 
    rownames_to_column("makenames") %>%
    inner_join(drugs) %>% 
    filter(p.val != -1)
  
  cors$Correlation <- cors$rho
  
  cors$`BH adj p.val` <- p.adjust(cors$p.val, method = "BH")
  cors
}