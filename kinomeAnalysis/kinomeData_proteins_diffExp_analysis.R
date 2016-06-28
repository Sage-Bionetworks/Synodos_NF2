library(synapseClient)
synapseLogin()
options(stringsAsFactors = F)

library(reshape2)
library('ggplot2')
library('data.table')
library('plyr')
library('dplyr')
library('tidyr')
library('reshape2')
library('gdata')
library('rbokeh')
library('NMF')
library(doMC)
registerDoMC(3)

source("~/dev/apRs/plotting_helpers.R")

get_overlapping_proteins_across_replicates <- function(df){
  gene_lists <- lapply(unique(df$replicate), function(x)  { filter(df,replicate == x) %>% select(Gene) %>% .$Gene })
  gene_lists <- lapply(gene_lists, unique)
  common_genes <- Reduce(intersect, gene_lists)
}

get_diffExp_proteins <- function(df){
  proteins_cond1 <- kinome_protein_data %>% filter(cellLine == df$cellLine1 & time == df$time1 & drug == df$drug)
  proteins_cond2 <- kinome_protein_data %>% filter(cellLine == df$cellLine2 & time == df$time2 & drug == df$drug)
  #allow only common proteins across replicates
  common_prot_rep_cond1 <- get_overlapping_proteins_across_replicates(proteins_cond1)
  common_prot_rep_cond2 <- get_overlapping_proteins_across_replicates(proteins_cond2)
  proteins_cond1 <- proteins_cond1 %>% filter(Gene %in% common_prot_rep_cond1)
  proteins_cond2 <- proteins_cond2 %>% filter(Gene %in% common_prot_rep_cond2)
  #taking a union of all proteins/genes 
  all_proteins <- union(proteins_cond1$Gene, proteins_cond2$Gene)
  #get diff exp statistic for each protein
  proteins_pvals <- ldply(all_proteins, .fun = get_protein_pval, .parallel = T, df)
  #adjust p-val
  proteins_pvals$pval_adj <-  p.adjust( proteins_pvals$pval, method='fdr')
  proteins_pvals
}

get_protein_pval <- function(protein, df){
  cond1_peptides <- filter(kinome_peptide_data, Gene == protein &
                             cellLine == df$cellLine1 & time == df$time1 & drug == df$drug)
  cond2_peptides <- filter(kinome_peptide_data, Gene == protein &
                             cellLine == df$cellLine2 & time == df$time2 & drug == df$drug)
  if(nrow(cond1_peptides) < 2 || nrow(cond2_peptides) < 2){ pval= NA 
  } else { 
    pval <- wilcox.test(cond2_peptides$log2NormRatio ,cond1_peptides$log2NormRatio)$p.value 
  }
  result <- data.frame(pval=pval, peptides_cond1 = nrow(cond1_peptides), peptides_cond2 = nrow(cond2_peptides),
                       medRatio_peptides_cond1 = median(cond1_peptides$log2NormRatio),
                       medRatio_peptides_cond2 = median(cond2_peptides$log2NormRatio),
                       meanRatio_peptides_cond1 = mean(cond1_peptides$log2NormRatio),
                       meanRatio_peptides_cond2 = mean(cond2_peptides$log2NormRatio),
                       se_peptides_cond1 = sd(cond1_peptides$log2NormRatio)/ sqrt(length(cond1_peptides$log2NormRatio)),
                       se_peptides_cond2 = sd(cond2_peptides$log2NormRatio)/ sqrt(length(cond2_peptides$log2NormRatio)),
                       drug = df$drug, protein = protein)
  #also add the average ratio at the protein level
  avgRatio_protein_cond1 <- filter(kinome_protein_data, Gene == protein & 
                                     cellLine == df$cellLine1 & time == df$time1 & drug == df$drug ) %>% summarise(m = mean(log2NormRatio)) %>% .$m
  avgRatio_protein_cond2 <- filter(kinome_protein_data, Gene == protein & 
                                     cellLine == df$cellLine2 & time == df$time2 & drug == df$drug ) %>% summarise(m = mean(log2NormRatio)) %>% .$m
  result['avgRatio_proteinRep_cond1'] = avgRatio_protein_cond1
  result['avgRatio_proteinRep_cond2'] = avgRatio_protein_cond2
  result
}


#############
# MAIN
#############
#get all the kinome data
kinome_protein_synid <- 'syn4951080'
kinome_protein_data <- synGet(kinome_protein_synid)
kinome_protein_data <- read.table(kinome_protein_data@filePath, sep="\t", header=T, check.names = F)
kinome_protein_data$Gene <- toupper(kinome_protein_data$Gene)
reqd_cols <- c('Accession', 'Gene', 'Family', 'num_Proteins', 'num_Peptides',  'cellLine', 'replicate',  'drug',
               'time', 'ratio', 'count', 'log2ratio', 'log2NormRatio', 'uniq_peptides', 'condition')
kinome_protein_data <- kinome_protein_data[,reqd_cols]


kinome_peptide_synid <- 'syn4951148'
kinome_peptide_data <- synGet(kinome_peptide_synid)
kinome_peptide_data <- read.table(kinome_peptide_data@filePath, sep="\t", header=T, check.names = F)
kinome_peptide_data$fileId <- NULL
#keep only those peptides which for which proteins are present #will reduce the data size
kinome_peptide_data <- kinome_peptide_data %>% filter(Accession %in% kinome_protein_data$Accession)
#add in the Gene Name for mapping protein data to peptide data
kinome_peptide_data <- merge(kinome_peptide_data, kinome_protein_data[,c('Accession', 'Gene', 'cellLine', 'replicate', 'drug',  'time')], sort=F, all.x=T)


#removing replicate1 for cellLines MS03 and MS12
#kinome_peptide_data <- kinome_peptide_data %>% filter( !(cellLine %in% c('MS03', 'MS12') & replicate %in% c('run1')))
#kinome_protein_data <- kinome_protein_data %>% filter( !(cellLine %in% c('MS03', 'MS12') & replicate %in% c('run1')))

#create the data frame for differential expression comparisons
drugs = c('Pano', 'CUDC', 'GSK458')

#diff expressed proteins across NULL AND WT at 2H and 24hours
comparisons_24h = data.frame(cellLine1 = c('MS03','HS01', 'MS02', 'Syn5', 'Syn6'),
                             cellLine2 = c('MS12','HS11', 'MS12', 'Syn1', 'Syn1'),
                             time1 = '24h', time2='24h')
comparisons_2h  = data.frame(cellLine1 = c('MS03','HS01', 'MS02', 'Syn5', 'Syn6'),
                             cellLine2 = c('MS12','HS11', 'MS12', 'Syn1', 'Syn1'),
                              time1 = '2h', time2='2h')
comparisons_across_same_cellLines  = data.frame(cellLine1 = c('MS03', 'MS12', 'HS11', 'HS01', 'MS02', 'Syn1', 'Syn5', 'Syn6'), 
                                                cellLine2 = c('MS03', 'MS12', 'HS11', 'HS01', 'MS02', 'Syn1', 'Syn5', 'Syn6'), 
                                                time1 = '2h', time2='24h')
comparisons <- rbind(comparisons_2h, comparisons_24h, comparisons_across_same_cellLines)
comparisons <- merge(comparisons, data.frame(drug=drugs))

library(doMC)
registerDoMC(3)
diff_exp_proteins <- ddply(comparisons, 
                           .variables = c('cellLine1','cellLine2', 'time1', 'time2', 'drug'),
                           .fun = get_diffExp_proteins, .parallel=T)


#upload to synapse
outfile <- "Synodos_kinome_diffExp_kinases_data.tsv"
write.table(diff_exp_proteins, file=outfile, sep="\t", col.names=T, row.names = F)
synStore(File(outfile, parentId = 'syn4259360'),
         used = c(kinome_peptide_synid, kinome_protein_synid),
         executed = "https://github.com/Sage-Bionetworks/Synodos_NF2/blob/master/kinomeAnalysis/kinomeData_proteins_diffExp_analysis.R")
unlink(outfile)


