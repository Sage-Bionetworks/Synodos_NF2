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
registerDoMC(4)

source("~/dev/apRs/plotting_helpers.R")

get_overlapping_proteins_across_replicates <- function(df){
  gene_lists <- lapply(unique(df$replicate), function(x)  { filter(df,replicate == x) %>% select(Gene) %>% .$Gene })
  gene_lists <- lapply(gene_lists, unique)
  common_genes <- Reduce(intersect, gene_lists)
}
get_diffExp_proteins <- function(df){
  replicate1 <- unlist(strsplit(df$replicate1, split=","))
  replicate2 <- unlist(strsplit(df$replicate2, split=","))
  proteins_cond1 <- kinome_protein_data %>% filter(cellLine == df$cellLine1 & time == df$time1 & replicate %in% replicate1 & drug == df$drug)
  proteins_cond2 <- kinome_protein_data %>% filter(cellLine == df$cellLine2 & time == df$time2 & replicate %in% replicate2 & drug == df$drug)
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
  replicate1 <- unlist(strsplit(df$replicate1, split=","))
  replicate2 <- unlist(strsplit(df$replicate2, split=","))
  cond1_peptides <- filter(kinome_peptide_data, Gene == protein &
                             cellLine == df$cellLine1 & time == df$time1 & drug == df$drug & replicate %in% replicate1)
  cond2_peptides <- filter(kinome_peptide_data, Gene == protein &
                             cellLine == df$cellLine2 & time == df$time2 & drug == df$drug & replicate %in% replicate2)
  if(nrow(cond1_peptides) < 2 || nrow(cond2_peptides) < 2){ pval= NA 
  } else { 
    pval <- wilcox.test(cond2_peptides$log2ratio ,cond1_peptides$log2ratio)$p.value 
  }
  result <- data.frame(pval=pval, peptides_cond1 = nrow(cond1_peptides), peptides_cond2 = nrow(cond2_peptides),
                       medRatio_peptides_cond1 = median(cond1_peptides$log2ratio),
                       medRatio_peptides_cond2 = median(cond2_peptides$log2ratio),
                       meanRatio_peptides_cond1 = mean(cond1_peptides$log2ratio),
                       meanRatio_peptides_cond2 = mean(cond2_peptides$log2ratio),
                       se_peptides_cond1 = sd(cond1_peptides$log2ratio)/ sqrt(length(cond1_peptides$log2ratio)),
                       se_peptides_cond2 = sd(cond2_peptides$log2ratio)/ sqrt(length(cond2_peptides$log2ratio)),
                       drug = df$drug, protein = protein)
  #also add the average ratio at the protein level
  avgRatio_protein_cond1 <- filter(kinome_protein_data, Gene == protein & 
                                     cellLine == df$cellLine1 & time == df$time1 & drug == df$drug & replicate %in% replicate1) %>% summarise(m = mean(log2ratio)) %>% .$m
  avgRatio_protein_cond2 <- filter(kinome_protein_data, Gene == protein & 
                                     cellLine == df$cellLine2 & time == df$time2 & drug == df$drug & replicate %in% replicate2) %>% summarise(m = mean(log2ratio)) %>% .$m
  result['avgRatio_proteinRep_cond1'] = avgRatio_protein_cond1
  result['avgRatio_proteinRep_cond2'] = avgRatio_protein_cond2
  result
}

get_waterfall_plot <- function(df, line1, line2, time1, time2, selected_drug, pval_adj = .05, title = NA){
  df <- filter(df, cellLine1 == line1 & cellLine2 == line2 & time1 == time1 & time2 == time2 & drug == selected_drug & pval_adj <= .05)
  direction <- factor(sign(df$medRatio_peptides_cond1) == sign(df$medRatio_peptides_cond2))
  direction <- revalue(direction, c('FALSE' = "opposite", 'TRUE' = 'same'))
  df['direction'] = direction
  protein_levels <- df %>% group_by(direction) %>% arrange(meanRatio_peptides_cond2) %>% .$protein
  df$protein <- factor(df$protein, levels=protein_levels)
  
  #select the kinase and mean peptide ratio
  df1 <- df %>% select(protein, meanRatio_peptides_cond1, meanRatio_peptides_cond2)
  df1 <- melt(df1, id.vars=c('protein'), variable.name = 'condition' , value.name = 'log2ratio')
  df1$condition <- revalue(df1$condition, c('meanRatio_peptides_cond1' = paste0(cellLine1,'_',selected_drug,'_',time1),
                                            'meanRatio_peptides_cond2' = paste0(cellLine2,'_',selected_drug,'_',time2)))
  #get the standard error
  df2 <- df %>% select(protein, se_peptides_cond1, se_peptides_cond2)
  df2 <- melt(df2, id.vars=c('protein'), variable.name = 'condition' , value.name = 'SE')
  df2$condition <- revalue(df2$condition, c('se_peptides_cond1' = paste0(cellLine1,'_',selected_drug,'_',time1),
                                            'se_peptides_cond2' = paste0(cellLine2,'_',selected_drug,'_',time2)))
  data <- merge(df1,df2)
  data['log2ratio_min'] = data$log2ratio - data$SE
  data['log2ratio_max'] = data$log2ratio + data$SE
  
  p <- ggplot(data=data,aes(x=protein, y=log2ratio, fill=condition)) + geom_bar(stat='identity',position="dodge") 
  p <- p + geom_errorbar(aes(ymax=log2ratio_max,ymin=log2ratio_min), width=.7, position="dodge", colour="grey50") 
  p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=10)) + ylab('log2 mean peptide ratio') 
  p <- p + xlab('kinases') +  theme(legend.position="top")
  #p <- p + ggtitle(paste('Differentially Expressed Kinases at', time1, 'vs', time2, cellLine1, 'vs', cellLine2 , 'drug: ', selected_drug)) 
  p <- p + ggtitle(paste(selected_drug,'-',title)) + theme(plot.title=element_text(family="Times", size=16, face="bold"))
  p <- p + coord_flip()
  p
}

plot_uniq_proteins_in_null <- function(null_line, wt_line, pval_adj = .05){
  results <- lapply(drugs, function(selected_drug) {
    #filter  
    null <- filter(diff_exp_proteins , pval_adj <= .05 & drug  == selected_drug & cellLine1 == null_line & cellLine2 == null_line )
    wt <- filter(diff_exp_proteins, pval_adj <= .05 & drug == selected_drug  & cellLine1 == wt_line & cellLine2 == wt_line )
    #only in NULL
    uniq_proteins_null <- setdiff(unique(null$protein), unique(wt$protein))
    #get the waterfall plot for kinase only in NULL
    null <- null %>% filter(protein %in% uniq_proteins_null)
    p <- get_waterfall_plot(null, null_line, null_line, '2h', '24h', selected_drug, title = 'unique to NULL')  
    p <- p + theme(axis.text.y=element_text(size=12))
    p
  })
  results
}


get_matrix_for_heatmap <- function(selected_drugs,cellLine1,cellLine2){
  # get signif differentially expressed proteins
  # a.) p-val adjusted
  # b.) get the proteins seen only in on condition and > 20 peptides at least
  # c.) get proteins for a choose set of drug/s
  sig_diff_exp_proteins <- diff_exp_proteins %>% filter( pval_adj < .10 | 
                                                           (is.na(pval) & (peptides_cond2 > 15 | peptides_cond1 > 15)) )
  sig_diff_exp_proteins <- sig_diff_exp_proteins %>% filter(drug %in% selected_drugs) %>% .$protein
  sig_diff_exp_proteins <- unique(sig_diff_exp_proteins)
  m <- diff_exp_proteins %>% filter(protein %in% sig_diff_exp_proteins & drug %in% selected_drugs)
  m$peptides_cond1 <- NULL
  m$peptides_cond2 <- NULL
  m$avgRatio_proteinRep_cond1 <- NULL
  m$avgRatio_proteinRep_cond2 <- NULL
  m <- melt(m, id.vars = c('pval', 'drug', 'pval_adj', 'protein', 'time'), variable.name = 'cellLine',
            value.name = 'log2ratio') 
  m <- m %>% transform(cellLine = revalue(cellLine, c('medRatio_peptides_cond1' = cellLine1, 'medRatio_peptides_cond2' = cellLine2)))
  m <- dcast(m,  protein ~ drug + cellLine + time, value.var = "log2ratio" )
  rownames(m) <- m$protein
  m$protein <- NULL
  #remove NA rows
  m <- m[apply(m, 1, function(x) sum(is.na(x)) == 0),]
  m <- t(scale(t(m)))
}





plots_prots_changing_dir_null_vs_wt <- function(null,wt,pval_adj = .05){
  plots <- lapply(drugs, function(selected_drug) {
    null1 <- filter(null, drug  == selected_drug)
    wt1 <- filter(wt, drug == selected_drug)
    #common kinase null and wt    
    common_proteins <- intersect(unique(null1$protein), unique(wt1$protein))
    #get the direction of protein expression in condition2 w.r.t condition1
    null1 <- null1 %>% mutate(direction = as.character(factor(sign(meanRatio_peptides_cond2 - medRatio_peptides_cond1)))) %>%
      mutate(direction = revalue(direction, c('-1' = 'down', '1' = 'up')) ) %>% filter(protein %in% common_proteins)
    null1 <- null1[match(common_proteins, null1$protein),]
    
    wt1 <- wt1 %>% mutate(direction = as.character(factor(sign(meanRatio_peptides_cond2 - medRatio_peptides_cond1)))) %>%
      mutate(direction = revalue(direction, c('-1' = 'down', '1' = 'up')) ) %>% filter(protein %in% common_proteins) %>% arrange(common_proteins)
    wt1 <- wt1[match(common_proteins, wt1$protein),]
    diff_direction_prots <- !wt1$direction  == null1$direction
    wt1 <- wt1[diff_direction_prots,]
    null1 <- null1[diff_direction_prots,]
    df = data.frame('WT_2h' = wt1$meanRatio_peptides_cond1, 'WT_24h' = wt1$meanRatio_peptides_cond2,
                    'NULL_2h' = null1$meanRatio_peptides_cond1, 'NULL_24h' = null1$meanRatio_peptides_cond2,
                    'WT_pval' = wt1$pval_adj, 'NULL_pval' = null1$pval_adj, 'protein' = null1$protein)
    #kinase changing direction
    x <- filter(df, WT_pval <= pval_adj | NULL_pval <= pval_adj) %>% select(protein,WT_2h,WT_24h,NULL_2h,NULL_24h) %>%
      melt(id.vars=c('protein')) %>% separate(variable, into=c('cellLine', 'time'))
    p <- ggplot(data=x, aes(x=factor(time,levels = c('2h', '24h')), y=value, color=cellLine)) + geom_point() + geom_line(aes(group=cellLine)) + facet_wrap( ~ protein )
    p <- p + ylab('log2ratio(mean)') + xlab('time') + ggtitle(paste(selected_drug,'-', 'kinases changing direction (FDR <= ', pval_adj, ' )'))
    p <- p + theme_bw()
    p
  })
  plots
}


#############
# MAIN
#############
#get all the kinome data
kinome_protein_synid <- 'syn4259365'
kinome_protein_data <- synGet(kinome_protein_synid)
kinome_protein_data <- read.table(kinome_protein_data@filePath, sep="\t", header=T, check.names = F)
kinome_protein_data$Gene <- toupper(kinome_protein_data$Gene)
reqd_cols <- c('Accession', 'Gene', 'Family', '# Proteins', '# Peptides',  'cellLine', 'replicate',  'drug',
               'time', 'ratio', 'count', 'log2ratio', 'uniq_peptides', 'condition')
kinome_protein_data <- kinome_protein_data[,reqd_cols]

kinome_peptide_synid <- 'syn4378913'
kinome_peptide_data <- synGet(kinome_peptide_synid)
kinome_peptide_data <- read.table(kinome_peptide_data@filePath, sep="\t", header=T, check.names = F)
kinome_peptide_data$fileId <- NULL
#keep only those peptides which for which proteins are present #will reduce the data size
kinome_peptide_data <- kinome_peptide_data %>% filter(Accession %in% kinome_protein_data$Accession)
#add in the Gene Name for mapping protein data to peptide data
kinome_peptide_data <- merge(kinome_peptide_data, kinome_protein_data[,c('Accession', 'Gene', 'cellLine', 'replicate', 'drug',  'time')], sort=F, all.x=T)

#create the data frame for comparisons
drugs = c('Pano', 'CUDC', 'GSK458')
data = data.frame(cellLine1 = c('MS03', 'MS12', 'HS11', 'HS01', 'MS02'), cellLine2 = c('MS03', 'MS12', 'HS11', 'HS01', 'MS02'), time1 = '2h', time2='24h', 
                  replicate1 = c('run2,run3', 'run2,run3', 'run1,run2,run3', 'run1,run2,run3', 'run2,run3'), 
                  replicate2 = c('run2,run3', 'run2,run3', 'run1,run2,run3', 'run1,run2,run3', 'run2,run3'))
temp_df <- expand.grid(data$cellLine1, drugs)
colnames(temp_df) <- c('cellLine1', 'drug')
data <- merge(data, temp_df)
library(doMC)
registerDoMC(12)
diff_exp_proteins <- ddply(data, .variables = c('cellLine1','cellLine2', 'time1', 'time2', 'replicate1', 'replicate2', 'drug'),
                           .fun = get_diffExp_proteins, .parallel=T)


#upload to synapse
outfile <- "kinome_data.tsv"
write.table(diff_exp_proteins, file=outfile, sep="\t", col.names=T, row.names = F)
synStore(File(outfile, parentId = 'syn4259360'),
         used = c(kinome_peptide_synid, kinome_protein_synid),
         executed = "")


#####################
# CASE 1: compare isogenic lines HS01(NULL) & HS11(WT)
#####################
#HS01 / NULL line
drug_plots_HS01_24h_vs_2h <- lapply(drugs, function(x){
  get_waterfall_plot(diff_exp_proteins, line1='HS01', line2='HS01', time1='2h', time2='24h', selected_drug=x, title = "NULL" )
})
#HS11 / WT line
drug_plots_HS11_24h_vs_2h <- lapply(drugs, function(x){
  get_waterfall_plot(diff_exp_proteins, line1='HS11', line2='HS11', time1='2h', time2='24h', selected_drug=x, title = "WT" )
})
#Uniq proteins in NULL and not in WT
drug_plots_proteins_uniq_to_HS01 <- plot_uniq_proteins_in_null(null_line = 'HS01', wt_line = 'HS11')
multiplot(drug_plots_HS11_24h_vs_2h[[1]], drug_plots_HS01_24h_vs_2h[[1]], drug_plots_proteins_uniq_to_HS01[[1]], cols=3)
multiplot(drug_plots_HS11_24h_vs_2h[[2]], drug_plots_HS01_24h_vs_2h[[2]], drug_plots_proteins_uniq_to_HS01[[2]], cols=3)
multiplot(drug_plots_HS11_24h_vs_2h[[3]], drug_plots_HS01_24h_vs_2h[[3]], drug_plots_proteins_uniq_to_HS01[[3]], cols=3)

#####################
# CASE 2: compare isogenic lines MS03(NULL) & MS12(WT)
#####################
#HS01 / NULL line
drug_plots_MS03_24h_vs_2h <- lapply(drugs, function(x){
  get_waterfall_plot(diff_exp_proteins, line1='MS03', line2='MS03', time1='2h', time2='24h', selected_drug=x, title = "NULL" )
})
#HS11 / WT line
drug_plots_MS12_24h_vs_2h <- lapply(drugs, function(x){
  get_waterfall_plot(diff_exp_proteins, line1='MS12', line2='MS12', time1='2h', time2='24h', selected_drug=x, title = "WT" )
})
#Uniq proteins in NULL and not in WT
drug_plots_proteins_uniq_to_MS03 <- plot_uniq_proteins_in_null(null_line = 'MS03', wt_line = 'MS12')
multiplot(drug_plots_MS12_24h_vs_2h[[1]], drug_plots_MS03_24h_vs_2h[[1]], drug_plots_proteins_uniq_to_MS03[[1]], cols=3)
multiplot(drug_plots_MS12_24h_vs_2h[[2]], drug_plots_MS03_24h_vs_2h[[2]], drug_plots_proteins_uniq_to_MS03[[2]], cols=3)
multiplot(drug_plots_MS12_24h_vs_2h[[3]], drug_plots_MS03_24h_vs_2h[[3]], drug_plots_proteins_uniq_to_MS03[[3]], cols=3)






direction <- factor(sign(flt$medRatio_peptides_cond1) == sign(flt$medRatio_peptides_cond2))
direction <- revalue(direction, c('FALSE' = "opposite", 'TRUE' = 'same'))
flt['direction'] = direction
new_protein_levels <- df %>% group_by(direction) %>% arrange(medRatio_peptides_cond2) %>% .$protein
length(new_protein_levels)




tmp <- melt(flt, id.vars=c('protein', 'drug', 'pval', 'pval_adj', 'peptides_cond1', 'peptides_cond2', 'avgRatio_proteinRep_cond1', 'avgRatio_proteinRep_cond2', 'direction'),
            variable.name = 'condition' , value.name = 'log2ratio')
tmp$condition <- revalue(tmp$condition, c('medRatio_peptides_cond1' = paste0(cellLine1,'_',selected_drug,'_',time1),
                                          'medRatio_peptides_cond2' = paste0(cellLine2,'_',selected_drug,'_',time2)))

p <- ggplot(data=tmp,aes(x=protein, y=log2ratio, fill=condition)) + geom_bar(stat='identity',position="dodge") 
p <- p + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) + ylab('log2ratio') 
p <- p + xlab('kinases') +  theme(legend.position="top")
p <- p + ggtitle(paste('Differentially Expressed Kinases at', time1, 'vs', time2, cellLine1, 'vs', cellLine2 , 'drug: ', selected_drug)) + theme(plot.title=element_text(family="Times", face="bold", size=12))
p + coord_flip()
outplot <- eval(paste0('plots/barplot_diffExp_kinase_',cellLine1, '_vs_', cellLine2 , '_drug_', selected_drug, '_time1_', time1, '_time2_', time2, '.png' ) )
ggsave(filename = outplot, width=12, height=8, units="in")
rm(new_protein_levels)




###############
# case 1 -> comparing kinases of HS01(NF2 -/- shRNA) and HS11(wt) 
###############
# At 24 hours
cellLine1 = 'HS11'    #Wild type
replicate1 = c('run1','run2', 'run3')
replicate2 = c('run1','run2', 'run3')
cellLine2 = 'HS01'   #of interest
time1 = '24h'
time2 = '24h'
drugs = c('Pano', 'CUDC', 'GSK458')
diff_exp_proteins_24h <- ldply(drugs, .fun = function(x){ get_diffExp_proteins(cellLine1, cellLine2, time1, time2, replicate1, replicate2, x) })
drug_plots <- lapply(drugs, function(x){
  get_waterfall_plot(diff_exp_proteins_24h, cellLine1=cellLine1, cellLine2=cellLine2, time1=time1, time2=time2, selected_drug=x )
})
multiplot(plotlist = drug_plots,cols=3)

# At 2 hours 
time1 = '2h'
time2 = '2h'
diff_exp_proteins_2h <- ldply(drugs, .fun = function(x){ get_diffExp_proteins(cellLine1, cellLine2, time1, time2, replicate1, replicate2, x) })
diff_exp_proteins_2h['time'] = '2h'
diff_exp_proteins_24h['time'] = '24h'
diff_exp_proteins <- rbind(diff_exp_proteins_2h, diff_exp_proteins_24h)










   ###Older analysis
             # 
             # results <- list()
             # for(i in 1:nrow(ids)){
             #   protein_id <- ids[i, 'protein_synid']
             #   peptide_id <- ids[i, 'peptide_synid']
             #   peptides <-   get_peptide_data(peptide_id, filter=T)
             #   proteins <-   get_protein_Data(protein_id)
             #   res <- ddply(proteins, .variables =c('Accession', "Gene", "Uniprot",'Unique_Peptides',
             #                                        'Peptides', 'cellLine',  'replicate', 'drug', 
             #                                        'count' ),
             #                .fun = get_statistics)
             #   res['pval_adj'] <- p.adjust(res$pval, method='fdr')
             #   results[[i]] = res
             # }
             # results <- do.call(rbind, results)
             # 
             # ###store the protein p-val results
             # write.table(results, file="Synodos_proteins_pvals_2h_vs_24h.csv",
             #             sep="\t")
             # synStore(File("Synodos_proteins_pvals_2h_vs_24h.csv",
             #               parentId = 'syn4259360'))
             # unlink("Synodos_proteins_pvals_2h_vs_24h.csv")
             # 
             
             
#              get_volcano_plot <- function(data,title=NA){
#                p <- figure(title=title) %>% 
#                  ly_points(log2foldchange, -log10(pval_adj), data=data,
#                            hover=list(Gene,log2foldchange, pval_adj, peptides_cond1,peptides_cond2),
#                            color=drug) %>%
#                  ly_abline(v=0, color='red', width=4)
#                p
#              }
#              
#              plot_significant_proteins <- function(data,title=NA){
#                #get the order of genes sorted by mean fold change across all the drugs
#                #   gene_order <- ddply(data, .variables=c('Gene'),
#                #                 .fun = function(x) mean(x$log2foldchange)) %>% 
#                #                 arrange(desc(V1)) %>% select(Gene)
#                #   ggene_order <- gene_order$Gene
#                p <- ggplot(data=data,aes(x=Gene, y=log2foldchange)) 
#                p <- p + geom_bar(stat='identity',position="dodge") 
#                p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#                p <- p + facet_grid(drug ~ .) 
#                p <- p + xlab('Top diff expressed kinases (FDR .10)') + ggtitle(title)
#                p
#              }
#              
#              
#              