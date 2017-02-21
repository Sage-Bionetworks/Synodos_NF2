library("synapseClient")
synapseLogin()

library("nplr")
library("plyr")
library("dplyr")
library("data.table")
library("plotly")
library("reshape2")
library("ggplot2")

get_drugResponse_stats <- function(conc,viability,...){
  res <- nplr(conc, viability,...)
  results <- getAUC(res)
  results['goodNess_of_fit'] <- getGoodness(res)
  results['stdErr'] <- getStdErr(res)
  ICx_est = getEstimates(res, targets= c(.10,.20,.30,.40,.50,.60,.70,.80,.90))
  results['IC10'] = ICx_est[1,'x']
  results['IC20'] = ICx_est[2,'x']
  results['IC30'] = ICx_est[3,'x']
  results['IC40'] = ICx_est[4,'x']
  results['IC50'] = ICx_est[5,'x']
  results['IC60'] = ICx_est[6,'x']
  results['IC70'] = ICx_est[7,'x']
  results['IC80'] = ICx_est[8,'x']
  results['IC90'] = ICx_est[9,'x']
  results['maxEfficacy'] = max(getYcurve(res)) #get the maximum efficacy of the drug
  results['bottom_asymptote'] = res@pars['bottom']
  results['top_asymptote'] = res@pars['top']
  results['hillSlope'] =  res@pars['scal']
  
  fittedVals <- data.frame(fittedX = getXcurve(res),
                           fittedY = getYcurve(res))
  results <- cbind(results,fittedVals)
  results
}


tmp_iterator <- function(df){
  tryCatch({
    stats <- get_drugResponse_stats(df$conc, df$normViability, useLog=T)  
  },error=function(e){
    print(dim(df))
    print(df$conc)
    print(df$normViability)
    print(unique(df$cellLine))
    print(unique(df$drug))
    print(unique(df$experiment))
    stop('stopped')
  })
}

UCF_normViab <- 'syn2773870'
UCF_normViab <- synGet(UCF_normViab)
UCF_normViab <- read.delim(UCF_normViab@filePath, check.names=F, sep="\t", header=T,stringsAsFactors = FALSE)
UCF_normViab$cellLine <- gsub("^ ", "", UCF_normViab$cellLine)
UCF_normViab$cellLine <- gsub("Nf2 --", "Nf2--", UCF_normViab$cellLine)

#drop unnecassary cols
drop_cols <- c('plate', 'medianDMSO', 'viability')
#MGH_normViab <- MGH_normViab[, !colnames(MGH_normViab) %in% drop_cols]
UCF_normViab <- UCF_normViab[, !colnames(UCF_normViab) %in% drop_cols]

# Drug response
doseResp <- ddply(.data=UCF_normViab, .variables = c('drug', 'cellLine', 'experiment'), 
                  .fun = tmp_iterator, .parallel = T)



doseResp_r <- ddply(.data=UCF_normViab, .variables = c('drug', 'cellLine', 'experiment','replicate'), 
                  .fun = tmp_iterator, .parallel = T)

doseResp_new <- ddply(.data=UCF_normViab, .variables = c('drug', 'cellLine'), 
                    .fun = tmp_iterator, .parallel = T)

drugs <- c("CUDC907","GSK2126458","Panobinostat")
###########

# Human
HS_normViab <- UCF_normViab[UCF_normViab$cellLine %in% c("HS01","HS11") & UCF_normViab$drug %in% drugs,]

HS_doseResp <- doseResp[doseResp$cellLine %in% c("HS01","HS11") & doseResp$drug %in% drugs,]
HS_doseResp <- HS_doseResp[order(HS_doseResp$IC50),]

labelVal <- as.integer(c(min(log10(HS_normViab$conc*(1e+6))), median(log10(HS_normViab$conc*(1e+6))), max(log10(HS_normViab$conc*(1e+6)))))

p1 <- ggplot(HS_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p1 <- p1 + geom_point(aes_string(color="cellLine")) 
p1 <- p1 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p1 <- p1 + geom_line(data = HS_doseResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p1 <- p1 + facet_grid(. ~ drug) + theme_bw(base_size = 15)
p1 <- p1 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p1 <- p1 + xlab('log10 micromolar conc') + ylab('cell viability %')
p1 <- p1 + scale_x_continuous(breaks = labelVal, labels = sapply(labelVal, function(x) format(10^x,scientific = T)))
p1
ggsave(filename = "HS_drug_response_2Dplot.pdf",width = 8, height = 4)


HS_doseResp_r <- doseResp_r[doseResp_r$cellLine %in% c("HS01","HS11") & doseResp_r$drug %in% drugs,]
HS_doseResp_r <- HS_doseResp_r[order(HS_doseResp_r$IC50),]
HS_doseResp_r$grp <- paste(HS_doseResp_r$cellLine,HS_doseResp_r$replicate)

p2 <- ggplot(HS_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p2 <- p2 + geom_point(aes_string(color="cellLine")) 
p2 <- p2 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p2 <- p2 + geom_line(data = HS_doseResp_r, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = grp))
p2 <- p2 + facet_grid(. ~ drug) + theme_bw(base_size = 15)
p2 <- p2 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p2 <- p2 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p2

ggsave(filename = "HS_rep_drug_response_2Dplot.pdf",width = 8, height = 4)

# Mouse
MS_normViab <- UCF_normViab[grep('^MS', UCF_normViab$cellLine),]
MS_normViab <- MS_normViab[MS_normViab$drug %in% drugs,]

MS_doseResp <- doseResp[grep('^MS', doseResp$cellLine),]
MS_doseResp <- MS_doseResp[MS_doseResp$drug %in% drugs,]
MS_doseResp <- MS_doseResp[order(MS_doseResp$IC50),]

p3 <- ggplot(MS_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p3 <- p3 + geom_point(aes_string(color="cellLine")) 
p3 <- p3 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p3 <- p3 + geom_line(data = MS_doseResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p3 <- p3 + facet_grid(experiment ~ drug) + theme_bw(base_size = 15)
p3 <- p3 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p3 <- p3 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p3
ggsave(filename = "MS_drug_response_2Dplot.pdf",width = 8, height = 6)

MS_doseResp_r <- doseResp_r[grep('^MS', doseResp_r$cellLine),]
MS_doseResp_r <- MS_doseResp_r[MS_doseResp_r$drug %in% drugs,]
MS_doseResp_r <- MS_doseResp_r[order(MS_doseResp_r$IC50),]
MS_doseResp_r$grp <- paste(MS_doseResp_r$cellLine,MS_doseResp_r$replicate)

p4 <- ggplot(MS_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p4 <- p4 + geom_point(aes_string(color="cellLine")) 
p4 <- p4 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p4 <- p4 + geom_line(data = MS_doseResp_r, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = grp))
p4 <- p4 + facet_grid(experiment ~ drug) + theme_bw(base_size = 15)
p4 <- p4 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p4 <- p4 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p4
ggsave(filename = "MS_rep_drug_response_2Dplot.pdf",width = 8, height = 6)

########################
# Cell line Comparision
########################

# set 1: HS01, HS11
comp_1_normViab <- UCF_normViab[UCF_normViab$cellLine %in% c("HS01","HS11"),]
comp_1_drugResp <- doseResp[doseResp$cellLine %in% c("HS01","HS11"),]

p5 <- ggplot(comp_1_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p5 <- p5 + geom_point(aes_string(color="cellLine"),size=0.5) 
p5 <- p5 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p5 <- p5 + geom_line(data = comp_1_drugResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p5 <- p5 + facet_wrap(~ drug, ncol = 4) + theme_bw(base_size = 14) + theme(legend.position=c(0.85,0.08))
p5 <- p5 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p5 <- p5 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p5
ggsave(filename = "compare1_drug_response_2Dplot.pdf",width = 6, height = 8)

# set 2: MS01, MS02, MS03
comp_2_normViab <- UCF_normViab[UCF_normViab$cellLine %in% c("MS01","MS02","MS03"),]
comp_2_drugResp <- doseResp_new[doseResp_new$cellLine %in% c("MS01","MS02","MS03"),]

p6 <- ggplot(comp_2_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p6 <- p6 + geom_point(aes_string(color="cellLine"),size=0.5) 
p6 <- p6 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p6 <- p6 + geom_line(data = comp_2_drugResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p6 <- p6 + facet_wrap( ~ drug, ncol = 4) + theme_bw(base_size = 14) + theme(legend.position=c(0.85,0.08))
p6 <- p6 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p6 <- p6 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p6
ggsave(filename = "compare2_drug_response_2Dplot.pdf",width = 6, height = 8)


# set 3: MS01, MS02, HS01
comp_3_normViab <- UCF_normViab[UCF_normViab$cellLine %in% c("MS01","MS02","HS01"),]
comp_3_drugResp <- doseResp_new[doseResp_new$cellLine %in% c("MS01","MS02","HS01"),]

p7 <- ggplot(comp_3_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p7 <- p7 + geom_point(aes_string(color="cellLine"),size=0.5) 
p7 <- p7 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p7 <- p7 + geom_line(data = comp_3_drugResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p7 <- p7 + facet_wrap( ~ drug, ncol = 4) + theme_bw(base_size = 14) + theme(legend.position=c(0.85,0.08))
p7 <- p7 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p7 <- p7 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p7
ggsave(filename = "compare3_drug_response_2Dplot.pdf",width = 6, height = 8)

# set 4: HS11, MS12
comp_4_normViab <- UCF_normViab[UCF_normViab$cellLine %in% c("HS11","MS12"),]
comp_4_drugResp <- doseResp[doseResp$cellLine %in% c("HS11","MS12"),]

p8 <- ggplot(comp_4_normViab, aes(x = log10(conc*(1e+6)), y = normViability*100)) 
p8 <- p8 + geom_point(aes_string(color="cellLine"),size=0.5) 
p8 <- p8 + scale_color_brewer(type = "qual", palette = 2, direction = 1)
p8 <- p8 + geom_line(data = comp_4_drugResp, aes(x = fittedX+6, y = fittedY*100, colour = cellLine, group = cellLine))
p8 <- p8 + facet_wrap( ~ drug, ncol = 4) + theme_bw(base_size = 14) + theme(legend.position=c(0.85,0.08))
p8 <- p8 + geom_hline(aes(yintercept=50), color='grey50', linetype='dashed')
p8 <- p8 + xlab('log10 micromolar conc') + ylab('cell viability %') 
p8
ggsave(filename = "compare4_drug_response_2Dplot.pdf",width = 6, height = 8)


# ###########
# # 3D plots
# ###########
# 
# # Human
# p1 <- plot_ly(HS_doseResp_orderIC50, x = fittedX, y = fittedY*100, z = drug, 
#              type = "scatter3d", mode = "lines", color = drug, 
#              #symbol = cellLine, symbols = c("cross", "square"),
#              hoverinfo = "text",marker = list(size=5),
#              text = paste("drug:", HS_doseResp_orderIC50$drug, 
#                           "</br>cell line:", HS_doseResp_orderIC50$cellLine)) %>% 
#   layout(#title = "Demo Plot",
#          scene = list(
#            xaxis = list(title = "micromolar conc"), 
#            yaxis = list(title = "cell viability %"), 
#            zaxis = list(title = "drug")))
# p1
# 
# # Mouse
# p2 <- plot_ly(MS_doseResp_orderIC50, x = fittedX, y = fittedY*100, z = drug, 
#               type = "scatter3d", mode = "markers+lines", color = drug, 
#               #symbol = cellLine, symbols = c("cross", "square"),
#               hoverinfo = "text",marker = list(size=2),
#               text = paste("drug:", MS_doseResp_orderIC50$drug, 
#                            "</br>cell line:", MS_doseResp_orderIC50$cellLine)) %>% 
#   layout(#title = "Demo Plot",
#     scene = list(
#       xaxis = list(title = "micromolar conc"), 
#       yaxis = list(title = "cell viability %"), 
#       zaxis = list(title = "drug")))
# p2
