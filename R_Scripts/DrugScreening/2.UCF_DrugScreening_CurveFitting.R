library(nplr)
library(plyr)
library(rCharts)
library(ggvis)
library(ggplot2)
library(grid)
library(reshape2)
library("synapseClient")

options(stringsAsFactors = F)
source("~/dev/apRs/expression_heatmap.R")



get_drugResponse_stats2 <- function(conc,viability,...){
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
  results['bottom_asymptote'] = res@pars['bottom']
  results['top_asymptote'] = res@pars['top']
  results['hillSlope'] =  res@pars['scal']
  
  fittedVals <- data.frame(fittedX = getXcurve(res),
                           fittedY = getYcurve(res))
  cbind(results,fittedVals)
}


tmp_iterator <- function(df){
  conc <- as.numeric(df[,'conc'])
  viability <- as.numeric(df[,'viability'])
  stats <- get_drugResponse_stats(conc, viability, useLog=F)
  stats
}

tmp_iterator2 <- function(df){
  conc <- as.numeric(df[,'conc'])
  viability <- as.numeric(df[,'viability'])
  stats <- get_drugResponse_stats2(conc, viability, useLog=F)
  stats
}


#calclate the drugResponse stats
# using viab normalized by meanDMSO
dr_fit_by_meanDMSO <- ddply(.data=ds_norm_by_meanDMSO, .variables = c('drug', 'cellLine', 'plate', 'run'), .fun=tmp_iterator)
dr_fit_by_meanDMSO_flt <- filter(dr_fit_by_meanDMSO, goodNess_of_fit > .50)
write.table(dr_fit_by_meanDMSO_flt, file="UCF_DrugScreen_ICVals.csv", col.names=T, sep="\t", quote=F, row.names=F)
synStore(File("UCF_DrugScreen_ICVals.csv", parentId="syn2753201"))

dr_fit_by_meanDMSO2 <- ddply(.data=ds_norm_by_meanDMSO, .variables = c('drug', 'cellLine', 'plate', 'run'), .fun=tmp_iterator2)
dim(dr_fit_by_meanDMSO)


#goodness of fit
hist(dr_fit_by_meanDMSO$goodNess_of_fit, breaks=50)
p <- ggplot(data=dr_res, aes(x=goodNess_of_fit)) + geom_histogram(binwidth=.005) + theme_bw() 
p <- p + ggtitle("hist: goodness of fit of dose response curve/drug ") + xlab("goodness of fit (p-value)")
p <- p + theme(plot.title=element_text(size=rel(.7)), 
               axis.title=element_text(size=rel(.7)), axis.text=element_text(size=rel(.5)))

drug_levels <- ddply(.data=dr_fit_by_meanDMSO_flt, .variables=c('drug'), .fun=function(x) mean(x$IC50,na.rm=T))
drug_levels <- arrange(drug_levels, desc(V1))
drug_levels <- drug_levels$drug

#IC50 across three cell lines across two runs
p <- ggplot(data=dr_fit_by_meanDMSO_flt, aes(x=factor(drug,levels=drug_levels), y=IC50, group=run)) + geom_line(aes(color=run)) 
p <- p + facet_grid(cellLine ~ .) + geom_point(aes(size=-IC50),color='grey50') + theme_bw()
p + theme(axis.text.x=element_text(angle=90, hjust=1)) + xlab('Drug')
ggsave("plots/IC50_var_across_dosages.png", width=10, height=8, dpi=300)


IC_cols <- grep('IC',colnames(dr_fit_by_meanDMSO_flt))
IC_variation <- melt(dr_fit_by_meanDMSO_flt, id.vars = c('drug', 'run', 'cellLine', 'run'), measure.vars = IC_cols)
ggplot(IC_variation, aes(x=variable, y=drug,)) + geom_tile(aes(fill=value)) + facet_grid( ~ cellLine)

?melt              




##############
## Test code
##############

par(mfrow=c(1,1))
row <- sample(1:nrow(ds),1)
s <- ds[row,]

x1 <- filter(ds, drug == s$drug & run == s$run , cellLine == s$cellLine)
norm_viability <- x1$viability / x1$T0 
res <- nplr(x1$conc, norm_viability, useLog=F)
plot(res)

x2 <- filter(ds_norm_by_meanUntreated, drug == s$drug & run == s$run & cellLine == s$cellLine)
res <- nplr(x2$conc, x2$viability,useLog=F)
plot(res)

x3 <- filter(ds_norm_by_meanDMSO, drug == s$drug & run == s$run & cellLine == s$cellLine)
res <- nplr(x3$conc, x3$viability,useLog=F)
plot(res, showEstim=T, showInfl=T)





tmp_iterator2(x3)

rnorm(1000,0,.5)

par(mfrow=c(1,1))
object <- res
x <- getX(object)
y <- getY(object)
newx <- getXcurve(object)
newy <- getYcurve(object)
gof <- round(getGoodness(object), 3)
plot(x, y, col = "aquamarine1", pch = 19, las = 1, cex.axis = 1.25, 
     cex.lab = 1.5)
points(x, y, pch = 1)
lines(newy ~ newx, col = "red3", lwd = 4)
gof <- round(getGoodness(object), 3)
stdErr <- getStdErr(object)
nplr::getStdErr

df <- data.frame(x=newx, y=newy)
ggplot(data=dr_fit_by_meanDMSO2, aes(x=fittedX, y=fittedY, group=cellLine)) + geom_line(aes(color=cellLine)) + facet_grid(run ~ drug)

colnames(dr_fit_by_meanDMSO2)
