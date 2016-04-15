#install.packages(c("WGCNA", "cluster"))
library(WGCNA) 
#Allowing 6 threads for WGCNA
allowWGCNAThreads(6)
library(cluster) 
library(synapseClient)
options(stringsAsFactors = FALSE) 

synapseLogin()

setwd("/external-data/DAT_123__Synodos_NF2//publicData//expressionData/networkAnalysis/")


#replace the following with the relevant study file id's from synapse
geo_id <- 'GSE43290'
study_synId <- 'syn2554761'
geneExp_synId <- 'syn2563614'
phenotype_synId <- 'syn2554764'


#download the normalized expression study
geneExp <- synGet(geneExp_synId)
geneExp <- read.table(geneExp@filePath, sep="\t", header=T)
phenotype <- synGet(phenotype_synId)
phenotype <- read.table(phenotype@filePath, sep="\t",header=T)


datExp <- t(geneExp)
rownames(datExp) <- gsub('_.*', '', rownames(datExp), perl=T)
nGenes = ncol(datExp)
nSamples = nrow(datExp)


# rownames(datExp)
# colnames(datExp)

gsg = goodSamplesGenes(datExp, verbose = 3);
gsg$allOK


sampleTree = flashClust(dist(datExp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
png(file = paste0(geo_id,"_sampleClustering.png"));
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()




#relevant clinical traits
datTraits <- phenotype[,c('who.grade', 'gender', 'age', 'cytogenetic.group', 'histology','sampleId')]
traitRows <- match(rownames(datExp), phenotype$sampleId)
datTraits <- datTraits[traitRows,]
rownames(datTraits) = datTraits$sampleId
datTraits$sampleId <- NULL


#convert character cols to  numeric
datTraits$who.grade = as.factor(datTraits$who.grade)
levels(datTraits$who.grade) = c(1:length(levels(datTraits$who.grade)))
datTraits$gender = as.factor(datTraits$gender)
levels(datTraits$gender) = c(1:length(levels(datTraits$gender)))
datTraits$cytogenetic.group = as.factor(datTraits$cytogenetic.group)
levels(datTraits$cytogenetic.group) = c(1:length(levels(datTraits$cytogenetic.group)))
datTraits$histology = as.factor(datTraits$histology)
levels(datTraits$histology) = c(1:length(levels(datTraits$histology)))


# Re-cluster samples
sampleTree2 = flashClust(dist(datExp), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraits <- apply(datTraits, 2, as.numeric)
traitColors = numbers2colors(datTraits, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
png(file = paste0(geo_id,"_samplesDendogram_traits.png"));
datTraits <- apply(datTraits, 2, as.numeric)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()




# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
random_cols <- sample(1:ncol(datExp),5000)

sft = pickSoftThreshold(datExp[,random_cols], powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
png(file = paste0(geo_id,"_softPowerCalc.png"), height=8, width=12, units="in", res=300);
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



#single block analysis for all genes
#takes longer time
net_GSE43290 = blockwiseModules(datExp,corType="bicor", 
                       maxBlockSize=21000, networkType="signed", power=14, minModuleSize=30,
                       mergeCutHeight=.25, numericLabels=TRUE, saveTOMs=TRUE, 
                       pamRespectsDendro=FALSE, saveTOMFileBase="GSE16581_network", verbose=5)

save(net_GSE43290, file="GSE43290_networkAnalysis.RData")
load("GSE43290_networkAnalysis.RData")

# Plot the network dendogram and modules
png(file = paste0(geo_id,"_networkDendogram_moduleColors.png"), width=12, height=10, units="in", res=300);
moduleColors <- labels2colors(net_GSE43290$colors)
plotDendroAndColors(net_GSE43290$dendrograms[[1]], moduleColors,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE43290: denovo network modules")
dev.off()                    

table(moduleColors, net_GSE43290$colors)


moduleTraitCor = cor(net_GSE43290$MEs, datTraits, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = colnames(net_GSE43290$MEs),
               ySymbols = colnames(net_GSE43290$MEs),
               colorLabels = T,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#######
#gene relationship to trait and importance modules
#######

cytogenetic.group = as.data.frame(datTraits[,'cytogenetic.group'])
names(cytogenetic.group) = "cytogenetic.group"

#gene module correlation
geneModuleMembership = as.data.frame(cor(datExp, net_GSE43290$MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

#gene trait correlation
geneTraitSignificance <- as.data.frame(cor(datExp, cytogenetic.group, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ncol(datExp)))

module = 4
moduleName = "ME4"
moduleGenes = net_GSE43290$colors == module
column = match(moduleName, colnames(geneModuleMembership) )

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", moduleName, "module"),
                   ylab = "Gene significance for cytogenetic group",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")





sum(moduleGenes)
sum(colnames(datExp) %in% c('NF2'))
paste(colnames(datExp)[moduleGenes], collapse=",")

table(net_GSE43290$colors)

table(phenotype['cytogenetic.group'])


library("org.Hs.eg.db")
colnames(datExp)[moduleGenes]
symbols <- keys(org.Hs.eg.db, keytype="SYMBOL")

x <- select(org.Hs.eg.db, keys=symbols, columns=c("CHR"), keytype="SYMBOL")
y <- x[x$SYMBOL %in% colnames(datExp)[moduleGenes],]
barplot(table(y$CHR), xlab="chr", ylab="#genes")
table(phenotype$who.grade)
