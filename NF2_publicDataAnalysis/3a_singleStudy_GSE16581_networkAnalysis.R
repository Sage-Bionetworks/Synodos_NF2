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
geo_id <- 'GSE16581'
study_synId <- 'syn2482470'
geneExp_synId <- 'syn2504985'
phenotype_synId <- 'syn2485136'


#download the normalized expression study
geneExp <- synGet(geneExp_synId)
geneExp <- read.table(geneExp@filePath, sep="\t", header=T)

phenotype <- synGet(phenotype_synId)
phenotype <- read.table(phenotype@filePath, sep="\t",header=T)


datExp <- t(geneExp)
rownames(datExp) <- gsub('.CEL.gz', '', rownames(datExp))


# rownames(datExp)
# colnames(datExp)

gsg = goodSamplesGenes(datExp, verbose = 3);
gsg$allOK


sampleTree = flashClust(dist(datExp), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



# Plot a line to show the cut
abline(h = 115, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 115, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExp = datExp[keepSamples, ]
nGenes = ncol(datExp)
nSamples = nrow(datExp)

#relevant clinical traits
datTraits <- phenotype[,c('who.grade', 'recurrence_frequency','vital.status', 'gender', 'age', 'sampleId')]
traitRows <- match(rownames(datExp), phenotype$sampleId)
datTraits <- datTraits[traitRows,]
rownames(datTraits) = datTraits$sampleId
datTraits$sampleId <- NULL
dim(datTraits)

#convert character cols to  numeric
datTraits$vital.status = as.factor(datTraits$vital.status)
levels(datTraits$vital.status) = c(1:2)
datTraits$gender = as.factor(datTraits$gender)
levels(datTraits$gender) = c(1:2)


# Re-cluster samples
sampleTree2 = flashClust(dist(datExp), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits[,c('who.grade', 'recurrence_frequency')] , signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
random_cols <- sample(1:ncol(datExp),5000)
sft = pickSoftThreshold(datExp[,random_cols], powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
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



#single block analysis for all genes
#takes longer time
net_GSE16581 = blockwiseModules(datExp,corType="bicor", 
                               maxBlockSize=21000, networkType="signed", power=6, minModuleSize=30,
                               mergeCutHeight=mergingThresh, numericLabels=TRUE, saveTOMs=TRUE, 
                               pamRespectsDendro=FALSE, saveTOMFileBase="GSE16581_network", verbose=5)


save(net_GSE16581, file="GSE16581_networkAnalysis.RData")

load("GSE16581_networkAnalysis.RData")


# Plot the network dendogram and modules
png(file = paste0(geo_id,"_networkDendogram_moduleColors.png"), width=12, height=10, units="in", res=300);
moduleColors <- labels2colors(net_GSE16581$colors)
plotDendroAndColors(net_GSE16581$dendrograms[[1]], moduleColors,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "GSE16581: denovo network modules")
dev.off()



###########
## Quantifying module-trait associations
##########
#calculating the correlation between module Eigen genes AND 
moduleTraitCor = cor(net_GSE16581$MEs, datTraits, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExp))


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(net_GSE16581$MEs),
               ySymbols = names(net_GSE16581$MEs),
               colorLabels = TRUE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#######
#gene relationship to trait and importance modules
#######
who.grade = as.data.frame(datTraits$who.grade)
names(who.grade) = "who.grade"

#gene module correlation
geneModuleMembership = as.data.frame(cor(datExp, net_GSE16581$MEs, use="p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), ncol(datExp)))
            
#gene trait correlation
geneTraitSignificance <- as.data.frame(cor(datExp, who.grade, use="p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ncol(datExp)))

moduleGenes = net_GSE16581$colors == 10
column = match("ME10", colnames(geneModuleMembership) )

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", "ME10", "module"),
                   ylab = "Gene significance for WHO GRADE",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "red")


paste0(colnames(datExp)[moduleGenes],collapse=', ')




###########
#Intramodular Analysis : identifying genes with high geneTrait significance score and Module membership score
###########



