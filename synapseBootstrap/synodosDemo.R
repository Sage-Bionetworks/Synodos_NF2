require(synapseClient)
require(rGithubClient)
require(affy)
require(simpleaffy)
require(limma)
require(ggplot2)

## GET THE GITHUB REPO WHERE THIS CODE IS STORED
repo <- getRepo("Sage-Bionetworks/Synodos_NF2")
thisCode <- getPermlink(repo, "synapseBootstrap/synodosDemo.R")

## PROVIDE THE FOLDER ID WHERE YOU WANT TO STORE YOUR FILES
folderId <- ""

## GRAB EXPRESSION AND CLINICAL DATA FOR GSE43290
exprFile <- synGet("syn2482541")
clinFile <- synGet("syn2485151")

exprDat <- read.table(getFileLocation(exprFile), as.is=T)
colnames(exprDat) <- sapply(strsplit(colnames(exprDat), "_", fixed=T), "[[", 1)
exprDat <- as.matrix(exprDat)
clinDat <- read.delim(getFileLocation(clinFile), as.is=T)
rownames(clinDat) <- clinDat$sampleId
clinDat <- clinDat[ colnames(exprDat), ]

## ONLY FOUR NORMALS - BUT LOOK AT DIFFENCES BETWEEN NORMAL MENINGE AND MENINGIOMAS
myFit <- lmFit(exprDat, model.matrix(~ clinDat$histology != "normal meninge"))
myFit <- eBayes(myFit)

## VOLCANO PLOT
plotDF <- data.frame(log2fc = myFit$coefficients[, 2],
                     mlog10pval = -1*log10(myFit$p.value[, 2]))

myPlot <- ggplot(plotDF, aes(log2fc, mlog10pval)) + 
  geom_point() + 
  xlim(c(-6,6)) +
  xlab('normal expression higher -- log2(FC) -- tumor expression higher') + 
  ylab('-log10(eBayes moderated p-value)')

plotPath <- file.path(tempdir(), "gse43290-tumorNormal-volcanoPlot.png")
png(plotPath)
show(myPlot)
dev.off()

plotFile <- synStore(File(path=plotPath, 
                          parentId=folderId), 
                     used=list(
                       list(entity=exprFile, wasExecuted=F),
                       list(entity=clinFile, wasExecuted=F),
                       list(url=thisCode, name=basename(thisCode), wasExecuted=T)))
