source("cBioPortalData.R")
library(synapseClient)
library(data.table)
synapseLogin()

ex2 <- "https://raw.githubusercontent.com/Sage-Bionetworks/Synodos_NF2/master/TCGA_comparison/getTCGAdata.R"

##TCGA function source from Sara Gosline
exp <- getDisExpressionData(study = "tcga")
fwrite(exp, "TCGAexpressionMatrix.txt", sep = "\t")
synStore(File("TCGAexpressionMatrix.txt", parentId= "syn9705598"), executed = c(ex2)) 
