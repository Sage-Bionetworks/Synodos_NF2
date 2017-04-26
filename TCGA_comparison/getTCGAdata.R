source("cBioPortalData.R")
library(synapseClient)
exp <- getDisExpressionData(study = "tcga")
write.table("expression.txt", sep = "t")
synStore(File("expression.txt", parentID= "syn9705598"), executed = c(this.file, 