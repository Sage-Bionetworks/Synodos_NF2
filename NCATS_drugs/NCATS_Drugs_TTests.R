library(tidyverse)
library(synapseClient)
synapseLogin()

ncats <-read.table(synGet("syn8314523")@filePath, header = T, sep = "\t") %>% filter(Cell.line %in% c("HS01", "HS12", "Syn5", "Syn1"))

res<-lapply(unique(ncats$Sample.Name), function(x){
  drugs.ncats <- ncats %>% filter(Sample.Name == x)
  null <- drugs.ncats %>% filter(NF2.status != "NF2 expressing")
  wt <- drugs.ncats %>% filter(NF2.status != "NF2 null")
  bar<-t.test(wt$AUC, null$AUC)
  c(bar$statistic, "pval" = bar$p.value)
})

names(res) <- unique(ncats$Sample.Name)
res <- ldply(res)

##nothing significant t.test across AUC ....

res<-lapply(unique(ncats$Sample.Name), function(x){
  drugs.ncats <- ncats %>% filter(Sample.Name == x)
  null <- drugs.ncats %>% filter(NF2.status != "NF2 expressing")
  wt <- drugs.ncats %>% filter(NF2.status != "NF2 null")
  bar<-t.test(wt$Max.Resp, null$Max.Resp)
  c(bar$statistic, "pval" = bar$p.value)
})

names(res) <- unique(ncats$Sample.Name)
res <- ldply(res)
##nothing significant t.test across Max Resp ....

