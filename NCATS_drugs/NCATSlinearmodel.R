library(synapseClient)
library(dplyr)
library(tidyr)
synapseLogin()

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

x<-synGet("syn8398744")@filePath
target_map<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)

drugdat<-filter(drugdat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% select(AUC, Hugo_Gene) %>% 

lm <- lm(drugdat$Hugo_Gene, drugdat$AUC)

View(lm)
