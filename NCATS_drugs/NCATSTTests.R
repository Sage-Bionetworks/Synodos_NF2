library(synapseClient)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
synapseLogin()
#######
#
#
# THIS USES OUTDATED MAPPING APPROACH FOR TARGETS (which skews the results) DO NOT USE THIS ANALYSIS 
#
#
#######
x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

drugdat$NF2.Status<-as.character(drugdat$NF2.Status)

drugdat$NF2.Status[drugdat$NF2.Status=="NF2 null "] <- 0 
drugdat$NF2.Status[drugdat$NF2.Status=="NF2 expressing "] <- 1

x<-synGet("syn8398744")@filePath
target_map<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)

drugdat<-filter(drugdat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(Max.Resp, AUC, Hugo_Gene, NF2.Status) 

drugdat<-dlply(drugdat, .variable = "Hugo_Gene")
drugdat<-drugdat[-1]

foo<-lapply(names(drugdat), function(i){
  j <- drugdat[[i]]
  p <- 1
  t <- t.test(AUC~NF2.Status, data = j)
  p <- t$p.value
  p
})

foo <- ldply(foo)
foo <- as.data.frame(cbind(foo, names(drugdat)))
colnames(foo) <- c("p.value","Gene.Name")
foo$bh <- p.adjust(foo$p.value)

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

drugdat<-filter(drugdat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(Max.Resp, AUC, Hugo_Gene, NF2.Status, Cell.Line, Drug) 


Syn6<-drugdat %>% filter(Cell.Line %in% c("Syn6")) %>%
  dplyr::select(Max.Resp, Drug, Hugo_Gene) 
colnames(Syn6) <- c("Syn6", "Drug", "Gene.Name")
syn1<-drugdat %>% filter(Cell.Line %in% c("Syn1")) %>%
  dplyr::select(Max.Resp, Drug, Hugo_Gene)
colnames(syn1) <- c("Syn1", "Drug", "Gene.Name")
dat<-full_join(Syn6, syn1)
dat$Max.Resp.Ratio<-dat$Syn6/dat$Syn1
dat<-full_join(dat, foo)
dat<-distinct(dat)

dat1 <- dat %>% select(Drug, Gene.Name, Syn1, bh)
dat1 <- cbind(dat1, c(rep("Syn1", nrow(dat1))))
colnames(dat1) <- c("Drug", "Gene.Name", "Max.Resp", "bh", "Cell.Line")

dat2 <- dat %>% select(Drug, Gene.Name, Syn6, bh)
dat2 <- cbind(dat2, c(rep("Syn6", nrow(dat2))))
colnames(dat2) <- c("Drug", "Gene.Name", "Max.Resp", "bh", "Cell.Line")

dat3 <- rbind(dat1,dat2)

ggplot(dat3 %>% filter(Gene.Name=="KAT2A"), aes(x=Cell.Line, y=Max.Resp)) +
  ggbeeswarm::geom_beeswarm()

