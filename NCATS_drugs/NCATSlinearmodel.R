library(synapseClient)
library(dplyr)
library(tidyr)
library(speedglm)
library(ggplot2)
library(viridis)
library(purrr)
library(FField)
synapseLogin()

##train and test split function, Steven Turner, UVA
splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/2))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
}

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

x<-synGet("syn8398744")@filePath
target_map<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)

drugdat<-filter(drugdat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(AUC, Hugo_Gene) 

lm <- speedlm(drugdat$AUC~drugdat$Hugo_Gene)
data <- summary(lm)$coefficients

## modeling difference in max response
x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

syn5<-drugdat %>% filter(Cell.Line %in% c("Syn5")) %>%
  dplyr::select(Max.Resp, Drug) 
colnames(syn5) <- c("Syn5", "Drug")
syn1<-drugdat %>% filter(Cell.Line %in% c("Syn1")) %>%
  dplyr::select(Max.Resp, Drug)
colnames(syn1) <- c("Syn1", "Drug")
dat<-full_join(syn5, syn1)
dat$Max.Resp.Ratio<-dat$Syn5/dat$Syn1 
dat<-dat %>%  dplyr::select(Drug, Max.Resp.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(Max.Resp.Ratio, Hugo_Gene, Drug) 

data.list<-splitdf(dat)

lm.gene <- speedlm(Max.Resp.Ratio~Hugo_Gene, data = data.list$trainset)
data<-as.data.frame(summary(lm.gene)$coefficients)

predict.speedlm(lm.gene, data.list$testset)

lm.drug <- speedlm(dat$Max.Resp.Ratio~dat$Drug)

## syn6/syn1

Syn6<-drugdat %>% filter(Cell.Line %in% c("Syn6")) %>%
  dplyr::select(Max.Resp, Drug) 
colnames(Syn6) <- c("Syn6", "Drug")
syn1<-drugdat %>% filter(Cell.Line %in% c("Syn1")) %>%
  dplyr::select(Max.Resp, Drug)
colnames(syn1) <- c("Syn1", "Drug")
dat<-full_join(Syn6, syn1)
dat$Max.Resp.Ratio<-dat$Syn1/dat$Syn6 
dat<-dat %>%  dplyr::select(Drug, Max.Resp.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(Max.Resp.Ratio, Hugo_Gene, Drug) 

lm.gene <- speedlm(Max.Resp.Ratio~Hugo_Gene, data = dat)
data<-as.data.frame(summary(lm.gene)$coefficients)

lm.drug <- speedlm(dat$Max.Resp.Ratio~dat$Drug)

## hs01vshs11
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

hs01<-drugdat %>% filter(Cell.Line %in% c("HS01")) %>%
  dplyr::select(Max.Resp, Drug) 
colnames(hs01) <- c("HS01", "Drug")
hs12<-drugdat %>% filter(Cell.Line %in% c("HS12")) %>%
  dplyr::select(Max.Resp, Drug)
colnames(hs12) <- c("HS12", "Drug")
dat<-full_join(hs01, hs12)
dat$Max.Resp.Ratio<-dat$HS01/dat$HS12 
dat<-dat %>%  dplyr::select(Drug, Max.Resp.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(Max.Resp.Ratio, Hugo_Gene, Drug) 

lm.gene <- speedlm(Max.Resp.Ratio~Hugo_Gene, data = dat)
data<-as.data.frame(summary(lm.gene)$coefficients)

lm.drug <- speedlm(dat$Max.Resp.Ratio~dat$Drug)
data<-as.data.frame(summary(lm.drug)$coefficients)


################################################################################ 
data$gene <- sub(".*Gene", "", rownames(data))

ggplot(data, aes(x= coef, y=-p.value)) +
  geom_point(stat = "identity", aes(color=(p.value))) + 
  scale_color_viridis() +
  geom_text(data=subset(data, p.value<0.05), aes(label = gene), size = 2)

################################################################################ 


## modeling difference in AUC
syn5<-drugdat %>% filter(Cell.Line %in% c("Syn5")) %>%
  dplyr::select(AUC, Drug) 
colnames(syn5) <- c("Syn5", "Drug")
syn1<-drugdat %>% filter(Cell.Line %in% c("Syn1")) %>%
  dplyr::select(AUC, Drug)
colnames(syn1) <- c("Syn1", "Drug")
dat<-full_join(syn5, syn1)
dat$AUC.Ratio<-dat$Syn5/dat$Syn1 
dat<-dat %>%  dplyr::select(Drug, AUC.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(AUC.Ratio, Hugo_Gene, Drug) 

lm.gene <- speedlm(AUC.Ratio~Hugo_Gene, data = dat)
data<-as.data.frame(coef(lm.gene))

lm.drug <- speedlm(dat$AUC.Ratio~dat$Drug)

##syn6 vs syn1
syn6<-drugdat %>% filter(Cell.Line %in% c("Syn5")) %>%
  dplyr::select(AUC, Drug) 
colnames(syn6) <- c("Syn6", "Drug")
syn1<-drugdat %>% filter(Cell.Line %in% c("Syn1")) %>%
  dplyr::select(AUC, Drug)
colnames(syn1) <- c("Syn1", "Drug")
dat<-full_join(syn6, syn1)
dat$AUC.Ratio<-dat$Syn1/dat$Syn6 
dat<-dat %>%  dplyr::select(Drug, AUC.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(AUC.Ratio, Hugo_Gene, Drug) 

lm.gene <- speedlm(AUC.Ratio~Hugo_Gene, data = dat)
data<-as.data.frame(summary(lm.gene)$coefficients)
lm.drug <- speedlm(dat$AUC.Ratio~dat$Drug)
data<-as.data.frame(summary(lm.drug)$coefficients)

## hs01vshs11
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)
colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Drug', 'AUC', 'AUC.Fit', 'Gene.Name') 

hs01<-drugdat %>% filter(Cell.Line %in% c("HS01")) %>%
  dplyr::select(AUC, Drug) 
colnames(hs01) <- c("HS01", "Drug")
hs12<-drugdat %>% filter(Cell.Line %in% c("HS12")) %>%
  dplyr::select(AUC, Drug)
colnames(hs12) <- c("HS12", "Drug")
dat<-full_join(hs01, hs12)
dat$AUC.Ratio<-dat$HS01/dat$HS12 
dat<-dat %>%  dplyr::select(Drug, AUC.Ratio)

dat<-filter(dat, Drug %in% target_map$Drug) %>% 
  full_join(target_map) %>% dplyr::select(AUC.Ratio, Hugo_Gene, Drug) 

lm.gene <- speedlm(AUC.Ratio~Hugo_Gene, data = dat)
data<-as.data.frame(summary(lm.gene)$coefficients)

lm.drug <- speedlm(dat$AUC.Ratio~dat$Drug)
data<-as.data.frame(summary(lm.drug)$coefficients)
