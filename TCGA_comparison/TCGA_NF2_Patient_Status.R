library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(ggplot2)

this.file = "x"

mutation.files<-list.files("mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("mutations_by_cancer/",x, sep=""), header = TRUE))

cancers<-sub("mutations.txt", "", mutation.files)
names(mutations) <- cancers

nf2.stat<-lapply(1:35, function(i){
  ##get samples with nf2 mutation
  nf2<-filter(mutations[[i]], gene=="NF2")
  nf2<-nf2[,-1]
  nf2<-as.data.frame(t(nf2))
  nf2$samples<-rownames(nf2)
  try(colnames(nf2) <- c("mutation", "samples"))
  cancer.type<-sub("mutations.txt", "", mutation.files[[i]])
  print(cancer.type)
  return(nf2)
})

names(nf2.stat) <- cancers
nf2.stat.df <- ldply(nf2.stat, .id = "cancer")
nf2.stat.df <- filter(nf2.stat.df, !is.na(mutation))

exp <- getDisExpressionData(study = "tcga")
exp <- select(exp, one_of(nf2.stat.df$samples))
nf2.stat.df <- filter(nf2.stat.df, samples %in% colnames(exp))
exp$gene_id <- rownames(exp)

##get expression data
HS01.1 <- read.table(synGet("syn5567197")@filePath, header = TRUE, sep = "\t") %>% select(gene_id, TPM)
colnames(HS01.1) <- c("gene_id", "HS01.1")
HS01.1 <- inner_join(HS01.1, exp)

HS01.1.names <- colnames(HS01.1)[-1]

x<-as.data.frame(sapply(HS01.1.names, function(i){
  x<-cor(HS01.1$HS01.1, HS01.1[i])
}))

colnames(x) <- "HS01.1"
x$samples <- rownames(x)

x<-inner_join(x, nf2.stat.df)
x$mutation <- sub(1, "mut", x$mutation)
x$mutation <- sub(0, "wt", x$mutation)

#####HS01.2
HS01.2 <- read.table(synGet("syn5567201")@filePath, header = TRUE, sep = "\t") %>% select(gene_id, TPM)
colnames(HS01.2) <- c("gene_id", "HS01.2")
HS01.2 <- inner_join(HS01.2, exp)

HS01.2.names <- colnames(HS01.2)[-1]

temp<-as.data.frame(sapply(HS01.2.names, function(i){
  t<-cor(HS01.2$HS01.2, HS01.2[i])
}))

colnames(temp) <- "HS01.2"
temp$samples <- rownames(temp)
x<-inner_join(x, temp)

#####Syn5
Syn5.1 <- read.table(synGet("syn5269658")@filePath, header = TRUE, sep = "\t") %>% select(gene_id, TPM)
colnames(Syn5.1) <- c("gene_id", "Syn5.1")
Syn5.1 <- inner_join(Syn5.1, exp)

Syn5.1.names <- colnames(Syn5.1)[-1]

temp<-as.data.frame(sapply(Syn5.1.names, function(i){
  t<-cor(Syn5.1$Syn5.1, Syn5.1[i])
}))

colnames(temp) <- "Syn5.1"
temp$samples <- rownames(temp)
x<-inner_join(x, temp)

#####Syn6
Syn6 <- read.table(synGet("syn5269682")@filePath, header = TRUE, sep = "\t") %>% select(gene_id, TPM)
colnames(Syn6) <- c("gene_id", "Syn6")
Syn6 <- inner_join(Syn6, exp)

Syn6.names <- colnames(Syn6)[-1]

temp<-as.data.frame(sapply(Syn6.names, function(i){
  t<-cor(Syn6$Syn6, Syn6[i])
}))

colnames(temp) <- "Syn6"
temp$samples <- rownames(temp)
x<-inner_join(x, temp)


ggplot(x, aes(mutation, HS01.1, fill = mutation)) +
  geom_boxplot() 

ggplot(x, aes(mutation, HS01.2, fill = mutation)) +
  geom_boxplot()

ggplot(x, aes(mutation, Syn5.1, fill = mutation)) +
  geom_boxplot()

ggplot(x, aes(mutation, Syn6, fill = mutation)) +
  geom_boxplot()

ggplot(x, aes(cancer, HS01.1, fill = cancer)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(x, aes(cancer, HS01.2, fill = cancer)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(x, aes(cancer, Syn5.1, fill = cancer)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(x, aes(cancer, Syn6, fill = cancer)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

