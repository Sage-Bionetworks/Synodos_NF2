library(synapseClient)
library(tidyverse)
synapseLogin()


dat <- read.table(synGet("syn8314523")@filePath, sep = "\t", header = T)

mat <- dat %>%
    select(Cell.line, Sample.Name, AUC) %>%
    filter(Sample.Name != "", Cell.line != "MS02") %>%
    group_by(Sample.Name, Cell.line) %>%
    summarize("meanAUC" = mean(AUC)) %>%
    ungroup() %>%
    spread(Sample.Name, meanAUC) %>%
    remove_rownames() %>%
    column_to_rownames("Cell.line") %>%
    t() %>%
    as.data.frame()

mat <- mat[complete.cases(mat),]

map <- select(dat, Sample.Name, Gene.Symbol)

pca <- prcomp(mat)
pca <- prcomp(mat)$x %>%
    as.data.frame() %>%
    rownames_to_column("Sample.Name") %>%
    left_join(map)

ggplot(data = pca, aes(x=PC1, y=PC2)) +
    geom_point(aes(color = Gene.Symbol)) +
    theme(legend.position="none")

hist(scale(dat$AUC, center = T))

library(limma)
library(edgeR)
dge=DGEList(counts = mat)

design <- as.data.frame(c(1:5))
design$Cell.Line <- colnames(mat)
design$genotype <- c("MUT", "WT", "WT", "MUT","MUT")
design$species <- c("Hu", "Hu", "Hu", "Hu", "Hu")
design$tumorType <- c("Sch", "Sch", "Men", "Men", "Men")

design2 <- model.matrix(~ design$genotype)
v <- voom(mat, design2, plot = T)

fit <- lmFit(v, design2)
fit <- eBayes(fit)
tab<-topTable(fit, coef=ncol(design2), number = 10000)

ggplot(data = dat %>% filter(Sample.Name=="IMD-0354", Cell.line != "MS02"), aes(x=Cell.line, y=AUC.Fit))+
    geom_point(aes(color = NF2.status))


