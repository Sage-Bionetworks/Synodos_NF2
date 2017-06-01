library(synapseClient)
library(plyr)
library(dplyr)
synapseLogin()

##counts and %mapped from bams using samtools
qc<- read.table(synGet("syn9904189")@filePath, sep = "\t", header = TRUE)

ggplot(data = qc, aes(x=run, y=reads)) +
  geom_bar(stat="identity") +
  labs(x="Run", y="Read Count") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("UNC_reseq_readcount.png")

ggplot(data = qc, aes(x=run, y=perc.map)) +
  geom_bar(stat="identity") +
  labs(x="Run", y="Percent Mapped (%)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("UNC_reseq_percentmapped.png")

ggplot(data = qc, aes(x=reads, y=perc.map)) +
  geom_point() +
  labs(x="Read Count", y="Percent Mapped (%)") +
  xlim(0, 110000000) + ylim(0, 100) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
ggsave("UNC_reseq_readcount_percentmapped_scatter.png")

##normalized reads - Salmon output from UNC
norm.reads<-read.table(synGet("syn9779231")@filePath, sep = "\t", header = TRUE)
norm.reads[,2:27]<- sapply(norm.reads[, 2:27], as.numeric)

clust<-hclust(dist(t(norm.reads[,-1])))
plot(clust)


i <- "HS01_panobinostat_Run3_S4"

for(i in colnames(norm.reads)[-1]){
  ggplot(data = norm.reads, aes(x = i)) +
    geom_histogram(stat = "count")
  }
