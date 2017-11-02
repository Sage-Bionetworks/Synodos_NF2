library(synapseClient)
library(VennDiagram)
library(ggplot2)
library(ggrepel)
synapseLogin()

sch.new <- read.table(synGet("syn9884855")@filePath, sep = "\t", header = T) %>% filter(BH<0.05)

sch.old <- read.table(synGet("syn9884855", version = 11)@filePath, sep = "\t", header = T) %>% filter(BH<0.05)


for(x in unique(sch.new$comparison)){
  print(x)
  sch.new.foo <- filter(sch.new, comparison == x)
  sch.old.foo <- filter(sch.old, comparison == x)
  ens<-list(na.omit(unique(sch.new.foo$ensembl)), na.omit(unique(sch.old.foo$ensembl)))
  names(ens) <- c("new", "old")
  venn.diagram(ens, filename = paste0(x,"_ensembl_venn.png"),
               imagetype = "png",
               compression = "lzw",
               height = 1200,
               width = 1200,
               resolution = 300,
               units = "px")
}

sch.new2 <- read.table(synGet("syn9884855")@filePath, sep = "\t", header = T) 

for(x in unique(sch.new2$comparison)){
  sch.new.foo <- filter(sch.new2, comparison == x)
  ggplot(data = sch.new.foo, aes(x = logFC, y = -log(BH))) +
    ggthemes::theme_few() +
    geom_point(aes(color = BH < 0.05)) +
    scale_color_manual(values=c("FALSE"="lightgrey","TRUE"="#586BA4")) +
    geom_label_repel(data = filter(sch.new.foo, BH < 0.05) %>% 
                       top_n(10, logFC),
                     aes(x = logFC, y = -log(BH), label = Hugo_Gene), fill = "#FF7780") +
    geom_label_repel(data = filter(sch.new.foo, BH < 0.05) %>% 
                       top_n(10, -logFC),
                     aes(x = logFC, y = -log(BH), label = Hugo_Gene), fill = "#60BAFF") +
    ggtitle(x)
   ggsave(paste0(x,"_VolcanoPlotsforUCF.png"))
}
