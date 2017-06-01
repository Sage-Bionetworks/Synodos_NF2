library(synapseClient)
library(dplyr)
library(ggplot2)
library(viridis)
library(plotrix)
library(openxlsx)
synapseLogin()
dat<-read.table(synGet("syn5840701")@filePath, sep = "\t", header = TRUE, comment.char = "")

##syn1-syn5
syn1syn5 <- dat %>% filter(cellLine=="Syn5", referenceSample=="Syn1") %>% 
  group_by(Gene) %>% summarize(mean = mean(log2ratio), se= std.error(log2ratio))  %>% 
  ungroup() %>% filter(se<=abs(mean))
  
syn1syn5.top <- distinct(rbind((syn1syn5 %>% top_n(10, mean)), (syn1syn5 %>% top_n(10, -mean))))

ggplot(syn1syn5.top, aes(x=Gene %>% reorder(mean), y=mean, fill = mean)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "C") +
  labs(x = "Kinase", y = paste("mean fold change (log2)"), title = "Syn5 vs Syn1") +
  coord_flip()  +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 18))

ggsave("syn1syn5baselineKinome.pdf", height = 8, width = 5)

##hs11-hs01
hs11hs01 <- dat %>% filter(cellLine=="HS01", referenceSample=="HS11") %>% 
  group_by(Gene) %>% summarize(mean = mean(log2ratio), se= std.error(log2ratio)) %>% 
  ungroup() %>% filter(se<=abs(mean))

hs11hs01.top <- distinct(rbind((hs11hs01 %>% top_n(10, mean)), (hs11hs01 %>% top_n(10, -mean))))

ggplot(hs11hs01.top, aes(x=Gene %>% reorder(mean), y=mean, fill = mean)) +
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "C") +
  labs(x = "Kinase", y = paste("mean fold change (log2)"), title = "HS01 vs HS11") +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 18))

ggsave("hs01hs11baselineKinome.pdf", height = 8, width = 5)

## MSO3 vs MS12
mouse<-read.table(synGet("syn5840664")@filePath, sep = "\t", header = TRUE, comment.char = "") %>% 
  select(Gene, MS03, MS12) %>% 
  mutate(log2ratio = log2(MS03/MS12)) %>% 
  filter(log2ratio != "Inf" & log2ratio != "-Inf" & log2ratio != "NA" & log2ratio != "NaN")

mouse.top <- distinct(rbind((mouse %>% top_n(10, log2ratio)), (mouse %>% top_n(10, -log2ratio))))

ggplot(mouse.top, aes(x=Gene %>% reorder(log2ratio), y=log2ratio, fill = log2ratio)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "C") +
  labs(x = "Kinase", y = paste("mean fold change (log2)"), title = "MS03 vs MS12") +
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "#FFFFFF"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position="none",
        plot.margin=unit(c(1,1,1,1), "cm"), text=element_text(size = 18))
        

ggsave("ms0baselineKinome.pdf", height = 8, width = 5)