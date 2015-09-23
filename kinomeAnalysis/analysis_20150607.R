library(synapseClient)
library(dplyr)
synapseLogin()

k <- synGet('syn4259365')
k <- read.table(k@filePath, header=T, sep="\t", check.names=F)

data <- k %>% filter(cellLine %in% c('Syn5', 'Syn1') & time  == '2h')

library("devtools")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")


#withtin each cellLine overlap
syn1 = data %>% filter(cellLine == 'Syn1')
syn1_cudc = (syn1 %>% filter(drug == 'CUDC') %>% select(Gene))$Gene
syn1_pano = (syn1 %>% filter(drug == 'Pano') %>% select(Gene))$Gene
syn1_gsk458 =  (syn1 %>% filter(drug == 'GSK458') %>% select(Gene))$Gene
library(gplots)
venn(list(syn1_gsk458=syn1_gsk458,
          syn1_pano=syn1_pano,syn1_cudc=syn1_cudc))

syn5 = data %>% filter(cellLine == 'Syn5')
syn5_cudc = (syn5 %>% filter(drug == 'CUDC') %>% select(Gene))$Gene
syn5_pano = (syn5 %>% filter(drug == 'Pano') %>% select(Gene))$Gene
syn5_gsk458 =  (syn5 %>% filter(drug == 'GSK458') %>% select(Gene))$Gene

#across two cell lines comparison
venn(list(syn5_gsk458=syn5_gsk458,syn1_gsk458=syn1_gsk458))


# counts 
syn1_cudc = (syn1 %>% filter(drug == 'CUDC') %>% select(count))$count
syn1_pano = (syn1 %>% filter(drug == 'Pano') %>% select(count))$count
png('plots/1.png')
plot(syn1_cudc, syn1_pano)
dev.off()




