library(synapseClient)
synapseLogin()

#get processed drug screen data
processed_drugScreen_data_synId <- "syn6138237"
data <- read.csv(synGet(processed_drugScreen_data_synId)@filePath, sep="\t")

rawData <- read.csv(synGet('syn6138251')@filePath, sep="\t")


View(rawData)

library("reshape2")

head(data)
d <- dcast(data=rawData, drug ~ cellLine)
rownames(d) <- d$drug
d$drug <- NULL
library(d3heatmap)
d3heatmap(d, cexCol=.7, brush_color="black",
          colors="Spectral")
?d3heatmap
View(d)
?melt
