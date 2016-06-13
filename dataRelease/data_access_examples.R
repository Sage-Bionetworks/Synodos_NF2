library(synapseClient)
#Login to synapse. username/password can be pre-set. see synapse R client page for details
synapseLogin() 

#get the processed drug screen data
processed_drugScreen_data_synId <- "syn6138237"
data_path <- synGet(processed_drugScreen_data_synId)@filePath
data <- read.csv(data_path, sep="\t")
print(paste('#drugs:', length(unique(data$drug))))


### download metadata for synodos drugs
synodos_drugs_table = 'syn6138291'
synodos_drugs <- synTableQuery('select * from syn6138291')
synodos_drugs <- synodos_drugs@values
