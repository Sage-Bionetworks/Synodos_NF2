library(synapseClient)
synapseLogin()

synIDs <- c("syn6138237", "syn6138251")
annos <- list(consortium="synodos-NF2", dataType="drugScreen", fileType="tsv", disease = "Neurofibromatosis type II")

sapply(synIDs, function(id){
  temp <- synGet(id, downloadFile = FALSE)
  synSetAnnotations(temp) <- annos
  temp <- synStore(temp, forceVersion = FALSE)
  return(NA)
})