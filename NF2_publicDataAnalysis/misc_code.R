#get all the entities under a project
df <- synQuery("select * from entity where benefactorId == 'syn2347420'")


gse_folders <- df[df$entity.concreteType == "org.sagebionetworks.repo.model.Folder",]
gse_folders = gse_folders[grepl('GSE',gse_folders$entity.name),]


#move the meningioma studies to diff folder
meningioma_folders <- gse_folders[!gse_folders$entity.name %in% c('GSE30563'),]
schwannoma_folders <- gse_folders[gse_folders$entity.name %in% c('GSE30563'),]

meningioma_studies_bucket <- "syn2537450"
schwannoma_studies_bucket <- "syn2537451"


#move folders
move_folder <- function(x,study,new_parentId){
  m <- synGet(x['entity.id'])
  m$properties$parentId = new_parentId
  annotations <- list(study=study)
  synSetAnnotations(m) <- annotations  
  m <- synStore(m)
}

#move meningioma folders
temp_ <- apply(meningioma_folders,1,move_folder,'meningioma',meningioma_studies_bucket)

#move schwannoma folders
temp_ <- apply(schwannoma_folders,1,move_folder,'schwannoma',schwannoma_studies_bucket)


