extractGEOdata <- function(f){
  if(file_ext(f) == "tar"){
    untar(f, exdir=dirname(f))
  }
}


get_GEO_ftpLink <- function(GEO){
  geotype <- toupper(substr(GEO, 1, 3))
  stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
  if (geotype == "GSM") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/", 
                   stub, GEO)
  } else if (geotype == "GSE") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/", 
                   stub, GEO)
  } else if (geotype == "GPL") {
    url <- sprintf("ftp://ftp.ncbi.nlm.nih.gov/geo/platform/%s/%s/suppl/", 
                   stub, GEO)
  } else {
    url <- None
  }
}

get_formatted_phenoData <- function(geo_study){
  df <- pData(geo_study)
  idx <- grep("characteristics", colnames(df), fixed=T, value=T)
  df <- df[,idx,drop=F]
  for(i in idx){
    tmp <- strsplit(as.character(df[[i]]), ": ", fixed=T)
    nm <- sapply(tmp, "[[", 1)
    df[[unique(nm)]] <- sapply(tmp, "[[", 2)
    df[[i]] <- NULL
  }
  df['sampleId'] <- rownames(df)
  return(df)
}
