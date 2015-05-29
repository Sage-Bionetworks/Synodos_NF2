#mode read delly code for intansv package



dellyCluster <- function(df) 
{
  maxReadPairSupp <- max(df$rp_support)
  dfFil <- df[df$rp_support>=(maxReadPairSupp/2), ]
  dfFilIrange <- IRanges(start=dfFil$start, end=dfFil$end)
  outTmp <- findOverlaps(dfFilIrange, reduce(dfFilIrange))
  dfFil$cludelly <- subjectHits(outTmp)
  dfFilRes <- ddply(dfFil, ("cludelly"), function(x){
    if(nrow(x)==1){
      return(x)
    } else {
      LeftMin <- min(x$start)
      RightMax <- max(x$end)
      RangeLength <- RightMax-LeftMin
      x$op <- (x$end - x$start)/RangeLength
      if(any(x$op<0.8)){
        return(NULL)
      } else {
        return(x[which.max(x$rp_support), ])
      }
    }
  })
}


readdelly_local <- function (dataDir = ".", regSizeLowerCutoff = 100, regSizeUpperCutoff = 1e+06, 
                             readsSupport = 3, method = "DELLY", pass = TRUE, minMappingQuality = 30) 
{
  dellyFileList <- list.files(dataDir, full.names = T)
  dellyFileCont <- lapply(dellyFileList, function(x) {
    dellyData <- try(read.table(x, as.is = T), silent = T)
    if (is.data.frame(dellyData)) {
      dellyData <- dellyData[, c(1:2, 7:8, 10)]
      names(dellyData) <- c("chr", "start", "pass", "info", 
                            "detail")
      return(dellyData)
    }
    else {
      return(NULL)
    }
  })
  dellyCont <- do.call(rbind, dellyFileCont)
  if (pass) {
    dellyCont <- dellyCont[dellyCont$pass == "PASS", ]
  }
  dellyCont$start <-  as.numeric(dellyCont$start)
  dellyCont$end <- as.numeric(gsub(".+;END=(\\d+);.+", "\\1", 
                                   dellyCont$info))
  dellyCont$type <- gsub(".+;SVTYPE=([A-Z]+);.+", "\\1", dellyCont$info)
  dellyCont$rp_support <- as.numeric(gsub(".+;PE=(\\d+);.+", 
                                          "\\1", dellyCont$info))
  dellyCont$map_quality <- as.numeric(gsub(".+;MAPQ=(\\d+).*", 
                                           "\\1", dellyCont$info))
  dellyCont$size <- dellyCont$end - dellyCont$start
  
  dellyCont <- dellyCont[dellyCont$map_quality >= minMappingQuality & 
                           dellyCont$rp_support >= readsSupport & dellyCont$size >= 
                           regSizeLowerCutoff & dellyCont$size <= regSizeUpperCutoff, 
                         ]
  dellyCont <- dellyCont[, c("chr", "start", "end", "size", 
                             "rp_support", "map_quality", "type")]
  dellyDelDf <- dellyCont[dellyCont$type == "DEL", ]
  if (!is.null(dellyDelDf)) {
    dellyDelIrange <- GRanges(seqnames = dellyDelDf$chr, 
                              ranges = IRanges(start = dellyDelDf$start, end = dellyDelDf$end))
    dellyDelIrangeRes <- findOverlaps(dellyDelIrange, reduce(dellyDelIrange))
    dellyDelDf$clu <- subjectHits(dellyDelIrangeRes)
    dellyDelDfFilMer <- ddply(dellyDelDf, ("clu"), dellyCluster)
    dellyDelDfFilMer <- dellyDelDfFilMer[, c("chr", "start", 
                                             "end", "size")]
    names(dellyDelDfFilMer) <- c("chromosome", "pos1", "pos2", 
                                 "size")
  }
  else {
    dellyDelDfFilMer <- NULL
  }
  dellyDupDf <- dellyCont[dellyCont$type == "DUP", ]
  if (!is.null(dellyDupDf)) {
    dellyDupIrange <- GRanges(seqnames = dellyDupDf$chr, 
                              ranges = IRanges(start = dellyDupDf$start, end = dellyDupDf$end))
    dellyDupIrangeRes <- findOverlaps(dellyDupIrange, reduce(dellyDupIrange))
    dellyDupDf$clu <- subjectHits(dellyDupIrangeRes)
    dellyDupDfFilMer <- ddply(dellyDupDf, ("clu"), dellyCluster)
    dellyDupDfFilMer <- dellyDupDfFilMer[, c("chr", "start", 
                                             "end", "size")]
    names(dellyDupDfFilMer) <- c("chromosome", "pos1", "pos2", 
                                 "size")
  }
  else {
    dellyDupDfFilMer <- NULL
  }
  dellyInvDf <- dellyCont[dellyCont$type == "INV", ]
  if (!is.null(dellyInvDf)) {
    dellyInvIrange <- GRanges(seqnames = dellyInvDf$chr, 
                              ranges = IRanges(start = dellyInvDf$start, end = dellyInvDf$end))
    dellyInvIrangeRes <- findOverlaps(dellyInvIrange, reduce(dellyInvIrange))
    dellyInvDf$clu <- subjectHits(dellyInvIrangeRes)
    dellyInvDfFilMer <- ddply(dellyInvDf, ("clu"), dellyCluster)
    dellyInvDfFilMer <- dellyInvDfFilMer[, c("chr", "start", 
                                             "end", "size")]
    names(dellyInvDfFilMer) <- c("chromosome", "pos1", "pos2", 
                                 "size")
  }
  else {
    dellyInvDfFilMer <- NULL
  }
  retuRes <- list(del = dellyDelDfFilMer, dup = dellyDupDfFilMer, 
                  inv = dellyInvDfFilMer)
  attributes(retuRes) <- c(attributes(retuRes), list(method = method))
  return(retuRes)
}



