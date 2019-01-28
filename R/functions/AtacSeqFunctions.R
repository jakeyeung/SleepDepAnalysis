CompareTwoTimes <- function(dat){
  test <- t.test(exprs ~ time, dat)
  out.df <- data.frame(statistic = test$statistic, pval = test$p.value, difference = diff(test$estimate))
  return(out.df)
}

SampToTime.merged <- function(s, time.prefix = "ZT"){
  # ZT03_NSD -> 3
  return(as.integer(gsub("ZT", "", strsplit(s, "_")[[1]][[1]])))
}

SampToTrt.merged <- function(s){
  # ZT03_NSD -> "NSD"
  return(strsplit(s, "_")[[1]][[2]])
}

SampToTime <- function(samp){
  # "ZT03_NSD_2" -> 3
  jtime <- strsplit(samp, "_")[[1]][[1]]
  jtime <- gsub("ZT", "", x = jtime)
  return(as.numeric(jtime))
}

SampToTrtmnt <- function(samp){
  # "ZT03_NSD_2" -> NSD
  trtmnt <- strsplit(samp, "_")[[1]][[2]]
  return(trtmnt)
}

SampToRep <- function(samp){
  # "ZT03_NSD_2" -> NSD
  jrep <- strsplit(samp, "_")[[1]][[3]]
  return(jrep)
}