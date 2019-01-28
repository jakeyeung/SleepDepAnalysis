# TimeToZT and PhaseToHsv need to be defined elsewhere
source("R/functions/TimeShiftFunctions.R")
library(PhaseHSV)

DoPca <- function(dat.long.shift, return.mat.and.pca = FALSE, input.is.long=TRUE, jcenter=TRUE, jscale=FALSE, sampname.cname = "sampname", left.formula = "gene ~ "){
  jform <- as.formula(paste0(left.formula, sampname.cname))
  if (input.is.long){
    dat.mat <- dcast(subset(dat.long.shift), formula = jform, value.var = "exprs")
    id.cname <- colnames(dat.mat)[[1]]  # assume first column is gene or transcript ID 
    rownames(dat.mat) <- dat.mat[[id.cname]]; dat.mat[[id.cname]] <- NULL
  } else {
    dat.mat <- dat.long.shift
  }
  
  # normalize rows
  dat.mat <- t(scale(t(dat.mat), center = jcenter, scale = jscale))
  dat.pca <- prcomp(dat.mat)
  if (return.mat.and.pca){
    return(list(pca=dat.pca, mat=dat.mat))
  } else {
    return(dat.pca)
  }
  return(dat.pca)
}

GetCols <- function(labs){
  times <- as.numeric(sapply(labs, function(x) strsplit(x, split = "_")[[1]][[2]]))
  phases <- sapply(times, TimeToZT, period = 24)
  phases.rad <- phases * 2 * pi / 24
  # color by phases
  cols <- hsv(PhaseToHsv(phases.rad, min.phase = 0, max.phase = 2 * pi), s = 1, v = 1)
  return(cols)
}

GetPchs <- function(labs, by.batch = TRUE){
  batches <- sapply(labs, function(l) strsplit(l, "_")[[1]][[3]])
  pchs <- sapply(batches, function(b){
    if (b == 1) return(".")
    if (b == 2) return("*")
  })
  return(pchs)
}
