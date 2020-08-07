WrangleDatForModelSelection <- function(dat.long.shift, zt0tozt24 = FALSE, removezt78 = TRUE){
  dat.long.shift <- dat.long.shift[order(dat.long.shift$time), ]
  dat.long.shift$exprs <- dat.long.shift$log2exprs; dat.long.shift$log2exprs <- NULL
  if (removezt78){
    dat.long.shift <- subset(dat.long.shift, time != 78)
  }
  if (zt0tozt24){
    # change ZT0 to ZT24
    dat.long.shift$time <- sapply(dat.long.shift$time, function(tt){
      if (tt == 0){
        return(24)
      } else {
        return(tt)
      }
    })
  }  
  return(dat.long.shift)
}

# MakeDatShift moved to TimeShiftFunctions

