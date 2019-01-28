MakeDatShift <- function(dat.long){
  dat.nsd <- subset(dat.long, trtmnt == "NSD")
  dat.sd.shift <- subset(dat.long, trtmnt == "SD") %>%
    group_by(time, gene) %>%
    mutate(time.shift = time + 24)
  dat.sd.shift$time <- dat.sd.shift$time.shift; dat.sd.shift$time.shift <- NULL
  dat.long.shift <- rbind(dat.nsd, dat.sd.shift)
  return(dat.long.shift)
}

MakeDatShift.nogene <- function(dat.long){
  # same wihtout gene column 
  dat.long <- as.data.frame(dat.long)
  dat.nsd <- subset(dat.long, trtmnt == "NSD")
  dat.sd.shift <- subset(dat.long, trtmnt == "SD") %>%
    group_by(time) %>%
    mutate(time.shift = time + 24)
  dat.sd.shift$time <- dat.sd.shift$time.shift; dat.sd.shift$time.shift <- NULL
  dat.long.shift <- rbind(as.data.frame(dat.nsd), as.data.frame(dat.sd.shift))
  return(dat.long.shift)
}

RenameTime <- function(dat.long, time.old, time.new, trt = "NSD", increment.samp = 0){
  # Rename timepoint, optionally increment sample if time.new already
  # has samples to prevent identical samples in new timepoints
  #
  # Replaces time only for specified trt
  row.i <- which(dat.long$time == time.old & dat.long$trtmnt == trt)
  dat.long$time[row.i] <- time.new
  # optionally increment samp
  if (increment.samp > 0){
    dat.long$samp[row.i] <- dat.long$samp[row.i] + increment.samp
  }
  return(dat.long)
}

TimeToZT <- function(x, period=24){
  dayshift <- floor(x / period)
  hourshift <- dayshift * period
  return(x - hourshift)
}

TimeToDay <- function(x, period=24){
  dayshift <- floor(x / period)
}
