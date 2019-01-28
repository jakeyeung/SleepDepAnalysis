LoadSleuthOutput <- function(inf){
  load(inf, v=T)
  # label Time
  bname <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(inf))
  jtime <- strsplit(bname, "_")[[1]][[3]]
  jtime <- as.numeric(gsub("[^\\d]+", "", jtime, perl=TRUE))
  sleuth_table$time.nsd <- jtime
  full_model$time.nsd <- jtime
  return(list(so=sleuth_table, model=full_model))
}

LoadSleuthOutput.cleaner <- function(inf){
  # cleaner includes sleuth_table.wt
  load(inf, v=T)
  # label Time
  bname <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(inf))
  jtime <- strsplit(bname, "_")[[1]][[3]]
  jtime <- as.numeric(gsub("[^\\d]+", "", jtime, perl=TRUE))
  sleuth_table$time.nsd <- jtime
  full_model$time.nsd <- jtime
  sleuth_table.wt$time.nsd <- jtime
  return(list(sleuth_table=sleuth_table, full_model=full_model, sleuth_table.wt=sleuth_table.wt, so=so))
}