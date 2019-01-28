# Jake Yeung
# Date of Creation: 2017-09-13
# File: ~/projects/sleep_deprivation/scripts/functions/MergingAssigningSamples.R
# Functions for merging technical replicates and assigning samples to time24

MergeRep <- function(dat.sub){
  # merge rows and make new SampID as concatenation of the samples
  # new batch name, take minimum of the batches
  if(nrow(dat.sub) == 1){
    warning("Input is of row 1, returning same df")
    return(dat.sub)
  }
  gene.new <- unique(dat.sub$gene)
  sampid.new <- paste(unique(dat.sub$SampID), collapse = ",")
  exprs.new <- mean(dat.sub$exprs)
  trtmnt.new <- unique(dat.sub$trtmnt)
  time.new <- unique(dat.sub$time)
  batch.new <- paste(unique(dat.sub$batch), collapse = ",")
  samp.new <- unique(dat.sub$samp)
  
  lengths <-sapply(list(gene.new, exprs.new, trtmnt.new, samp.new, batch.new, samp.new), length)
  
  if(any(lengths > 1)){
    warning("Found length greater than 1")
  }
  dat.merged <- data.frame(gene = gene.new,
                           SampID = sampid.new,
                           exprs = exprs.new,
                           trtmnt = trtmnt.new,
                           time = time.new,
                           batch = batch.new,
                           samp = samp.new)
  return(dat.merged)
}

MergeTechnicalReplicates <- function(dat.long){
  # For each gene, time, trtmnt, merge technical replicates (same samp, different batch)
  # in Round 2: Time 0 contains samp 5 which is in both batch 1 and batch 2.
  # Take average exprs between the two samps
  # 
  # This function is not general, for Round 2 only.
  dat.nochange <- subset(dat.long, time != 0 | samp != 5 | trtmnt != "NSD")
  dat.tochange <- subset(dat.long, time == 0 & samp == 5 & trtmnt == "NSD")
  
  # handle zt0 separately, then merge back into dataframe
  dat.merged <- dat.tochange %>%
    group_by(gene, trtmnt, time, samp) %>%
    do(MergeRep(.))
  dat.new <- bind_rows(dat.merged, dat.nochange)
  return(dat.new)
}

AssignToZT24 <- function(dat.sub, samps.to.zt24){
  # dat.sub should be all same time, trtmnt, but different samps
  # If samp s in samps.to.zt24, move to time 24
  if (any(dat.sub$time != 0)){
    warning("Function not written for times not 0")
  }
  row.i <- dat.sub$samp %in% samps.to.zt24
  dat.sub$time[row.i] <- 24
  return(dat.sub)
}

AssignZT0ToZT24Batching <- function(dat.long){
  # Split ZT0 to ZT24s evenly and deterministically
  # Randomize afterwards as an option
  dat.nochange <- subset(dat.long, time != 0 | trtmnt != "NSD")
  dat.tochange <- subset(dat.long, time == 0 & trtmnt == "NSD")
  
  samps.to.zt24 <- c(3, 5, 9, 72)  # these samps at ZT0 get send to ZT24
  
  dat.changed <- dat.tochange %>%
    group_by(gene, exprs, trtmnt, time) %>%
    do(AssignToZT24(., samps.to.zt24))
  
  dat.merged <- bind_rows(dat.changed, dat.nochange)
  return(dat.merged)
}

AssignZT0ToZT24Batching.nogene <- function(dat.long){
  # Split ZT0 to ZT24s evenly and deterministically
  # Randomize afterwards as an option
  # identical to AssignZT0ToZT24Batching() but no gene. Used to find samples that were asigned to 0 or 24
  dat.nochange <- subset(dat.long, time != 0 | trtmnt != "NSD")
  dat.tochange <- subset(dat.long, time == 0 & trtmnt == "NSD")
  
  samps.to.zt24 <- c(3, 5, 9, 72)  # these samps at ZT0 get send to ZT24
  
  dat.changed <- dat.tochange %>%
    group_by(trtmnt, time) %>%
    do(AssignToZT24(., samps.to.zt24))
  
  dat.merged <- bind_rows(dat.changed, dat.nochange)
  return(dat.merged)
}
