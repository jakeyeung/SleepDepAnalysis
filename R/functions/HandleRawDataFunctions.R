ShiftTime <- function(time, period = 24, shift.even=TRUE, shift.final = 144){
  # NSD has 3 samples, one has 4 samples.
  # shift 3 samples by 24 hrs, leave 4th one alone
  # 
  time[which(time == 0)] <- 24
  for (i in 1:length(time)){
    if (shift.even){
      time[i] <- time[i] - 2 * i * period
    } else {
      time[i] <- time[i] - (2 * i - 1) * period
    }
  }
  return(time + shift.final)
}

ReadProcessMetadata <- function(meta.path, keep.original.samps = FALSE){
  if (missing(meta.path)){
    meta.path <- "data/SDrecovery_design_formatted.txt"
  }
  dat.meta <- read.table(meta.path, header = TRUE, sep = "\t")
  if (keep.original.samps){
    dat.meta$SampleID.orig <- dat.meta$SampleID
  }
  dat.meta$SampleID <- StandardizeColnames(dat.meta$SampleID)
  dat.meta$SD_NSD <- as.character(dat.meta$SD_NSD)
  dat.meta <- AddTimeToMeta(dat.meta)
  dat.meta <- AddSampleToMeta(dat.meta)
  return(dat.meta)
}

StandardizeColnames <- function(cnames, add.suffix=NULL, fix.batch2.cnames=FALSE){
  # If begins with "ZT" change to T
  # If samples begin with ABCDE -> change to 12345
  # if samples begin with ABCDE, it doesnot have N or S, add an N (for normal)
  #
  # add.suffix: add vector for suffix (batch number for batch2)
  # fix.batch2.cnames: TRUE/FALSE: change ZT0E2 -> ZT0E, J70 -> T168, J76 -> T174 distinguish the two batches by adding suffix
  
  cnames <- as.character(cnames)
  badnames.i <- grep("^ZT", cnames)
  

  # fix ZTs first, before changing ZT to T
  if (fix.batch2.cnames){
    # ZT0E2 -> ZT0E
   cnames <- gsub("ZT0E2", "ZT0E", cnames) 
    if (is.null(add.suffix)){
      warning("Changed ZT0E2 -> ZT0E, but add.suffix is null. Distinguish ZT0E from batch1 and batch2 by adding suffix.")
    }
   # J70 -> T168, 7 days after SD, ZT0
   cnames <- gsub("^J70", "T168", cnames) 
   # J76 -> T174, 7 days after SD, ZT6
   cnames <- gsub("^J76", "T174", cnames) 
  }
  
  
  for (i in badnames.i){
    cname <- cnames[i]
    cname <- chartr(old = "ABCDEFGHI", "123456789", cname)
    # insert N into second last position of string
    cname <- insert_str(cname, insert = "N", index = nchar(cname))  # function from StringFunctions.R
    # ZT to T
    cname <- gsub("^ZT", "T", cname)
    cnames[i] <- cname
  }
  # add suffix at end
  if (!is.null(add.suffix)){
    cnames <- paste(cnames, add.suffix, sep = "_")
  }
  return(cnames)
}

SampIDToTime <- function(sampid){
  # T12N1 -> 12 (as numeric)
  # ID begins with T, ends with N or S. Extract everything inbetween
  # coerce to numeric as a check
  
  # str_extract requires stringr library
  time <- str_extract(sampid, "(?<=[A-Z]).*(?=[A-Z])")
  return(as.numeric(time))
}

AddTimeToMeta <- function(dat.meta){
  # from sample ID, grep between two letters (T and [N,S] and get the numbers)
  # handle both T0N1 and T0N1_1 equivalently by splitting by underscore
  sampIDs <- sapply(dat.meta$SampleID, function(s) strsplit(s, "_")[[1]][[1]])
  time <- sapply(sampIDs, SampIDToTime)
  dat.meta$time <- time
  return(dat.meta)
}

AddSampleToMeta <- function(dat.meta){
  # T0N1 -> 1, should to be standardized (no ABCDEs)
  # take last character from sample name.
  # Split sample by "_", in case of T0N1_1
  samps <- sapply(dat.meta$SampleID, function(s){
    s <- strsplit(s, "_")[[1]][[1]]
    sampid <- substr(s, nchar(s), nchar(s))
    return(as.numeric(sampid))
  })
  dat.meta$samp <- samps
  return(dat.meta)
}

LoadHtseq <- function(jpath="data/SDrecovery.htseqcounts.txt", metapath = "data/SDrecovery_design.txt", to.long = FALSE){
  dat.meta <- read.table("data/SDrecovery_design.txt", header = TRUE, sep = "\t")
  
  rownames(dat.htseq) <- EnsemblGene2Gene(rownames(dat.htseq))
  # fix colnames T3N4.counts.txt -> T3N4
  colnames(dat.htseq) <- sapply(colnames(dat.htseq), function(c) strsplit(c, "\\.")[[1]][[1]])
  # make into long
  if (!to.long){
    # annotate ZT and SD and NSD
    return(dat.htseq)
  } else {
    dat.htseq.long <- data.frame(gene = rownames(dat.htseq), 
                                 exprs = unlist(dat.htseq), 
                                 sampname = rep(colnames(dat.seq), each = nrow(dat.seq)))  
    return(dat.htseq.long)
  }
}

DatToLong <- function(dat, dat.meta, use.rownames=FALSE, contains.gene.cname = TRUE, has.transcript=FALSE){
  if (missing(dat.meta)){
    dat.meta <- ReadProcessMetadata()
  }
  # expects $gene and all others are samples
  if (contains.gene.cname){
    genes <- dat$gene
    if (!has.transcript){
      exprs <- subset(dat, select = -(gene))
    } else {
      exprs <- subset(dat, select = -c(gene, transcript))
    }
  } else {
    genes <- rownames(dat)
    exprs <- dat
  }
  if (has.transcript){
    trans <- dat$transcript
  } else {
    trans <- NULL
  }
  
  time.hash <- hash(as.character(dat.meta$SampleID), as.numeric(dat.meta$time))
  trtmnt.hash <- hash(as.character(dat.meta$SampleID), dat.meta$SD_NSD)
  samp.hash <- hash(as.character(dat.meta$SampleID), as.numeric(dat.meta$samp))
  
  time <- sapply(dat.meta$SampleID, function(s) time.hash[[s]])
  trtmnt <- sapply(dat.meta$SampleID, function(s) trtmnt.hash[[s]])
  samp <- sapply(dat.meta$SampleID, function(s) samp.hash[[s]])
  
  if (has.transcript){
    dat.long <- data.frame(gene = genes,
                           transcript = trans,
                           exprs = unlist(exprs),
                           samp = rep(samp, each = nrow(exprs)),
                           time = rep(time, each = nrow(exprs)),
                           trtmnt = rep(trtmnt, each = nrow(exprs)))
  } else {
    dat.long <- data.frame(gene = genes,
                           exprs = unlist(exprs),
                           samp = rep(samp, each = nrow(exprs)),
                           time = rep(time, each = nrow(exprs)),
                           trtmnt = rep(trtmnt, each = nrow(exprs)))
  }
  if (!is.null(dat.meta$Batch)){
    batch <- dat.meta$Batch
    dat.long$batch <- rep(batch, each = nrow(exprs))
  }
  return(dat.long)
}



RemoveMtGenes <- function(dat.long){
  dat.long <- dat.long[!grepl("^mt", dat.long$gene), ]
  return(dat.long)
}

Make48Hours <- function(dat.sub){
  # NSD has 4 replicates so lets split into 2 replicates and do an average 
  # to get 0 to 42 hours
  warning("Not implemented")
}

MakeNegHours <- function(time, trtmnt, reverse=FALSE){
  # Make NSD to negative hours
  if (reverse){
    if (trtmnt == "NSD" & time <= 0){
      time <- time + 24
    } else {
      time <- time
    }
  } else {
    if (trtmnt == "NSD" & time >= 0){
      time <- time - 24
    } else {
      time <- time
    } 
  }
  return(time)
}

SpreadRepsOverTime <- function(dat.sub, triple.time = c(0, 6)){
  # Take "NSD" treated mice and spread the 4 replicates over time.range (if triple.time = 6, then 6, 30, 54 hours)
  # triple.time: time at which we should get triples over 4 replicates (take 2 replicates at time ZT6 to be same at ZT54)
  # check trtmnt
  
#   trtmnt <- unique(dat.sub$trtmnt)
#   if (length(trtmnt) > 1) warning("Must handle a single treatment")
  
#   time <- unique(dat.sub$time)
#   if (length(time) > 1) warning("Must handle a single time point")
  n.reps <- nrow(dat.sub)
  time <- dat.sub$time[[1]]
#   if (trtmnt == "NSD"){

  if (time %in% triple.time){
    # take last two replicates and add 24 to each 
    indx <- (n.reps - 1) : n.reps
  } else {
    # take last replicate and add 24 hours
    indx <- n.reps
  }
  for (i in indx){
    dat.sub$time[i] <- dat.sub$time[i - 1] + 24
  }
  return(dat.sub)
}

DoubleRepsOverTime <- function(dat.long, triple.time = c(0, 6)){
  # Instead of spreading 4 replicates into 2, just double plot the 4 replicates
  # triple time: triple up 0 24 48 or 6 30 54
  dat.sub <- subset(dat.long, trtmnt == "NSD")   # needs doubling
  dat.sd <- subset(dat.long, trtmnt == "SD")  # keep 
  
  dat.sub$time <- as.numeric(as.character(dat.sub$time))
  dat.trip <- subset(dat.sub, time %in% triple.time)
  dat.trip$time <- subset(dat.sub, time %in% triple.time)$time + 48
  dat.dble <- dat.sub
  dat.dble$time <- dat.sub$time + 24
  
  dat.nsd <- rbind(dat.sub, dat.dble, dat.trip)
  dat.dbl <- rbind(dat.nsd, dat.sd)
  return(dat.dbl)
}

DuplicateTime <- function(dat.sub, time.old, time.new){
  # Duplicate ZT0 to ZT24
  dat.dbl <- subset(dat.sub, time == time.old)
  dat.dbl$time <- time.new
  return(rbind(dat.sub, dat.dbl))
}

ReadMetaMatchSampOrder <- function(metapath, dat.sampnames){
  # Read meta path, rearrange rows to match order of dat.sampnames of expression matrix
  dat.meta <- read.table(metapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # order by target vector:
  # http://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-a-target-vector-that-specifies-the-desired-or
  order.i <- match(dat.sampnames, dat.meta$SampleID)
  return(dat.meta[order.i, ])
}

ReadMetaDataRound2 <- function(metapath){
  if (missing(metapath)){
    metapath <- "/home/yeung/data/sleep_deprivation/SDrecovery_design_formatted.round2.txt"
  }
  dat.meta <- read.table(metapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  jadd.suffix <- dat.meta$Batch
  jfix.batch2.cnames <- TRUE
  
  dat.meta$OrigSampID <- dat.meta$SampleID
  dat.meta$SampleID <- StandardizeColnames(dat.meta$SampleID, jadd.suffix, jfix.batch2.cnames)
  dat.meta$SD_NSD <- as.character(dat.meta$SD_NSD)
  dat.meta <- AddTimeToMeta(dat.meta)
  dat.meta <- AddSampleToMeta(dat.meta)
  
  return(dat.meta)
}

LoadProcessData <- function(inpath, metapath, ftype = "kallisto", has.transcript=FALSE, make.neg.hours=TRUE, return.edata.and.meta=FALSE){
  if (missing(metapath) & ftype != "kallisto2"){
    metapath <- "/home/yeung/projects/sleep_deprivation/data/SDrecovery_design_formatted.txt"
  } else if (missing(metapath) & ftype == "kallisto2"){
    metapath <- "/home/yeung/data/sleep_deprivation/SDrecovery_design_formatted.round2.txt"
  }
  if (missing(inpath)){
    # infer path from ftype
    if (ftype == "kallisto"){
      inpath <- "/home/yeung/projects/sleep_deprivation/data/SDrecovery_kallisto_tpm.txt"  
      dat <- read.table(inpath, header = TRUE, sep = "\t")
    } else if (ftype == "htseq"){
      inpath <- "data/SDrecovery.htseqcounts.txt"
      dat <- read.table(inpath)
    } else if (ftype == "kallisto2"){
      # kallisto round two 2017-02-20
      inpath <- "/home/yeung/data/sleep_deprivation/kallisto_round_two/SDrecovery_kallisto_tpm_round2.txt"
      dat <- read.table(inpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else {
      warning("ftype must be kallisto, kallisto2, or htseq")
    }
  } else {
    dat <- read.table(inpath)
  }
  dat.meta <- ReadMetaMatchSampOrder(metapath, dat.sampnames = colnames(dat)[-1])  # no first column
  # fix samp colnames differently if round2 samps vs in round1
  if (ftype == "kallisto2"){
    jadd.suffix <- dat.meta$Batch
    jfix.batch2.cnames <- TRUE
  } else {
    jadd.suffix <- NULL
    jfix.batch2.cnames <- FALSE
  }
  if (return.edata.and.meta){
    return(list(edata = dat, meta=dat.meta))
  }
  # if (ftype == "kallisto2"){
  #   dat.meta <- read.table(metapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # } else {
  #   dat.meta <- read.table(metapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # }
  dat.meta$SampleID <- StandardizeColnames(dat.meta$SampleID, jadd.suffix, jfix.batch2.cnames)
  dat.meta$SD_NSD <- as.character(dat.meta$SD_NSD)
  dat.meta <- AddTimeToMeta(dat.meta)
  dat.meta <- AddSampleToMeta(dat.meta)
  print(dat.meta)
  if (ftype == "kallisto" | ftype == "kallisto2"){
    gene <- Transcript2Gene(dat$transcript)
    dat$gene <- gene
    colnames(dat) <- StandardizeColnames(colnames(dat))
    if (!has.transcript){
      dat$transcript <- NULL
    }
  } else if (ftype == "htseq"){
    gene.htseq <- EnsemblGene2Gene(rownames(dat))
    dat$gene <- gene.htseq
    # fix colnames
    colnames(dat) <- sapply(colnames(dat), function(cname) strsplit(cname, "\\.")[[1]][[1]])
    colnames(dat) <- StandardizeColnames(colnames(dat))
  } else {
    warning("ftype must be kallisto, kallisto2, or htseq")
  }
  dat.long <- DatToLong(dat, dat.meta, has.transcript = has.transcript)
  if (make.neg.hours){
    dat.long$time <- mapply(MakeNegHours, dat.long$time, dat.long$trtmnt)
  }
  return(dat.long)
}

SumOverSamples <- function(dat.long, eps = 0.1, remove.mt.genes = TRUE, has.transcript=FALSE){
  # sum over samples and do log2 in mean, variance 
  if (has.transcript){
    dat.long <- dat.long %>%
        group_by(gene, time, trtmnt, transcript) %>%  
      summarise(exprs.mean = mean(exprs), exprs.var = var(exprs), 
                exprs.linear = exprs.mean, 
                exprs.var.log2 = var(log2(exprs + eps)),
                exprs = mean(log2(exprs + eps)))
  } else {
    dat.long <- dat.long %>%
        group_by(gene, time, trtmnt) %>%
      summarise(exprs.mean = mean(exprs), exprs.var = var(exprs), 
                exprs.linear = exprs.mean, 
                exprs.var.log2 = var(log2(exprs + eps)),
                exprs = mean(log2(exprs + eps)))
    
  }
  if (remove.mt.genes){
    dat.long <- RemoveMtGenes(dat.long)
  }
  return(dat.long)
}
