# Jake Yeung
# Date of Creation: 2018-01-10
# File: ~/projects/sleep_deprivation/scripts/functions/LoadData.R
# Load data from primetime, get the necessary 

LoadInits <- function(){
  # rename cluster according to Charlotte email
  LoadClstrConstants()
  
  prefilt.sleuth <- TRUE
  counts.cutoff <- 5.5  # arbitrary
  
  if (!prefilt.sleuth){
    plotmain <- paste0("/home/yeung/projects/sleep_deprivation/plots/primetime.cutoff.", counts.cutoff)
  } else {
    plotmain <- paste0("/home/yeung/projects/sleep_deprivation/plots/primetime.sleuth.cutoff.", counts.cutoff)
  }
  plotmain <- paste0(plotmain, "_", Sys.Date())
  dir.create(plotmain)
  
  
  min.time <- 0
  max.time <- 78
  downsamp <- 0.005
  tstep <- 4/3600  # 4seconds in units of hour
  do.linear <- FALSE
  
  low.pass.filt.time <- rep(0.1, length = 780) * 1:780
  # low.pass.filt.time <- c(0, low.pass.filt.time)  # causes bug when runnig PlotFit() "S.out did not get filtered properly into evenly spaced timepoints"
  filt.time <- seq(min.time, max.time)
  filt.time <- filt.time[!filt.time %in% c(31, 32)]
  
  var.lst <- list(clstr.hash, clstr.order, prefilt.sleuth, counts.cutoff, 
                  plotmain, min.time, max.time, tstep, do.linear, low.pass.filt.time, filt.time)
  names(var.lst) <- c("clstr.hash", "clstr.order", "prefilt.sleuth", "counts.cutoff", "plotmain", "min.time", "max.time", "tstep", "do.linear", "low.pass.filt.time", "filt.time")
  list2env(var.lst, envir = globalenv())
  return(NULL)
}


GetCountsMatrix <- function(prefilt.sleuth, max.time, filter.max.time = TRUE){
  if (!prefilt.sleuth){
    load("Robjs/combat/dat.long.cleaned.techmerged_zt24assigned.Robj", v=T)
    load("data/dat.long.cleaned.techmerged_zt24assigned.Robj", v=T)
    # dat.long.shift <- dat.long.cleaned; rm(dat.long.cleaned)
    # do linear!
    do.linear <- FALSE
    if (do.linear){
      dat.long.cleaned$exprs <- 2^dat.long.cleaned$exprs - 1
    }
    dat.long.cleaned <- dat.long.cleaned %>%
      group_by(gene, trtmnt) %>%
      arrange(time)
    # dat.long.shift <- subset(dat.long.shift, time <= max.time) do it within IF loop
  } else {
    load("/home/yeung/projects/sleep_deprivation/Robjs/combat_sleuth_filt/dat_combat_countsfilt.5.5.Robj", v=T)
  }
  # ignore timepoints: 192 and 198 SD???
  if (filter.max.time){
    dat.long.shift <- subset(dat.long.cleaned, time <= max.time) %>%
      group_by(gene) %>%
      arrange(time)
  } else {
    dat.long.shift <- dat.long.cleaned
  }
  return(dat.long.shift)
}

LoadClstrConstants <- function(){
  clstr.hash <- hash(c(4, 8, 3, 9, 6, 1, 2, 10, 5, 7),
                     c(1, 2, 3, 4, 9, 10, 5, 6, 7, 8))
  clstr.order <- c(1, 5, 2, 6, 3, 7, 4, 8, 9, 10) 
  var.lst <- list(clstr.hash, clstr.order)
  names(var.lst) <- c("clstr.hash", "clstr.order")
  list2env(var.lst, envir = globalenv())
}

GetGenesKeep <- function(prefilt.sleuth, counts.cutoff, fits){
  if (!prefilt.sleuth){
    # Necessary if running on TPM-filt counts
    # inf <- "/home/yeung/data/sleep_deprivation/sleuth_outputs.norm_all/sleuth_norm_all.Rdata"
    inf <- "data/sleuth_norm_all.Rdata"
    so <- sleuth_load(inf)
    sr <- sleuth_results(so, 'reduced:full', 'lrt')
    t2g <- hash(sr$target_id, sr$ext_gene)
    mat.long <- kallisto_table(so) %>%
      mutate(est_counts.log2 = so$transform_fun(est_counts),
             tpm.log2 = so$transform_fun(tpm))
    mat.long$gene <- sapply(mat.long$target_id, function(tx) t2g[[tx]])
    mat.long <- mat.long %>%
      dplyr::filter(!is.na(gene))
    mat.long.gene <- mat.long %>%
      group_by(gene, sample, ZT, trtmnt, batch, samp.raw, time, samp) %>%
      summarise(est_counts = sum(est_counts),
                tpm = sum(tpm)) %>%
      mutate(est_counts.log2 = so$transform_fun(est_counts),
             tpm.log2 = so$transform_fun(tpm))
    mat.gene.sum <- mat.long.gene %>%
      group_by(gene) %>%
      # summarise(est_counts.log2.min = quantile(est_counts.log2, probs = 0.25))
      summarise(counts.smry = mean(est_counts.log2))
    genes.keep <- as.character(subset(mat.gene.sum, counts.smry > counts.cutoff)$gene)
  } else {
    # If counts-filt genes, then just use all genes from fits
    genes.keep <- fits$gene  # take all
  }
  return(genes.keep)
}

LoadPCAMat <- function(dat.long.cleaned){
  # Load dat for PCA --------------------------------------------------------
  
  if (!prefilt.sleuth){
    # Old
    # load("/home/yeung/projects/sleep_deprivation/Robjs/combat/dat.long.kallisto.combat.mean_only.FALSE.exprs_cutoff.2.5.Robj"); dat.after <- dat.long.c
    load("data/dat.long.kallisto.combat.mean_only.FALSE.exprs_cutoff.2.5.Robj"); dat.after <- dat.long.c
  } else {
    # New with counts filt
    dat.after <- dat.long.cleaned  ## unfiltered 
  }
  return(dat.after)
  # assign("dat.after", dat.after, envir = globalenv())
}

LoadClusteringObjs <- function(){
  # Load Clustering output. Kmeans sometimes assigns clusters differently so if we reorder clusters by hard coding it is problematic
  load("data/kmeans_objs.Robj", verbose=TRUE)
  vars.lst <- list(dat.long.means, clstrs, cent.long.base, cent.long.full, mat.means.exprs, mat.means, wss, max.n, n.centers)
  names(vars.lst) <- c("dat.long.means", "clstrs", "cent.long.base", "cent.long.full", "mat.means.exprs", "mat.means", "wss", "max.n", "n.centers")
  list2env(vars.lst, envir = globalenv())
}

LoadClusteringObjs.runKmeans <- function(dat.long.shift, clstr.hash, clstr.order, n.centers = 10, qval.cutoff = 1e-3){
  # n.centers <- 10
  # qval.cutoff <- 1e-3
  
  # Fancier diff analysis with shrinkage
  # inf <- "/home/yeung/data/sleep_deprivation/sleuth_outputs.norm_all.gene_aggreg/sleuth_norm_all.aggreg_gene.TRUE.Rdata"
  inf <- "data/sleuth_norm_all.aggreg_gene.TRUE.Rdata"
  so <- sleuth_load(inf)
  sr <- sleuth_results(so, 'reduced:full', 'lrt')
  sr <- dplyr::rename(sr, gene = "target_id")
  sr <- subset(sr, select = c(-transcript_version, -ens_gene, -description, -transcript_biotype))
  sr <- sr %>%
    group_by(gene) %>%
    filter(rank(test_stat, ties.method="first")==1)  # https://stackoverflow.com/questions/21308436/dplyr-filter-get-rows-with-minimum-of-variable-but-only-the-first-if-multiple
  print(sr)
  
  # they should also be sufficiently expressed!
  signif.genes <- unique(subset(sr, qval < qval.cutoff)$gene)
  
  dat.long.shift.hits <- subset(dat.long.shift, gene %in% signif.genes)
  
  dat.long.means <- dat.long.shift.hits %>%
    group_by(gene, time) %>%
    summarise(exprs = mean(exprs))
  
  do.normalize <- TRUE
  if (do.normalize){
    dat.long.means <- dat.long.means %>%
      group_by(gene) %>%
      mutate(exprs = scale(exprs, center = TRUE, scale = TRUE))
  }
  
  mat.means <- dcast(dat.long.means, formula = gene ~ time, value.var = "exprs")
  mat.means.exprs <- as.matrix(mat.means[, !grepl("gene", colnames(mat.means))])
  
  # add pval from LRT analysis
  mat.means <- dplyr::inner_join(mat.means, subset(sr, select = c(gene, pval, qval)))
  
  # Plot within groups sum of squares over k
  max.n <- 30
  wss <- (nrow(mat.means.exprs)-1)*sum(apply(mat.means.exprs,2,var))
  for (i in 2:max.n) wss[i] <- sum(kmeans(mat.means.exprs, centers=i)$withinss)
  
  clstrs <- kmeans(mat.means.exprs, centers = n.centers)
  times.vec <- as.numeric(colnames(mat.means.exprs))
  
  # Order clusters in a meaningful way
  
  hclst <- hclust(dist(as.matrix(clstrs$centers)))
  
  # reorder clstrs and then rename 
  cc.reord <- clstrs$centers[hclst$order, ]
  reord.hash <- hash(rownames(cc.reord), as.character(seq(nrow(cc.reord))))
  
  # rename clusters
  clstrs$cluster <- sapply(as.character(clstrs$cluster), function(x) reord.hash[[x]])
  rownames(clstrs$centers) <- sapply(as.character(rownames(clstrs$centers)), function(x) reord.hash[[x]])
  
  # Plot all clusters
  clstr.gene.map <- hash(mat.means$gene, clstrs$cluster)
  dat.long.means$clstr <- as.numeric(sapply(dat.long.means$gene, function(g) clstr.gene.map[[g]]))
  
  cent.long.full <- data.frame(clstr = as.numeric(rownames(clstrs$centers)),
                               time = as.numeric(rep(colnames(clstrs$centers), each = nrow(clstrs$centers))),
                               exprs = unlist(as.data.frame(clstrs$centers), use.names = FALSE), 
                               stringsAsFactors = FALSE)
  cent.long <- subset(cent.long.full, time < 24)  # take first day and duplicate 2 times (3 day total)
  cent.long.base <- dplyr::bind_rows(cent.long, cent.long %>% mutate(time = time + 24), cent.long %>% mutate(time = time + 48))
  # add NSD for time 72 and 78 (ZT0 and ZT6 day 4)
  cent.long.base <- dplyr::bind_rows(cent.long.base, cent.long %>% dplyr::filter(time <= 6) %>% mutate(time = time + 72))
  
  # rename clusters according to Charlotte's email 2017-12017-10-18
  
  dat.long.means$clstr <- RearrangeCluster(dat.long.means$clstr, clstr.hash, clstr.order)
  cent.long.full$clstr <- RearrangeCluster(cent.long.full$clstr, clstr.hash, clstr.order)
  cent.long.base$clstr <- RearrangeCluster(cent.long.base$clstr, clstr.hash, clstr.order)
  clstrs$cluster <- RearrangeCluster(clstrs$cluster, clstr.hash, clstr.order = NA)
  rownames(clstrs$centers) <- RearrangeCluster(rownames(clstrs$centers), clstr.hash, clstr.order = NA)
  
  cent.long.full$clstr.num.order <- factor(as.character(cent.long.full$clstr), levels = sort(as.numeric(unique(cent.long.full$clstr))))
  
  mat.means$cluster <- clstrs$cluster
  vars.lst <- list(dat.long.means, clstrs, cent.long.base, cent.long.full, mat.means.exprs, mat.means, wss, max.n, n.centers)
  names(vars.lst) <- c("dat.long.means", "clstrs", "cent.long.base", "cent.long.full", "mat.means.exprs", "mat.means", "wss", "max.n", "n.centers")
  list2env(vars.lst, envir = globalenv())
  return(NULL)
}

LoadPrimetimeObjs <- function(){
  # setwd("/home/yeung/projects/sleep_deprivation")
  
  library(dplyr)
  library(reshape2)
  library(wordcloud)
  library(hash)
  library(ggplot2)
  library(sleuth)
  
  source("R/functions/RoundTwoPCAScripts.R")
  source("R/functions/PlotFunctions.R")
  source("R/functions/PlotFunctions.auxiliary.R")
  source("R/functions/EegFunctions.R")
  source("R/functions/ModelSelectionFunctions.R")
  source("R/functions/FitFunctions.R")
  source("R/functions/SummarizeModels.R")
  source("R/functions/GetTFs.R")
  source("R/functions/FitFunctions_Downstream.R")
  source("R/functions/MaraDownstream.R")
  source("R/functions/SleuthFunctions.R")
  source("R/functions/ClusteringFunctions.R")
  
  
  # Inits -------------------------------------------------------------------
  
  LoadInits()
  
  # Load RNASeq ---------------------------------------------------------------
  
  
  # # Load old TPM filtered genes
  dat.long.shift <- GetCountsMatrix(prefilt.sleuth, max.time, filter.max.time = TRUE)
  dat.long.cleaned <- GetCountsMatrix(prefilt.sleuth, max.time, filter.max.time = FALSE)

  
  # Load eeg ----------------------------------------------------------------
  
  # dat.eeg.plot <- GetSubsampledEeg(inf="Robjs/eeg_data_merged_72_and_78_hr_mice/wake.df.method.mode.Robj", max.time=78)
  # dat.eeg.plot <- GetSubsampledEeg(inf="Robjs/eeg_data_merged_72_and_78_hr_mice/wake.df.method.mode.Robj", max.time=78)
  dat.eeg.plot <- GetSubsampledEeg(inf="dummyVar", max.time=78)
  wake.df <- GetWakeCollapsed(inf="dummyVar", tstep = tstep, max.time = max.time, return.wake.df.only = TRUE)
  wake.collapsed <- CollapseWake(wake.df, tstep, filter.time = low.pass.filt.time)
  
  # Load fits ---------------------------------------------------------------
  # 
  if (!prefilt.sleuth){
    # Load old TPM filtered fits
	load("data/fits_Rcpp_bugfix.maxiter.2000.LogLBicCorrected.model.sleep_mix_mixedaf_ampfreestep.lambda.1000.dolinear.FALSE.minhl.0.33.RData")
	fits <- fits.penal.fixed; rm(fits.penal.fixed)
  } else {
    # Load new counts filtered fits
    load("data/fits_Rcpp_sleuthfilt.maxiter.3000.model.LogLBicCorrected.model.sleep_mix_mixedaf.lambda.1000.dolinear.FALSE.minhl.0.33.RData")
    fits <- fits.all.fixed; rm(fits.all.fixed)
  }
  
  # Load sleuth to filter low counts ----------------------------------------
  genes.keep <- GetGenesKeep(prefilt.sleuth, counts.cutoff, fits)

  # Load PCA matrix, dat.after
  dat.after <- LoadPCAMat(dat.long.cleaned)  # dat.after
  
  # Filter out genes ---------------------------------------------------------
  
  dat.after <- subset(dat.after, gene %in% genes.keep)
  dat.long.shift <- subset(dat.long.shift, gene %in% genes.keep)
  fits <- subset(fits, gene %in% genes.keep)
  
  # Load Cluster objs -------------------------------------------------------
  # LoadClusteringObjs.runKmeans(dat.long.shift, clstr.hash, clstr.order)  # runs Kmeans
  LoadClusteringObjs()  # load from Robj
  mat.means$cluster <- clstrs$cluster
  mat.means <- mat.means %>%
    arrange(pval)
  models.all <- hash(fits$gene, fits$model)
  mat.means.withmodel <- mat.means
  mat.means.withmodel$model <- sapply(mat.means.withmodel$gene, function(g) models.all[[g]])
  
  # assign to environment
  vars.lst <- list(dat.long.shift, fits, dat.after, dat.eeg.plot, wake.collapsed, mat.means.withmodel)
  names(vars.lst) <- c("dat.long.shift", "fits", "dat.after", "dat.eeg.plot", "wake.collapsed", "mat.means.withmodel")
  list2env(vars.lst, envir = globalenv())
  return(NULL)
}


