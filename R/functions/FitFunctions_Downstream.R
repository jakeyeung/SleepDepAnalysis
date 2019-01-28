source("R/functions/PlotFunctions.auxiliary.R")

SummarizeConvergenceByModel <- function(fits){
  jmodels <- colnames(fits)[grepl("^fit.", colnames(fits))]
  jtables <- lapply(jmodels, function(model){
    jparams <- bind_rows(fits[[model]])
    jparams$model <- model
    jtable <- table(jparams$convergence)
    if (nrow(jtable) > 0){
      # jtable <- cbind(jtable, "model" = model)
      jtable <- as.data.frame(t(as.matrix(jtable)))
      jtable <- cbind(jtable, model = model)
      return(jtable)
    }
  })
  jtables <- Filter(Negate(is.null), jtables)
  return(jtables)
}

LoadTwoFitsAndMerge <- function(inf.mixedaf, inf.original){
  # new  mixedafs
  if (missing(inf.mixedaf)){
    inf.mixedaf <- "/home/yeung/projects/sleep_deprivation/Robjs/combat/fits.sleep.sleepLP.circadian.flat.mix.step.maxAmpInf.tswitch.33.dolinear.FALSE.MixedAFOnlyBugFix.FixInitLowPass.Robj"
  }
  load(inf.mixedaf, verbose = T)
  fits.mixedafs <- fits
  
  # original model 
  if (missing(inf.original)){
    inf.original <- "/home/yeung/projects/sleep_deprivation/Robjs/combat/fits.sleep.circadian.flat.mix.step.maxAmpInf.tswitch.33.Robj"
  }
  load(inf.original, verbose=T)
  
  if (identical(fits$gene, fits.mixedafs$gene)){
    # discard mixedaf in fits, replace it with fits.mixedafs
    print("Replacing lowpass models with new fits")
    fits$fit.mixedaf <- fits.mixedafs$fit.mixedaf
    fits$bic.mixedaf <- fits.mixedafs$bic.mixedaf
    
    fits$fit.sleep <- fits.mixedafs$fit.sleep
    fits$bic.sleep <- fits.mixedafs$bic.sleep
    
    fits$fit.mix <- fits.mixedafs$fit.mix
    fits$bic.mix <- fits.mixedafs$bic.mix
  } else {
    stop("Fits are not identical, doing nothing...")
  }
  return(fits)
}

PredictFlat <- function(fits.sub, time.vec){
  col.i <- which(grepl("fit.sleep", colnames(fits.sub)))
  cname <- colnames(fits.sub)[col.i]
  params.flat <- GetParams.flat(fits.sub, cname)
  S.pred <- rep(params.flat[[1]], length(time.vec))
  dat.pred.flat <- data.frame(time = time.vec, exprs = S.pred, model = "Flat")
  return(dat.pred.flat)
}

PredictAnyProcessS <- function(fits.sub, time.vec, wake.collapsed, low.pass.filter.times = NA, cname.grep = "fit.sleep", jmodel = "Process S"){
  # predict either PredictSleep, PredictMix, or PredictMixedAF based on cname.grep
  if (cname.grep == "fit.sleep"){
    dat.pred <- PredictSleep(fits.sub, time.vec, wake.collapsed, low.pass.filter.times, paste0(cname.grep, "$"), jmodel)
  } else if (cname.grep == "fit.mix"){
    dat.pred <- PredictMix(fits.sub, time.vec, wake.collapsed, low.pass.filter.times, paste0(cname.grep, "$"), jmodel)
  } else if (cname.grep == "fit.mixedaf"){
    dat.pred <- PredictMixedAF(fits.sub, time.vec, wake.collapsed, low.pass.filter.times, paste0(cname.grep, "$"), jmodel)
  } else {
    warning(paste("cname.grep should be fit.sleep, fit.mix, fit.mixedaf:", cname.grep))
  }
  return(dat.pred)
}

PredictSleep <- function(fits.sub, time.vec, wake.collapsed, low.pass.filter.times = NA, cname.grep = "fit.sleep", jmodel = "Process S"){
  col.i <- which(grepl(cname.grep, colnames(fits.sub)))
  cname <- colnames(fits.sub)[col.i] 
  
  params.sleep <- GetParams.sleep(fits.sub, cname)
  params.title.sleep <- GetTitle.sleep(params.sleep)
  if ("tau.mrna" %in% names(params.sleep)){
    jdo.lowpass <- TRUE
    jlow.pass.filter.times <- low.pass.filter.times
  } else {
    jdo.lowpass <- FALSE
    jlow.pass.filter.times <- NA
  }
  S.pred <- S.process.collapsed(params.sleep, wake.collapsed, time.vec, do.lowpass = jdo.lowpass, low.pass.filter.times = low.pass.filter.times)
  dat.pred.sleep <- data.frame(time = time.vec, exprs = S.pred, model = jmodel)
  return(dat.pred.sleep)
}

PredictCircadian <- function(fits.sub, time.vec){
  col.i <- which(grepl("fit.circadian", colnames(fits.sub)))
  cname <- colnames(fits.sub)[col.i] 
  
  params.sine <- GetParams.circadian(fits.sub, cname)
  params.title.sine <- GetTitle.circadian(params.sine)
  
  # make dataframes for Process S and Circadian
  C.pred <- CosSine(params.sine, time = time.vec, ampphase=TRUE)
  dat.pred.circ <- data.frame(time = time.vec, exprs = C.pred, model = "Circadian")
  return(dat.pred.circ)
}

PredictMix <- function(fits.sub, time.vec, wake.collapsed, low.pass.filter.times = NA, cname.grep = "fit.mix$", jmodel = "S+Circ"){
  col.i <- which(grepl(cname.grep, colnames(fits.sub)))
  cname <- colnames(fits.sub)[col.i] 
  
  params.mix <- GetParams.mixed(fits.sub, remove.na=TRUE, mix.name = cname)
  if ("tau.mrna" %in% names(params.mix)){
    jdo.lowpass <- TRUE
    jlow.pass.filter.times <- low.pass.filter.times
  } else {
    jdo.lowpass <- FALSE
    jlow.pass.filter.times <- NA
  }
  params.title.mix <- GetTitle.mixed(fits.sub, params.mix)
  
  has.ampphase <- AmpPhaseOrCosSin(names(params.mix))
  include.mix.weight <- "mix.weight" %in% names(params.mix)
  include.intercept <- "intercept" %in% names(params.mix)
  
  Mix.pred <- WeightedSCircadian(params.mix, wake.collapsed, time.vec, ampphase = has.ampphase, 
                                 include.mix.weight = include.mix.weight, include.intercept = include.intercept, 
                                 do.lowpass = jdo.lowpass, low.pass.filter.times = jlow.pass.filter.times)
  # mix.basename <- jmodel
  dat.pred.mix <- data.frame(time = time.vec, exprs = Mix.pred, model = jmodel)
  return(dat.pred.mix)
}

PredictAmpFree <- function(fits.sub, time.vec, jmodel = "fit.ampfree"){
  # get all possible colnames containig ampfree, including other suffices
  # such as ampfree.step, ... 
  ampfree.cnames <- colnames(fits.sub)[grepl("fit.ampfree", colnames(fits.sub))]
  for (ampfree.cname in ampfree.cnames){
    params.ampfree <- GetParams.ampfree(fits.sub, ampfree.cname)
    params.title.ampfree <- GetTitle.ampfree(fits.sub, params.ampfree, ampfree.cname)
    
    form <- GetForm(fits.sub, ampfree.cname)
    tswitch <- GetTswitch(fits.sub, ampfree.cname)
    # make pretty
    # ampfree.basename <- strsplit(ampfree.cname, split = "fit.")[[1]][[2]]
    # ampfree.basename <- SimpleCap(ampfree.basename)
    # ampfree.basename <- "Circ Free"
    
    ampfree.pred <- Rhythmic.FreeAmp(params.ampfree, time.vec, AmpFunc, 
                                     T.period = 24, tswitch = tswitch, form = form, include.intercept = TRUE,
                                     amp.phase = TRUE)
    dat.pred.ampfree <- data.frame(time = time.vec, exprs = ampfree.pred, model = jmodel)
    return(dat.pred.ampfree)
  }
}

PredictMixedAF <- function(fits.sub, time.vec, wake.collapsed, low.pass.filter.times = NA, cname.grep = "fit.mixedaf", jmodel = "Mixed2"){
  filt.time <- time.vec  # can remove weird kinks in prediction at ZT31 and ZT32 if necessary
  
  col.i <- which(grepl("fit.mixedaf", colnames(fits.sub)))
  cname <- colnames(fits.sub)[col.i] 
  
  params.mixedaf <- GetParams.mixedaf(fits.sub, cname)
  include.intercept <- "intercept" %in% names(params.mixedaf)
  if ("tau.mrna" %in% names(params.mixedaf)){
    jdo.lowpass <- TRUE
    jlow.pass.filter.times <- low.pass.filter.times
  } else {
    jdo.lowpass <- FALSE
    jlow.pass.filter.times <- NA
  }
  params.title.mixedaf <- GetTitle.mixedaf(fits.sub, params.mixedaf, "fit.mixedaf")
  form <- GetForm(fits.sub, cname)
  tswitch <- GetTswitch(fits.sub, cname)
  # make pretty
  # mixedaf.basename <- "Mixed2"  # jmodel
  mixedaf.pred <- WeightedSFreeAmp(params.mixedaf, wake.collapsed, filt.time, tswitch, form, ampphase = TRUE, 
                                   include.intercept = include.intercept, 
                                   do.lowpass = jdo.lowpass, low.pass.filter.times = jlow.pass.filter.times)
  dat.pred.mixedaf <- data.frame(time = time.vec, exprs = mixedaf.pred, model = jmodel)
  return(dat.pred.mixedaf)
}

PredictExprs <- function(jgene, fits.sub, time.vec, wake.collapsed, low.pass.filter.times, jmodel="auto"){
  # Predict gene expression from fits from a variety of hardcoded models
  # jmodel can be supplied or automatically determined
  # hardcode model to function
  
  # define some constants
  low.pass.filter.times <- rep(0.1, length = 780) * 1:780
  
  # fits.sub <- subset(fits, gene == jgene)
  if (jmodel == "auto"){
    jmodel <- fits.sub$model
    # print(paste("Using model", jmodel))
  }
  if (jmodel == "circadian") {
    dat.pred <- PredictCircadian(fits.sub, time.vec)
  } else if (jmodel == "ampfree.step"){
    dat.pred <- PredictAmpFree(fits.sub, time.vec) 
  } else if (jmodel == "sleep"){
    dat.pred <- PredictSleep(fits.sub, time.vec, wake.collapsed, low.pass.filter.times) 
  } else if (jmodel == "mixedaf"){
    dat.pred <- PredictMixedAF(fits.sub, time.vec, wake.collapsed, low.pass.filter.times) 
  } else if (jmodel == "mix"){
    dat.pred <- PredictMix(fits.sub, time.vec, wake.collapsed, low.pass.filter.times) 
  } else if (jmodel == "flat"){
    dat.pred <- PredictFlat(fits.sub, time.vec) 
  } else {
    warning(paste("Model:", jmodel, "not yet coded"))
  }
  return(dat.pred)
}

PlotBestFit <- function(dat.sub, fits.sub, time.vec, jgene, wake.collapsed, low.pass.filter.times, dat.pred=NA, dat.eeg = NA, add.sd.period = TRUE){
  # can supply your own dat.pred or get it from prediction
  # if dat.eeg should be dataframe if you want to plot eeg data
  if (is.na(dat.pred)){
    dat.pred <- PredictExprs(jgene, fits.sub, time.vec, wake.collapsed, low.pass.filter.times, jmodel = "auto")
  }
  if (missing(jgene)){
    jgene <- dat.sub$gene[[1]]
  }
  n.models <- length(unique(as.character(dat.pred$model)))
  jxlim <- range(time.vec)
  jtitle <- paste0(jgene, "\n", dat.pred$model[[1]])
  m <- ggplot() + geom_point(data = dat.sub, aes(x = time, y = exprs), alpha = 0.5) + geom_line(data = dat.pred, aes(x = time, y = exprs, colour = model, linetype = model)) + 
    ggtitle(jtitle) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                        legend.position = c(0.925, 0.1)) + 
    scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) +  
    # scale_y_continuous(limits = jylim) + 
    geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black") + 
    # geom_rect(aes(xmin=24, xmax=30, ymin=-Inf, ymax=Inf), alpha = 0.1) +  
    AddDarkPhases(dark.start = c(12, 36, 60), alpha = 0.05) + 
    ylab("log2 mRNA Abundance")
  if (add.sd.period){
    m <- AddSDPeriod(m)
  }
  # make colors black if only one fit
  if (n.models == 1){
    m <- m + scale_color_manual(values = "black")
  }
  if ("data.frame" %in% class(dat.eeg)){
    labsize = 20
    ylabs <- GetYMajorSource(m)
    # ylabs <- ggplot_build(m)$layout$panel_params[[1]]$y.major_source  # after ggplot2 2.2.0
    # if (is.null(ylabs)){
    #   warning("Ylabs is NULL perhaps your hacky script (ggplot_build(m)$layout$panel_params[[1]]$y.major_source) has failed")
    # }
    ndecs <- max(sapply(ylabs, NDecimals), na.rm = TRUE)
    m.eeg <- PlotEeg(dat.eeg, jxlim, n.decs = ndecs) + 
      AddDarkPhases(dark.start = c(12, 36, 60), alpha = 0.05)
    jlay <- matrix(c(rep(1, 4), 2), ncol = 1)
    m <- m + xlab("") + scale_x_continuous(breaks = NULL) + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize), title=element_text(size=0.4*labsize))
    m.eeg <- m.eeg + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize))
    multiplot(m, m.eeg, layout = jlay) 
    return(NULL)  # can change this later
  }
  return(m)
}

GetResiduals <- function(dat.sub, fits.sub, jgene){
  if (missing(jgene)){
    jgene <- dat.sub$gene[[1]]
  }
  dat.means <- dat.sub %>%
    group_by(gene, time) %>%
    summarise(exprs.means = mean(exprs))
  dat.preds <- fits.sub %>%
    group_by(gene) %>%
    do(PredictExprs(jgene, ., time.vec = time.dat, wake.collapsed = wake.collapsed, jmodel = "auto"))
  if (identical(dat.preds$gene, dat.means$gene)){
    dat.merged <- bind_cols(dat.means, dat.preds) %>%
      mutate(residuals = exprs.means - exprs)
  } else {
    warning("Expected columns to be identical")  
  }
  return(dat.merged)
}

PlotResiduals <- function(dat.sub, fits.sub, time.vec, jgene){
  if (missing(jgene)){
    jgene <- dat.sub$gene[[1]]
  }
  
  dat.resids <- GetResiduals(dat.sub, fits.sub) 
  
  jxlim <- range(time.vec)
  
  m <- ggplot() + 
    geom_line(data = dat.resids, aes(x = time, y = residuals)) + 
    ggtitle(jgene) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                         legend.position = c(0.925, 0.1)) + 
    scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) +  
    # scale_y_continuous(limits = jylim) + 
    geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black") + 
    # geom_rect(aes(xmin=24, xmax=30, ymin=-Inf, ymax=Inf), alpha = 0.1) +  
    AddDarkPhases(dark.start = c(12, 36, 60), alpha = 0.05) + 
    ylab("Delta of log2 mRNA accumulation") + 
    geom_hline(yintercept = 0, linetype = "dotted")
  return(m)
}
