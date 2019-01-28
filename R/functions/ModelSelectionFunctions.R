FixLogLBic.3models <- function(fits.sub, dat.sub, wake.collapsed, low.pass.filter.times){
  # run FixLogLBic on "sleep", "mix", and "mixedaf"
  # jgene <- fits.sub$gene[[1]]
  # dat.sub <- subset(dat.long, gene == jgene)
  models <- list("sleep", "mix", "mixedaf")
  for (model in models){
    fits.sub <- FixLogLBic(fits.sub, dat.sub, wake.collapsed, low.pass.filter.times, model)
  }
  return(fits.sub)
}

FixLogLBic <- function(fits.sub, dat.sub, wake.collapsed, low.pass.filter.times, model = "auto"){
  # LogL for sleep, mix, and mixedaf include lambda in logL. Remove it and then recalculate logL.
  if (model == "auto"){
    if ("model" %in% colnames(fits.sub)){
      model <- fits.sub$model
    } else {
      warning("If model is auto, colnames of fits.sub must include 'model'")
      return(NA)
    }
  }
  
  # inputs for model
  jlambda <- 0  # can set to 1000 for debugging purposes to see if you get back what you had before
  x <- unlist(dat.sub$exprs)
  filter.times <- dat.sub$time
  
  model.cname <- paste0("fit.", model)
  jparams.all <- fits.sub[[model.cname]][[1]]
  # get index of logL and beyond, indices before this are params for model
  jparams.logL.i <- which(names(jparams.all) == "logL")
  jparams.bic.i <- jparams.logL.i + 1  # the next index should be BIC
  
  
  if (model == "sleep"){
    jparams <- GetParams.sleep(fits.sub, sleep.name = model.cname, append.bic=FALSE)
    loss.nopenal <- logL.Scollapsed(th = jparams, x, 
                                    wake.collapsed, filter.times, 
                                    do.lowpass = TRUE, low.pass.filter.times, 
                                    estimate.var = TRUE, has.reps = TRUE, lambda = jlambda)
  } else if (model == "mix"){
    jparams <- GetParams.mixed(fits.sub, mix.name = model.cname, as.cos.sin = TRUE)
    loss.nopenal <- logL.WeightedSCircadian(th = jparams, x = x, wake.collapsed = wake.collapsed, include.mix.weight = FALSE, include.intercept = FALSE,
                                            filter.times = filter.times, do.lowpass = TRUE, 
                                            low.pass.filter.times = low.pass.filter.times, 
                                            estimate.var = TRUE, has.reps = TRUE, lambda = jlambda)
  } else if (model == "mixedaf"){
    jparams <- GetParams.mixedaf(fits.sub, cname = model.cname, as.cos.sin = TRUE)
    form <- GetForm(fits.sub, ampfree.cname = model.cname)
    tswitch <- GetTswitch(fits.sub, ampfree.cname = model.cname)
    loss.nopenal <- logL.WeightedSFreeAmp(th = jparams, x = x, wake.collapsed = wake.collapsed, include.intercept = FALSE, 
                                          filter.times = filter.times, tswitch = tswitch, 
                                          form = form, do.lowpass = TRUE, 
                                          low.pass.filter.times = low.pass.filter.times, has.reps = TRUE, lambda = jlambda)
  } else {
    warning(paste("Model must be sleep, mix, or mixedaf:", model))
    return(NA)
  }
  logL.nopenal <- -1 * loss.nopenal
  
  bic.nopenal <- GetBIC.logL(logL.nopenal, n = length(x), k = length(jparams))
  
  # # Debugging options
  # print(jparams)
  # print(paste("model:", model, "OldLogL:", jparams.all[jparams.logL.i], "NewLogL:", logL.nopenal))
  # print(paste("model:", model, "OldBic:", jparams.all[jparams.bic.i], "NewBic:", bic.nopenal))
  # print(paste("n:", length(x), "k:", length(jparams))) 
  
  jparams.all[jparams.logL.i] <- logL.nopenal
  jparams.all[jparams.bic.i] <- bic.nopenal
  
  fits.sub[[model.cname]] <- list(jparams.all)
  return(fits.sub)
}

MergeLists <- function(jlist, jname1, jname2){
  merged.name <- paste(jname1, jname2, sep = "_")
  glist1 <- jlist[[jname1]]
  glist2 <- jlist[[jname2]]
  if (length(glist1) == 0 | length(glist2) == 0){
    stop("Empyt gene list")
  }
  jlist[[merged.name]] <- c(glist1, glist2)
  return(jlist)
}

RemoveLists <- function(jlist, jname){
  glist <- jlist[[jname]]
  if (length(glist) == 0){
    stop("Empyt gene list")
  }
  jlist[[jname]] <- NULL
  return(jlist)
}


MergeLowPassFits <- function(inf1, inf2){
  # Low pass fits were performed on only relevant models, merge with non-sleep related models
  
  if (missing(inf1)){
    inf1 <- "/home/yeung/projects/sleep_deprivation/Robjs/combat/fits.sleep.sleepLP.circadian.flat.mix.step.maxAmpInf.tswitch.33.dolinear.FALSE.MixedAFOnlyBugFix.FixInitLowPass.Robj"
  }
  if (missing(inf2)){
    inf2 <- "/home/yeung/projects/sleep_deprivation/Robjs/combat/fits.sleep.circadian.flat.mix.step.maxAmpInf.tswitch.33.Robj"
  }
  load(inf1, verbose = T)
  fits.mixedafs <- fits
  
  load(inf2, verbose=T)
  if (identical(fits$gene, fits.mixedafs$gene)){
    # discard mixedaf in fits, replace it with fits.mixedafs
    print("Replacing lowpass models with new fits")
    fits$fit.mixedaf <- fits.mixedafs$fit.mixedaf
    fits$bic.mixedaf <- fits.mixedafs$bic.mixedaf
    
    fits$fit.sleep <- fits.mixedafs$fit.sleep
    fits$bic.sleep <- fits.mixedafs$bic.sleep
    
    fits$fit.mix <- fits.mixedafs$fit.mix
    fits$bic.mix <- fits.mixedafs$bic.mix
  }
  # Add best model
  fits$model <- apply(fits, 1, function(row) SelectBestModel(row, colnames(fits)))
  return(fits)
}


BindParameters <- function(fit.orig, modelname, fit.cname){
  fit.orig <- subset(fit.orig, model == modelname)
  fit <- bind_rows(fit.orig[[fit.cname]])
  fit$gene <- fit.orig$gene
  fit <- fit %>% arrange(bic)
  return(fit)
}


GetBICEle <- function(x){
  return(x[[1]]$bic)
}
FixCname <- function(cname){
  # rename cnames: fit.sleep_bic -> bic.sleep
  # fix only if ends with _bic
  # if (endsWith(cname, "_bic")){
  if (grepl("_bic$", cname)){
    s <- strsplit(cname, "fit.")[[1]][[2]]
    model <- strsplit(s, "_")[[1]][[1]]
    cname <- paste0("bic.", model)
    return(cname)
  } else {
    # do nothing
    return(cname)
  }
}

AddBICToDat <- function(fits){
  jvars <- colnames(fits)[which(colnames(fits) != "gene")]
  fits.bics <- fits %>%
    group_by(gene) %>%
    mutate_at(.vars = jvars,
              .funs = funs(bic = GetBICEle(.)))
  if (length(jvars) == 1){
    # handle special case here: the colname will be named "bic", change it 
    # to fit.sleep_bic
    if ("bic" %in% colnames(fits.bics)){
      new.cname <- paste(jvars, "bic", sep = "_")
      fits.bics[[new.cname]] <- fits.bics$bic
      fits.bics$bic <- NULL
    } else {
      warning("Expected bic to be in colnames")
    }
  }
  colnames(fits.bics) <- sapply(colnames(fits.bics), FixCname)
  return(fits.bics)
}




ProcessBIC <- function(fits){
  if (is.null(fits$bic.sleep)){
    fits$bic.sleep <- sapply(fits$fit.sleep, function(s) s[["bic"]])
    fits$bic.circadian <- sapply(fits$fit.circadian, function(s) s[["bic"]])
    fits$bic.flat <- sapply(fits$fit.flat, function(s) s[["bic"]])
    # fits$bic.delta <- mapply(function(sl, cir) sl - cir, fits$bic.sleep, fits$bic.circadian)
    # fits$bic.mean <- mapply(function(x, y) mean(c(x, y)), fits$bic.sleep, fits$bic.circadian)
  }
  cnames <- c("sleep", "circadian", "flat")
  fits$model <- apply(fits, 1, function(row){
    bic.s <- as.numeric(row[[5]])
    bic.c <- as.numeric(row[[6]])
    bic.f <- as.numeric(row[[7]])
    bic.vec <- c(bic.s, bic.c, bic.f)
    return(cnames[[which.min(bic.vec)]])
  })
  return(fits)
}

GetBICWeights <- function(cname, subrow, cname.prefix=NA){
  # cnames: "bic.sleep", "bic.circadian", "bic.flat"
  # automatically detect if "bic.mix" is present
  if (!is.na(cname.prefix)){
    # add suffix to cname
    cname <- paste0(cname.prefix, cname)
  }
  if (is.null(subrow[[cname]])){
    warning(paste("Unknown cname:", cname))
  }
  weight.exp <- exp(-0.5 * subrow[[cname]])
  # grep all columns beginning with "bic.", use those for weight sum
  # cnames.all <- colnames(subrow)[grepl("^bic.", colnames(subrow))]
  cnames.i <- grep("^bic.", colnames(subrow))
  if (is.null(colnames(subrow))){
    # cnames.all <- names(subrow)[grepl("^bic.", names(subrow))]
    cnames.i <- grep("^bic.", names(subrow))
  }
  weight.sum <- sum(exp(-0.5 * unlist(subrow[cnames.i])))
  return(weight.exp / weight.sum)
}

BICToLong <- function(fits){
  fits.bic <- melt(data = subset(fits, select = c(-fit.circadian, -fit.flat)), 
                   id.vars = c("gene", "model", "fit.sleep"), 
                   measure.vars = c("bic.sleep", "bic.circadian", "bic.flat"), 
                   variable.name = "bic.model", 
                   value.name = "bic")
  
  fits.bic <- fits.bic %>%
    group_by(gene) %>%
    mutate(weight = exp(-0.5 * bic) / sum(exp(-0.5 * bic)))
}

SelectBestModel <- function(fit.row, cnames){
  # Find best model by grepping "bic." and
  # returning the model that has the lowest bic
  col.i <- grepl("^bic.", cnames)
  fit.rowsub <- fit.row[col.i]
  cnames.sub <- cnames[col.i]
  best.i <- which.min(fit.rowsub)
  models.sub <- sapply(cnames.sub, function(cname) strsplit(cname, split = "bic.")[[1]][[2]])
  model.best <- models.sub[best.i]
  return(model.best)
}

FitsToLong <- function(){
  # TODO?
}
