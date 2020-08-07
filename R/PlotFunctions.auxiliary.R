# 2017-02-27
# Jake Yeung
# Companion functions fot PlotFits() and possibly other functions

ChangeModelNames <- function(mdls){
  # Change model names
  # flat = "Flat"
  # sleep = "Sleep-wake (SW)"
  # circadian = "Sinusoidal"
  # mix = "Combined"
  # ampfree.step = "Sine with amp change (SA)"
  # mixedaf = "Combined with amp change (CA)"
  jhash <- hash(
    flat = "Flat",
    sleep = "Sleep wake (SW)",
    circadian = "Sinusoidal", 
    mix = "Combined", 
    ampfree.step = "Sine + \namp change (SA)", 
    mixedaf = "Combined + \namp change (CA)"
  )
  mdls.new <- sapply(as.character(mdls), function(mdl) jhash[[mdl]])
  return(mdls.new)
}

GetYMajorSource <- function(plotobj){
  ylabs <- ggplot_build(plotobj)$layout$panel_params[[1]]$y.major_source
  if (is.null(ylabs)){
    warning("Ylabs is NULL perhaps your hacky script (ggplot_build(m)$layout$panel_params[[1]]$y.major_source) has failed")
  }
  return(ylabs)
}

AppendTauMrna <- function(fits.sub.model, params){
  # check colnames in fits.sub.model, if contains tau.mrna, add it to params
  if ("tau.mrna" %in% colnames(fits.sub.model)){
    tau.mrna <- unlist(fits.sub.model[["tau.mrna"]])
    params <- c(params, "tau.mrna"=tau.mrna)
  }
  return(params)
}

GetParams.flat <- function(fits.sub, flat.name = "fit.flat"){
  params <- unlist(fits.sub[[flat.name]][[1]][1])  # intercept
  # add BIC to end of params
  # get BIC cname by replace "fit" with "bic" in sleep.name
  bic.name <- gsub(pattern = "fit", replacement = "bic", x = flat.name)
  params <- c(params, signif(GetBICWeights(bic.name, fits.sub), digits = 2))
  return(params)
}

GetTitle.flat <- function(params){
  params.lab <- c("intercept", "w")
  params.title <- paste0(paste(params.lab, signif(params, 2), sep = "="), collapse = ", ")
}

GetParams.sleep <- function(fits.sub, sleep.name = "fit.sleep", append.bic=TRUE){
  params <- unlist(fits.sub[[sleep.name]][[1]][1:5])
  # detect tau.mrna, add it to end of params
  params <- AppendTauMrna(fits.sub[[sleep.name]][[1]], params)
  # add BIC to end of params
  # get BIC cname by replace "fit" with "bic" in sleep.name
  if (append.bic){
    bic.name <- gsub(pattern = "fit", replacement = "bic", x = sleep.name)
    params <- c(params, signif(GetBICWeights(bic.name, fits.sub), digits = 2))
  }
  return(params)
}

GetTitle.sleep <- function(params){
  if (!"tau.mrna" %in% names(params)){
    params.lab <- c("i", "U", "tau.w", "L", "tau.s", "w")
  } else {
    params.lab <- c("i", "U", "tau.w", "L", "tau.s", "tau.mrna", "w")
  }
  params.title <- paste0(paste(params.lab, signif(params, 2), sep = "="), collapse = ", ")
}

GetParams.lowpass <- function(fits.sub){
  params <- unlist(fits.sub$fit.sleep.lowpass[[1]][1:6])
  # add BIC to end of params
  params <- c(params, signif(GetBICWeights("bic.sleep.lowpass", fits.sub), digits = 2))
  return(params)
}

GetTitle.lowpass <- function(params){
  params.lab <- c("i", "U", "tau.w", "L", "tau.s", "tau.mrna", "w")
  params.title <- paste0(paste(params.lab, signif(params, 2), sep = "="), collapse = ", ")
}

GetParams.circadian <- function(fits.sub, circ.name = "fit.circadian"){
  # add BIC to end of params
  params.sine <- unlist(fits.sub[[circ.name]][[1]][c(4, 1, 2)])
  bic.name <- gsub(pattern = "fit", replacement = "bic", x = circ.name)
  params.sine <- c(params.sine, signif(GetBICWeights(bic.name, fits.sub), digits = 2))
}

GetTitle.circadian <- function(params.sine){
  params.lab <- c("mu", "amp", "phase", "w")
  params.title <- paste0(paste(params.lab, signif(params.sine, 2), sep = "="), collapse = ", ")
  return(params.title)
}  

GetParams.mixed <- function(fits.sub, remove.na=TRUE, mix.name = "fit.mix", as.cos.sin = FALSE){
  # is it amp-phase or is it cos.part sin.part? auto-detect
  params.mix.all <- fits.sub[[mix.name]][[1]]
  params.mix <- params.mix.all[1:9]
  if (remove.na){
    # remove all parameters with NA: handles the mix.weight
    params.mix <- params.mix[which(!is.na(params.mix))]
  }
  # add tau.mrna if necessary
  if (!"tau.mrna" %in% names(params.mix)){
    params.mix <- AppendTauMrna(fits.sub[[mix.name]][[1]], params.mix)
  }
  if (as.cos.sin){
    has.ampphase <- AmpPhaseOrCosSin(names(params.mix))
    if (has.ampphase){
      amp.i <- which(names(params.mix) == "amp")
      phase.i <- which(names(params.mix) == "phase")
      cos.sin <- AmpPhaseToCosSine(params.mix[amp.i], 
                                   params.mix[phase.i])
      cos.part <- cos.sin$cos.part
      sin.part <- cos.sin$sin.part
      
      # reassignad and rename
      params.mix[amp.i] <- cos.part; names(params.mix)[amp.i] <- "cos.part"
      params.mix[phase.i] <- sin.part; names(params.mix)[phase.i] <- "phase.part"
      
    } else {
      print("Already cos and sin, no conversion...")
    }
  }
  return(params.mix)
}

GetTitle.mixed <- function(fits.sub, params.mix){
  # add new title 
  if (!"amp" %in% names(params.mix) & !"phase" %in% names(params.mix)){
    params.mix$amp = CosSineToAmpPhase(params.mix$cos.part, params.mix$sin.part, T.period = 24)$amp
    params.mix$phase = CosSineToAmpPhase(params.mix$cos.part, params.mix$sin.part, T.period = 24)$phase
  }
  params.lab.mix <- names(params.mix)
  # add BIC to end of params
  params.lab.mix <- c(params.lab.mix, "w")
  params.mix.with.w <- unlist(c(params.mix, signif(GetBICWeights("bic.mix", fits.sub), digits = 2)))
  
  # get indices
  cnames <- c("U", "tau.w", "L", "tau.s", "amp", "phase", "tau.mrna", "w")
  indx <- c()
  for (cname in cnames){
    indx.add <- which(params.lab.mix == cname)
    if (length(indx.add) == 0) warning(paste0("Unknown colname: ", cname))
    indx <- c(indx, indx.add)
  }
  
  params.title.mix <- paste0(paste(params.lab.mix[indx], 
                                   signif(params.mix.with.w, 2)[indx], sep = "="), 
                             collapse = ", ")
  return(params.title.mix)
}

GetParams.ampfree <- function(fits.sub, ampfree.cname){
  params.ampfree.all <- fits.sub[[ampfree.cname]][[1]]
  form.i <- which(colnames(params.ampfree.all) == "form")
  tswitch.i <- which(colnames(params.ampfree.all) == "tswitch")
  
  form <- params.ampfree.all[form.i]  # step/switch, decay, or others
  tswitch <- params.ampfree.all[tswitch.i][[1]]  # could be integer or "free"
  
  params.circ <- unlist(params.ampfree.all[1:3])  # intercept, amp, phase
  amp.i.begin <- 4
  amp.i.end <- form.i - 1
  params.amp <- unlist(params.ampfree.all[amp.i.begin:amp.i.end])
  
  if (is.character(params.amp)){
    # if a string, assume it is comma separated 
    params.amp <- as.numeric(strsplit(params.amp, split = ",")[[1]])
  }
  
  params.ampfree <- c(params.circ, params.amp)
  return(params.ampfree)
}

GetTitle.ampfree <- function(fits.sub, params.ampfree, ampfree.cname){
  ampfree.basename <- strsplit(ampfree.cname, split = "fit.")[[1]][[2]]
  # add new title 
  params.lab.ampfree <- names(params.ampfree)
  # add BIC to end of params
  params.lab.ampfree <- c(params.lab.ampfree, "w")
  params.lab.ampfree.with.w <- unlist(c(params.ampfree, signif(GetBICWeights(paste0("bic.", ampfree.basename), fits.sub), digits = 2)))
  params.i <- -1  # dont need to display intercept
  params.title.ampfree <- paste0(paste(params.lab.ampfree[params.i], 
                                       signif(params.lab.ampfree.with.w, 2)[params.i], sep = "="), 
                                 collapse = ", ")
  return(params.title.ampfree)
}

GetParams.mixedaf <- function(fits.sub, cname="fit.mixedaf", as.cos.sin = FALSE){
  fits.sub.model <- fits.sub[[cname]][[1]]
  # S parameters: tau.mrna should be appended here
  params.S <- fits.sub.model[1:5]
  params.S <- AppendTauMrna(fits.sub[[cname]][[1]], params.S)
  
  # get rhythmic parameters
  include.intercept <- ifelse(is.na(fits.sub.model["intercept"]), FALSE, TRUE)
  if (!include.intercept){
    params.sine <- unlist(fits.sub.model[7:8])
    amp.i.begin <- 9
  } else {
    params.sine <- unlist(fits.sub.model[6:8])
    amp.i.begin <- 9
  }
  form.i <- GetForm(fits.sub, cname, return.index=TRUE)
  amp.i.end <- form.i - 1
  params.amp <- unlist(fits.sub.model[amp.i.begin:amp.i.end])
  params.mixedaf <- c(unlist(params.S), unlist(params.sine), unlist(params.amp))
  # remove NAs from parameters: intercept NA causes bug when you put this into WeightedSFreeAmp()
  params.mixedaf <- params.mixedaf[!is.na(params.mixedaf)]
  
  if (as.cos.sin){
    has.ampphase <- AmpPhaseOrCosSin(names(params.mixedaf))
    if (has.ampphase){
      amp.i <- which(names(params.mixedaf) == "amp")
      phase.i <- which(names(params.mixedaf) == "phase")
      cos.sin <- AmpPhaseToCosSine(params.mixedaf[amp.i], 
                                   params.mixedaf[phase.i])
      cos.part <- cos.sin$cos.part
      sin.part <- cos.sin$sin.part
      
      # reassignad and rename
      params.mixedaf[amp.i] <- cos.part; names(params.mixedaf)[amp.i] <- "cos.part"
      params.mixedaf[phase.i] <- sin.part; names(params.mixedaf)[phase.i] <- "phase.part"
      
    } else {
      print("Already cos and sin, no conversion...")
    }
  }
  
  return(params.mixedaf)
}

GetTitle.mixedaf <- function(fits.sub, params.mixedaf, cname="fit.mixedaf"){
  mixedaf.basename <- strsplit(cname, split = "fit.")[[1]][[2]]
  params.lab.mixedaf <- names(params.mixedaf)
  params.lab.mixedaf <- c(params.lab.mixedaf, "w")
  params.lab.mixedaf.with.w <- unlist(c(params.mixedaf, signif(GetBICWeights(paste0("bic.", mixedaf.basename), fits.sub), digits = 2)))
  params.i <- -1  # dont need to display "init"
  params.title.mixedaf <- paste0(paste(params.lab.mixedaf[params.i], 
                                       signif(params.lab.mixedaf.with.w, 2)[params.i], sep = "="), 
                                 collapse = ", ")
  return(params.title.mixedaf)
}

GetForm <- function(fits.sub, ampfree.cname, return.index=FALSE){
  params.ampfree.all <- fits.sub[[ampfree.cname]][[1]]
  form.i <- which(colnames(params.ampfree.all) == "form")
  form <- params.ampfree.all[form.i]  # step/switch, decay, or others
  if (!return.index){
    return(form)
  } else {
    return(form.i)
  }
}



GetTswitch <- function(fits.sub, ampfree.cname){
  params.ampfree.all <- fits.sub[[ampfree.cname]][[1]]
  tswitch.i <- which(colnames(params.ampfree.all) == "tswitch")
  tswitch <- params.ampfree.all[tswitch.i][[1]]  # could be integer or "free"
  return(tswitch)
}

AmpPhaseOrCosSin <- function(param.names){
  if (any(param.names == "cos.part")){
    # has cos and sin
    has.ampphase <- FALSE
  } else if (any(param.names == "amp")){
    has.ampphase <- TRUE
  } else {
    print("Warning: params must have cos.part or amp")
    print(params.names)
  }
  return(has.ampphase)
}



# Aux ---------------------------------------------------------------------

SimpleCap <- function(x, split = "\\.") {
  s <- strsplit(x, split)[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

NDecimals <- function(x){
  # 1.0 = 1
  # 2.25 = 2
  x <- as.character(x)
  if (!grepl("\\.", x)){
    return(0)
  }
  s <- strsplit(x, split = "\\.")[[1]][[2]]
  return(nchar(s))
}

