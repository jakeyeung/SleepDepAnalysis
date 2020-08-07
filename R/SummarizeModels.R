# Jake Yeung
# Date of Creation: 2017-08-22
# File: ~/projects/sleep_deprivation/scripts/functions/SummarizeModels.R
# Functions for summarizing and plotting summary plots for each model

SummarizeParameters <- function(fits, jmodel, ...){
  jfit.cname <- paste0("fit.", jmodel)
  if (jmodel == "sleep"){
    # use UL.delta.thres = 1.5 by default?
    params.sleep <- BindParameters(fits, modelname = jmodel, fit.cname = jfit.cname)
    params.sleep$min.tau <- as.numeric(apply(params.sleep, 1, function(row) min(row[[3]], row[[5]])))
    params.sleep$max.tau <- as.numeric(apply(params.sleep, 1, function(row) max(row[[3]], row[[5]])))
    params.sleep$UL.delta <- params.sleep$U - params.sleep$L
    UL.delta.i <- which(colnames(params.sleep) == "UL.delta")
    gene.delta.i <- which(colnames(params.sleep) == "gene")
    UL.delta.thres <- list(...)$UL.delta.thres
    if (is.null(UL.delta.thres)){
      print("Setting UL.delta.thres = 1.5 by default")
      UL.delta.thres <- 1.5
    }
    params.sleep$genelab <- apply(params.sleep, 1, function(row, UL.delta.i, gene.delta.i){
      UL.delt <- as.numeric(row[[UL.delta.i]])
      ifelse(abs(UL.delt) > UL.delta.thres, yes = row[[gene.delta.i]], no = NA)
    }, UL.delta.i, gene.delta.i)
    # sort by abs(max-min)
    params.sleep <- params.sleep %>% 
      arrange(desc(abs(UL.delta)))
    return(params.sleep)
  } else if (jmodel == "circadian"){
    params.circ <- BindParameters(fits, modelname = jmodel, fit.cname = jfit.cname)
    params.circ <- params.circ %>%
      arrange(desc(amp)) 
    return(params.circ)
  } else if (jmodel == "mix"){
    params.mixed <- BindParameters(fits, modelname = jmodel, fit.cname = jfit.cname)
    tau.s.i <- which(colnames(params.mixed) == "tau.s")
    tau.w.i <- which(colnames(params.mixed) == "tau.w")
    params.mixed$min.tau <- as.numeric(apply(params.mixed, 1, function(row) min(row[[tau.s.i]], row[[tau.w.i]])))
    params.mixed$UL.delta <- params.mixed$U - params.mixed$L
    gene.delta.i.m <- which(colnames(params.mixed) == "gene")
    UL.delta.i.m <- which(colnames(params.mixed) == "UL.delta")
    params.mixed$genelab <- apply(params.mixed, 1, function(row, UL.delta.i.m, gene.delta.i.m){
      UL.delt <- as.numeric(row[[UL.delta.i.m]])
      ifelse(abs(UL.delt) > 1, yes = row[[gene.delta.i.m]], no = NA)
    }, UL.delta.i.m, gene.delta.i.m)
    # add amp and phase
    if (c("cos.part", "sin.part") %in% colnames(params.mixed)){
      params.mixed$amp = CosSineToAmpPhase(params.mixed$cos.part, params.mixed$sin.part, T.period = 24)$amp
      params.mixed$phase = CosSineToAmpPhase(params.mixed$cos.part, params.mixed$sin.part, T.period = 24)$phase
    }
    # sort by abs min max
    params.mixed <- params.mixed %>% 
      arrange(desc(abs(UL.delta)))
    return(params.mixed)
  } else if (jmodel == "ampfree.step"){
    params.circstep <- BindParameters(fits, modelname = jmodel, fit.cname = jfit.cname)
    params.circstep$amp.final <- params.circstep$amp * params.circstep$amp.param 
    m1 <- PlotAmpPhase(subset(params.circstep), ampscale = 2, constant.amp = 5)
    gene.i.m <- which(colnames(params.circstep) == "gene")
    amp.i.m <- which(colnames(params.circstep) == "amp")
    amp2.i.m <- which(colnames(params.circstep) == "amp.final")
    params.circstep$genelab <- apply(params.circstep, 1, function(row, amp.i.m, amp2.i.m, gene.i.m){
      amp1 <- as.numeric(row[[amp.i.m]])
      amp2 <- as.numeric(row[[amp2.i.m]])
      gene <- row[[gene.i.m]]
      ifelse(amp1 > 0.45 | amp2 > 0.3, yes = gene, no = NA)
    }, amp.i.m, amp2.i.m, gene.i.m)
    # sort by change in amplitude
    params.circstep <- params.circstep %>% 
      arrange(desc(amp -amp.final))
    return(params.circstep)
  } else if (jmodel == "mixedaf"){
    params.mixedaf <- BindParameters(fits, modelname = jmodel, fit.cname = "fit.mixedaf")  # TODO jfit.cname
    params.mixedaf$amp.final <- params.mixedaf$amp * params.mixedaf$amp.param 
    tau.s.i <- which(colnames(params.mixedaf) == "tau.s")
    tau.w.i <- which(colnames(params.mixedaf) == "tau.w")
    params.mixedaf$min.tau <- as.numeric(apply(params.mixedaf, 1, function(row) min(row[[tau.s.i]], row[[tau.w.i]])))
    params.mixedaf$UL.delta <- params.mixedaf$U - params.mixedaf$L
    gene.i.m <- which(colnames(params.mixedaf) == "gene")
    amp.i.m <- which(colnames(params.mixedaf) == "amp")
    amp2.i.m <- which(colnames(params.mixedaf) == "amp.final")
    params.mixedaf$genelab <- apply(params.mixedaf, 1, function(row, amp.i.m, amp2.i.m, gene.i.m){
      amp1 <- as.numeric(row[[amp.i.m]])
      amp2 <- as.numeric(row[[amp2.i.m]])
      gene <- row[[gene.i.m]]
      ifelse(amp1 > 0.3 | amp2 > 0.2, yes = gene, no = NA)
    }, amp.i.m, amp2.i.m, gene.i.m)
    return(params.mixedaf)
  } else if (jmodel == "flat"){
    params.flat <- BindParameters(fits, modelname = jmodel, fit.cname = jfit.cname)
    return(params.flat)
  } else {
    warning(paste(jmodel, "model not yet coded"))
  }
} 

