#library(dplyr)
#library(Rcpp)


# Fit weighted circadian --------------------------------------------------


FitWeightedSCircadian <- function(dat, wake.collapsed, exprs.cname = "exprs", time.cname = "time", 
                                  condensed=FALSE, pseudo = 1, min.hl = 0.5, max.hl = 24,
                                  include.mix.weight = FALSE,
                                  include.intercept = FALSE,
                                  do.lowpass=FALSE, min.hl.mrna=NA, max.hl.mrna=NA, low.pass.filter.times=NA,
                                  jlambda = 0, jmaxit = 100){
  # fit data to weighted Process S and Circadian
  # return ML estimates, variance of estimates (from Hessian), and BIC
  # pseudo: additional constant to allow for max/min expression in upper/lower limits
  # include.mix.weight: add a mix weight to weigh sleep with circ model (useless I think)
  
  # get inits
  max.exprs <- max(dat[[exprs.cname]])
  min.exprs <- min(dat[[exprs.cname]])
  
  if (!do.lowpass){
    inits.novar.S <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8)
    lims.S <- GetLims.S(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl)
  } else {
    inits.novar.S <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8, 2)  # add mRNA lifetime
    lims.S <- GetLims.S.LowPass(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl, min.hl.mrna = min.hl.mrna, max.hl.mrna)
  }
  # inits for circadian fit by doing linear fit as initial guess
  fit.circ <- FitRhythmic(dat, T.period = 24, use.weights = FALSE, get.bic = FALSE, condensed = FALSE)  # params for circadian
  
  if (include.intercept){
    inits.novar.circ <- c(fit.circ$intercept, fit.circ$cos.part, fit.circ$sin.part)
  } else {
    inits.novar.circ <- c(fit.circ$cos.part, fit.circ$sin.part)
  }
  
  if (include.mix.weight){
    # concatenate and add weights
    w.init <- 0.5
    inits.novar <- c(inits.novar.S, inits.novar.circ, w.init)
  }  else {
    inits.novar <- c(inits.novar.S, inits.novar.circ)
  }

  # get lims
  # sleep lims independent of whether or not we include mix weights
  # circadian lims depend on whether or not we include mix weights?
  lims.circ <- GetLims.circ(Inf, -Inf, include.intercept = include.intercept, lowerpart = -Inf)  # max and min exprs only takes into effect if include.intercept is TRUE
  
  if (include.mix.weight){
    # concatenate and add weights
    w.min <- 0; w.max <- 1
    lims.lower <- c(lims.S$lower, lims.circ$lower, w.min)
    lims.upper <- c(lims.S$upper, lims.circ$upper, w.max)
  }  else {
    lims.lower <- c(lims.S$lower, lims.circ$lower)
    lims.upper <- c(lims.S$upper, lims.circ$upper)
  }
  if (any(is.na(c(inits.novar, lims.lower, lims.upper)))){
    print(c(inits.novar, lims.lower, lims.upper))
    stop("Inits, limits, cannot contain NAs")
  }
  fit.novar.optim <- tryCatch({
    fit.novar.optim <- optim(inits.novar, fn = logL.WeightedSCircadian, 
                             x = dat[[exprs.cname]], wake.collapsed = wake.collapsed, filter.times = dat[[time.cname]], 
                             has.reps = TRUE,
                             estimate.var=TRUE, 
                             include.mix.weight = include.mix.weight,
                             include.intercept = include.intercept,
                             do.lowpass=do.lowpass, low.pass.filter.times = low.pass.filter.times,
                             lambda = jlambda,
                             hessian=TRUE, method = "L-BFGS-B",
                             lower = lims.lower, upper = lims.upper,
                             control = list(maxit = jmaxit))
    # lower = rep(0, length(inits.novar)), upper = rep(Inf, length(inits.novar)))
    fit.novar.optim
  }, error = function(e) {
    # print(paste("problem with:", dat$gene[[1]]))
    fit.novar.optim <- list()
    fit.novar.optim$par <- rep(NA, length(inits.novar))
    fit.novar.optim$value <- NA
    fit.novar.optim$convergence <- paste(dat$gene[[1]], e$message)
    fit.novar.optim$hessian <- diag(NA, length(inits.novar))
    fit.novar.optim
  })
  params <- fit.novar.optim$par  # 5 params from sleep, 3 or 2 from circ, 1 from mixing (optional)
  logL <- -1 * fit.novar.optim$value  # optim minimizes -logL
  
  bic <- GetBIC.logL(logL, n = nrow(dat), k = length(params))
  
  # get param estimates: S
  init <- params[[1]]; U <- params[[2]]; tau.w <- params[[3]]; L <- params[[4]]; tau.s <- params[[5]]
  if (do.lowpass){
    tau.mrna <- params[[6]]
    circ.begin.i <- 7
  } else {
    circ.begin.i <- 6
  }
  # get param estimates: circadian
  if (include.intercept){
    intercept <- params[[circ.begin.i]]; cos.part <- params[[circ.begin.i + 1]]; sin.part <- params[[circ.begin.i + 2]]
    mix.weight.i <- 9
  } else {
    cos.part <- params[[circ.begin.i]]; sin.part <- params[[circ.begin.i + 1]]
    mix.weight.i <- 8
    intercept <- NA
  }
  
  if (include.mix.weight){
    # mixing weights
    mix.weight <- params[[mix.weight.i]]
  } else {
    mix.weight <- NA
  }
  if (!condensed){
    # get variance estimates
    params.var <- diag(fit.novar.optim$hessian)  # standard error is sqrt of this
    # var estimates: S
    init.var <- params.var[[1]]; U.var <- params.var[[2]]; tau.w.var <- params.var[[3]]; L.var <- params.var[[4]]; tau.s.var <- params.var[[5]]
    if (do.lowpass){
      tau.mrna.var <- params.var[[6]]
    }
    # var estimates: circ
    if (include.intercept){
      intercept.var <- params.var[[circ.begin.i]]; cos.part.var <- params.var[[circ.begin.i + 1]]; sin.part.var <- params.var[[circ.begin.i + 2]]
    } else {
      cos.part.var <- params.var[[circ.begin.i]]; sin.part.var <- params.var[[circ.begin.i + 1]]
      intercept.var <- NA
    }
  }
  
  if (!condensed){
    if (!do.lowpass){
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            intercept = intercept,
                            cos.part = cos.part,
                            sin.part = sin.part,
                            mix.weight = mix.weight,
                            logL = logL,
                            init.var = init.var,
                            U.var = U.var,
                            tau.w.var = tau.w.var,
                            L.var = L.var,
                            tau.s.var = tau.s.var,
                            intercept.var = intercept.var,
                            cos.part.var = cos.part.var,
                            sin.part.var = sin.part.var,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    } else {
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            tau.mrna = tau.mrna,
                            intercept = intercept,
                            cos.part = cos.part,
                            sin.part = sin.part,
                            mix.weight = mix.weight,
                            logL = logL,
                            init.var = init.var,
                            U.var = U.var,
                            tau.w.var = tau.w.var,
                            L.var = L.var,
                            tau.s.var = tau.s.var,
                            intercept.var = intercept.var,
                            cos.part.var = cos.part.var,
                            sin.part.var = sin.part.var,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    }
  } else {
    if (!do.lowpass){
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            intercept = intercept,
                            amp = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$amp,
                            phase = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$phase,
                            mix.weight = mix.weight,
                            logL = logL,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message) 
    } else {
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            tau.mrna = tau.mrna,
                            intercept = intercept,
                            amp = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$amp,
                            phase = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$phase,
                            mix.weight = mix.weight,
                            logL = logL,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    }
  }
  return(dat.out)
}

logL.WeightedSCircadian <- function(th, x, wake.collapsed, filter.times, estimate.var=FALSE, 
                                    has.reps = TRUE, include.mix.weight=FALSE, include.intercept=FALSE,
                                    do.lowpass=FALSE, low.pass.filter.times=NA, lambda=0){
  # include.mix.weight: whether or not add additional parameter for weighing between sleep and circ (superfluous)
  # include.intercept: whether we include another intercept, but with S(t) it should suffice?
  if (!estimate.var){
    sig2 <- th[1]
    params <- th[-1]  # first parameter is variance of the signal: replicates of gene expression?
  } else {
    params <- th
  }
  mu <- WeightedSCircadian(params, wake.collapsed, unique(filter.times),  
                           include.mix.weight=include.mix.weight, include.intercept=include.intercept,
                           do.lowpass=do.lowpass, low.pass.filter.times = low.pass.filter.times)
  if (has.reps){
    # handle replicates because we fit with unique(filter.times)
    mu <- DuplicateVector(mu, filter.times)
  }  
  if (estimate.var){
    # biased estimator of variance
    sig2 <- sum((x - mu) ^ 2) / length(x)
  }
  if (length(x) != length(mu)) stop("x and mu must be equal lengths")
  loss <- -sum(dnorm(x, mean=mu, sd=sqrt(sig2), log=T))
  if (lambda != 0){
    penal.factor <- PenaltyFactor(mu, filter.times, L=2)
    loss <- PenalizeLoss(loss, lambda, penal.factor)
  }
  return(loss)
}

WeightedSCircadian <- function(params, wake.collapsed, filt.time, ampphase=FALSE, 
                               include.mix.weight=FALSE, include.intercept=FALSE,
                               do.lowpass=FALSE, low.pass.filter.times=NA){
  # Do weighted mixing of params.S and params.sine
  # do params is concatenated vector of S-process params and sine wave params with weighted param
  ## S process params
  # params[1] <- init
  # params[2] <- U
  # params[3] <- tau.w
  # params[4] <- L
  # params[5] <- tau.s
  ## Sine wave params
  # params[6] <- intercept (optional)
  # params[7] <- amp (if ampphase), cos.part if not
  # params[8] <- phase if ampphase, sin.part if not
  ## Finally add weighted parameter
  # params[9] <- mixing.weight (optional)
  # 
  # Optionally pass S(t) though lowpass filter, with appropriate mRNA levels
  if (include.mix.weight){
    w <- unlist(params[9])  # mixing weight and optional
  }
  
  if (!do.lowpass){
    params.s.end.i <- 5
  } else {
    params.s.end.i <- 6
  }
  params.circ.start.i <- params.s.end.i + 1
  if (include.intercept){
    params.circ.end.i <- params.circ.start.i + 2
  } else {
    params.circ.end.i <- params.circ.start.i + 1
  }
  params.S <- unlist(params[1:params.s.end.i])
  params.sine <- unlist(params[params.circ.start.i:params.circ.end.i])
  
  # if (!do.lowpass){
  #   params.S <- unlist(params[1:5])
  #   if (include.intercept){
  #     params.sine <- unlist(params[6:8])
  #   }  else {
  #     # no intercept
  #     params.sine <- unlist(params[6:7])
  #   }
  # } else {
  #   params.S <- unlist(params[1:6])
  #   if (include.intercept){
  #     params.sine <- unlist(params[7:9])
  #   }  else {
  #     # no intercept
  #     params.sine <- unlist(params[7:8])
  #   }
  # }

  y.S <- unlist(S.process.collapsed(params.S, wake.collapsed, filt.time, do.lowpass, low.pass.filter.times), use.names = FALSE)
  y.sine <- unlist(CosSine(params.sine, 
                           time = filt.time, 
                           ampphase = ampphase, 
                           include.intercept = include.intercept), 
                   use.names = FALSE)
  if (include.mix.weight){
    y <- w * y.S + (1 - w) * y.sine 
  } else {
    y <- y.S + y.sine 
  }
  return(y)
}




# Fit FreeAmp with S process ----------------------------------------------


FitWeightedSFreeAmp <- function(dat, AmpFunc, wake.collapsed, exprs.cname = "exprs", time.cname = "time", 
                                T.period = 24, tswitch = 30, form = "switch", 
                                pseudo = 0, min.hl = 0.5, max.hl = 24, include.intercept = FALSE,
                                do.lowpass = FALSE, min.hl.mrna = NA, max.hl.mrna = NA, low.pass.filter.times = NA, 
                                jlambda = 0, jmaxit = 100){
  # fit data to weighted Process S and Free Amp
  # return ML estimates, variance of estimates (from Hessian), and BIC
  # pseudo: additional constant to allow for max/min expression in upper/lower limits
  # include.intercept: dont need to include generally if you use S.process as your mean
  
  if (include.intercept){
    warning("Warning: including intercept in mixed S model. Consider setting it to False and let S process define mean")
  }
  # get inits
  max.exprs <- max(dat[[exprs.cname]])
  min.exprs <- min(dat[[exprs.cname]])
  
  if (!do.lowpass){
    inits.novar.S <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8)
  } else {
    inits.novar.S <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8, 2)
  }
  # inits for circadian fit by doing linear fit as initial guess
  fit.circ <- FitRhythmic(dat, T.period = 24, use.weights = FALSE, get.bic = FALSE, condensed = FALSE)  # params for circadian
  if (include.intercept){
    inits.novar.circ <- c(fit.circ$intercept, fit.circ$cos.part, fit.circ$sin.part)
  } else {
    inits.novar.circ <- c(fit.circ$cos.part, fit.circ$sin.part)
  }
  # init amplitude functions
  
  params.lst <- ParseFormGetParams(form)
  amp.params <- params.lst$amp.params
  amp.lims.lower <- params.lst$amp.lims.lower
  amp.lims.upper <- params.lst$amp.lims.upper
  

  inits.novar <- c(inits.novar.S, inits.novar.circ, amp.params)
  
  # get lims
  # sleep lims independent of whether or not we include mix weights
  if (!do.lowpass){
    lims.S <- GetLims.S(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl)
  } else {
    lims.S <- GetLims.S.LowPass(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl, min.hl.mrna = min.hl.mrna, max.hl.mrna = max.hl.mrna)
  }
  # circadian lims depend on whether or not we include mix weights?
  lims.circ <- GetLims.circ(Inf, -Inf, include.intercept = include.intercept, lowerpart = -Inf)  # max and min exprs only takes into effect if include.intercept is TRUE
  # lims for amp parameter defined in if loop checking form
  
  lims.lower <- c(lims.S$lower, lims.circ$lower, amp.lims.lower)
  lims.upper <- c(lims.S$upper, lims.circ$upper, amp.lims.upper)
  
  # check limits are sane
  if (any(abs(lims.upper - lims.lower) < 1e-5)){
    warning("Warninig, upper and lower limits are very similar. Consider redefining your limits")
  }
  
  fit.novar.optim <- tryCatch({
    fit.novar.optim <- optim(inits.novar, fn = logL.WeightedSFreeAmp, 
                             x = dat[[exprs.cname]], wake.collapsed = wake.collapsed, filter.times = dat[[time.cname]], 
                             tswitch = tswitch,
                             form = form,
                             has.reps = TRUE,
                             include.intercept = include.intercept,
                             do.lowpass=do.lowpass,
                             low.pass.filter.times=low.pass.filter.times,
                             lambda = jlambda, 
                             hessian=TRUE, method = "L-BFGS-B",
                             lower = lims.lower, upper = lims.upper,
                             control = list(maxit = jmaxit))
    # lower = rep(0, length(inits.novar)), upper = rep(Inf, length(inits.novar)))
    fit.novar.optim
  }, error = function(e) {
    # print(paste("problem with:", dat$gene[[1]]))
    fit.novar.optim <- list()
    fit.novar.optim$par <- rep(NA, length(inits.novar))
    fit.novar.optim$value <- NA
    fit.novar.optim$convergence <- paste(dat$gene[[1]], e$message)
    fit.novar.optim$hessian <- diag(NA, length(inits.novar))
    fit.novar.optim
  })
  params <- fit.novar.optim$par  # 5 params from sleep, 3 or 2 from circ, 1 from mixing (optional)
  logL <- -1 * fit.novar.optim$value  # optim minimizes -logL
  
  bic <- GetBIC.logL(logL, n = nrow(dat), k = length(params))
  
  # get param estimates: S
  init <- params[[1]]; U <- params[[2]]; tau.w <- params[[3]]; L <- params[[4]]; tau.s <- params[[5]]
  if (do.lowpass){
    tau.mrna <- params[[6]]
    circ.begin.i <- 7
  } else {
    circ.begin.i <- 6
  }
  # get param estimates: circadian
  if (include.intercept){
    intercept <- params[[circ.begin.i]]; cos.part <- params[[circ.begin.i + 1]]; sin.part <- params[[circ.begin.i + 2]]
    amp.param.i <- circ.begin.i + 3
  } else {
    cos.part <- params[[circ.begin.i]]; sin.part <- params[[circ.begin.i + 1]]
    intercept <- NA
    amp.param.i <- circ.begin.i + 2  # after sin part
  }
  if (form == "switch"){
    amp.param <- params[[amp.param.i]]  # amp.new
  } else if (form == "decay"){
    amp.param <- params[[amp.param.i]]  # amp.decay
  } else {
    stop("Form must be switch or ecay")
  }
  
  # # get variance estimates: TODO
  # params.var <- diag(fit.novar.optim$hessian)  # standard error is sqrt of this
  # # var estimates: S
  # init.var <- params.var[[1]]; U.var <- params.var[[2]]; tau.w.var <- params.var[[3]]; L.var <- params.var[[4]]; tau.s.var <- params.var[[5]]
  # # var estimates: circ
  # if (include.intercept){
  #   intercept.var <- params.var[[6]]; cos.part.var <- params.var[[7]]; sin.part.var <- params.var[[8]]
  # } else {
  #   cos.part.var <- params.var[[6]]; sin.part.var <- params.var[[7]]
  #   intercept.var <- NA
  # }
  # # var estimates: amp parameters
  
  if (!do.lowpass){
    dat.out <- data.frame(init = init,
                          U = U,  # upperobund
                          tau.w = tau.w,  # when awake
                          L = L,  # lowerbound
                          tau.s = tau.s,  # when asleep
                          intercept = intercept,
                          amp = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$amp,
                          phase = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$phase,
                          amp.param = amp.param, 
                          form = form,
                          tswitch = tswitch,
                          logL = logL,
                          bic = bic,
                          convergence = fit.novar.optim$convergence,
                          msg = fit.novar.optim$message)
  } else {
    dat.out <- data.frame(init = init,
                          U = U,  # upperobund
                          tau.w = tau.w,  # when awake
                          L = L,  # lowerbound
                          tau.s = tau.s,  # when asleep
                          tau.mrna = tau.mrna,
                          intercept = intercept,
                          amp = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$amp,
                          phase = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$phase,
                          amp.param = amp.param, 
                          form = form,
                          tswitch = tswitch,
                          logL = logL,
                          bic = bic,
                          convergence = fit.novar.optim$convergence,
                          msg = fit.novar.optim$message)
  }

  return(dat.out)
}

logL.WeightedSFreeAmp <- function(th, x, wake.collapsed, filter.times, tswitch, form, 
                                  has.reps = TRUE, include.intercept=FALSE,
                                  do.lowpass=FALSE, low.pass.filter.times=NA,
                                  lambda = 0){
  # include.intercept: whether we include another intercept, but with S(t) it should suffice?
  mu <- WeightedSFreeAmp(th, 
                         wake.collapsed, 
                         unique(filter.times), 
                         tswitch,
                         form,
                         ampphase=FALSE,
                         include.intercept=include.intercept,
                         do.lowpass=do.lowpass, low.pass.filter.time=low.pass.filter.times)
  if (has.reps){
    # handle replicates because we fit with unique(filter.times)
    # expand vector to handle replicates
    mu <- DuplicateVector(mu, filter.times)
  }  
  # biased estimator of variance (always estimate)
  sig2 <- sum((x - mu) ^ 2) / length(x)
  if (length(x) != length(mu)) stop("x and mu must be equal lengths")
  loss <- -sum(dnorm(x, mean=mu, sd=sqrt(sig2), log=T))
  if (lambda != 0){
    penal.factor <- PenaltyFactor(mu, filter.times, L=2)
    loss <- PenalizeLoss(loss, lambda, penal.factor)
  }
  return(loss)
}

WeightedSFreeAmp <- function(params, wake.collapsed, filt.time, tswitch, form, ampphase=FALSE, include.intercept=FALSE, 
                             do.lowpass=FALSE, low.pass.filter.times=NA){
  # Do weighted mixing of params.S and params.sine
  # do params is concatenated vector of S-process params and sine wave params with weighted param
  ## S process params
  # params[1] <- init
  # params[2] <- U
  # params[3] <- tau.w
  # params[4] <- L
  # params[5] <- tau.s
  ## Sine wave params
  # params[6] <- intercept (optional)
  # params[7] <- amp (if ampphase), cos.part if not
  # params[8] <- phase if ampphase, sin.part if not
  ## Free amp params
  # everything beyond last sine wave parameter (8 or 7 depending on include.intercept)
  # tswitch: integer or possibly a string to indicate time at which amplitude can "float"
  # form: form at which amplitude floats after tswitch (switch, decay)
  
  if (!do.lowpass){
    params.s.end.i <- 5
  } else {
    params.s.end.i <- 6
  }
  params.circ.start.i <- params.s.end.i + 1
  if (include.intercept){
    params.circ.end.i <- params.circ.start.i + 2
  } else {
    params.circ.end.i <- params.circ.start.i + 1
  }
  params.amp.start.i <- params.circ.end.i + 1
  
  params.S <- unlist(params[1:params.s.end.i])
  params.sine <- unlist(params[params.circ.start.i:params.circ.end.i])
  params.amp <- unlist(params[params.amp.start.i:length(params)])
 
  y.S <- unlist(S.process.collapsed(params.S, wake.collapsed, filt.time, do.lowpass, low.pass.filter.times), use.names = FALSE)
  
  y.sine <- unlist(Rhythmic.FreeAmp(th = c(params.sine, params.amp), 
                             times = filt.time, 
                             AmpFunc = AmpFunc, 
                             T.period = 24, 
                             tswitch = tswitch, 
                             form = form, 
                             include.intercept = include.intercept, 
                             amp.phase = ampphase),
                   use.names=FALSE)
  y <- y.S + y.sine 
  return(y)
}


# Fit Rhythmic with Free Amp ----------------------------------------------

FitRhythmic.FreeAmp <- function(dat, AmpFunc, exprs.cname = "exprs", time.cname = "time", 
                                T.period = 24, tswitch = 27, form = "switch", include.intercept=TRUE,
                                jmaxit = 100){
  # Fit rhythmic mode but with a general parametric amplitude function
  # use lapply to vectorize this code on long dat
  # 
  # Always estimate variance from the fit
  # 
  # include.intercept: whether we fit the mean or not, perhaps useful if we include sleep?
  
  # Init parameters with standard sinusoid fit
  fit.circ <- FitRhythmic(dat, T.period = T.period, use.weights = FALSE, get.bic = FALSE, condensed = FALSE)  # params for circadian
  
  if (include.intercept){
    inits.novar.circ <- c(fit.circ$intercept, fit.circ$cos.part, fit.circ$sin.part)
  } else {
    inits.novar.circ <- c(fit.circ$cos.part, fit.circ$sin.part)
  }
  # init amplitude functions

  params.lst <- ParseFormGetParams(form)
  amp.params <- params.lst$amp.params
  amp.lims.lower <- params.lst$amp.lims.lower
  amp.lims.upper <- params.lst$amp.lims.upper
  
  # if (form == "switch"){
  #   A.new <- 0.5  # switch to half amplitude at tswitch
  #   amp.params <- A.new
  #   amp.lims.lower <- 0
  #   amp.lims.upper <- 3
  # } else if (form == "decay"){
  #   tau <- 12  # 12 hours after tswitch, amplitude goes to 1/e (~0.367)
  #   amp.params <- tau
  #   amp.lims.lower <- 1/3 / log(2)  # 1/3h half life converted to natural scale
  #   amp.lims.upper <- 24 / log(2)  # 24h half life converted to natural scale
  # } else if (form == "stepfree"){
  #   # like switch, but the time at which the switch happens is not predefined
  #   A.new <- 0.5
  #   tswitch.param <- 35  # at ZT 30
  #   amp.params <- c(A.new, tswitch.param)
  #   # time switch must happen somehwere in 2nd day
  #   amp.lims.lower <- c(0, 24)
  #   amp.lims.upper <- c(3, 48)
  # } else if (form == "decayfree"){
  #   tau <- 12  # 12 hours after tswitch, amplitude goes to 1/e (~0.367)
  #   tswitch.param <- 35  # time at which decay happens
  #   amp.min <- 0.5  # decays to a minimum?
  #   amp.params <- c(tau, tswitch.param, amp.min)
  #   amp.lims.lower <- c(1/3 / log(2), 24, 0)  # 1/3h half life converted to natural scale
  #   amp.lims.upper <- c(24 / log(2), 48, 3)  # 24h half life converted to natural scale
  # } else if (form == "decaymin"){
  #   # similar to "decay" but include a minimum (does not decay to 0)
  #   tau <- 9
  #   amp.min <- 0.5
  #   amp.params <- c(tau, amp.min)
  #   amp.lims.lower <- c(1/3 / log(2), 0)  # 1/3h half life converted to natural scale
  #   amp.lims.upper <- c(24 / log(2), 3)  # 24h half life converted to natural scale
  # } else {
  #   stop(paste("form:", form, "not yet coded"))
  # }
  inits.novar <- c(inits.novar.circ, amp.params)
  
  lims.circ <- GetLims.circ(Inf, -Inf, include.intercept = include.intercept, lowerpart = -Inf)
  lims.lower <- c(lims.circ$lower, amp.lims.lower)
  lims.upper <- c(lims.circ$upper, amp.lims.upper)
   
  fit.novar.optim <- tryCatch({
    fit.novar.optim <- optim(inits.novar, fn = logL.Rhythmic.FreeAmp, 
                             x = dat[[exprs.cname]], times = dat[[time.cname]], 
                             AmpFunc = AmpFunc,
                             T.period = T.period, 
                             tswitch = tswitch, 
                             form = form, 
                             include.intercept = include.intercept,
                             hessian=TRUE, method = "L-BFGS-B",
                             lower = lims.lower, upper = lims.upper,
                             control = list(maxit = jmaxit))
    # lower = rep(0, length(inits.novar)), upper = rep(Inf, length(inits.novar)))
    fit.novar.optim
    # # diagnostics
    # x <- dat[[exprs.cname]]
    # times <- dat[[time.cname]]
    # xhat <- Rhythmic.FreeAmp(fit.novar.optim$par, times, AmpFunc, T.period = 24, tswitch = tswitch, form = form, amp.phase=FALSE, include.intercept=TRUE)
    # # check final plot
    # plot(times, x); lines(times, xhat)
    # plot(times, xhat); lines(times, x)
    # # check inits
    # plot(times, x)
    # lines(lines(times, Rhythmic.FreeAmp(inits.novar, x, times, AmpFunc, T.period, tswitch, form, include.intercept=TRUE), type = "o"))
    # logL.Rhythmic.FreeAmp(inits.novar, x, times, AmpFunc, T.period, tswitch, form, include.intercept=TRUE)
  }, error = function(e) {
    # print(paste("problem with:", dat$gene[[1]]))
    fit.novar.optim <- list()
    fit.novar.optim$par <- rep(NA, length(inits.novar))
    fit.novar.optim$value <- NA
    fit.novar.optim$convergence <- paste(dat$gene[[1]], e$message)
    fit.novar.optim$hessian <- diag(NA, length(inits.novar))
    fit.novar.optim
  })
  params <- fit.novar.optim$par  # 5 params from sleep, 3 or 2 from circ, 1 from mixing (optional)
  logL <- -1 * fit.novar.optim$value  # optim minimizes -logL
  
  bic <- GetBIC.logL(logL, n = nrow(dat), k = length(params))
  
  # get param estimates: circadian
  if (include.intercept){
    intercept <- params[[1]]; cos.part <- params[[2]]; sin.part <- params[[3]]
    amp.param.i <- 4  # after sin part
  } else {
    cos.part <- params[[1]]; sin.part <- params[[2]]
    intercept <- NA
    amp.param.i <- 3  # after sin part
  }
  
  if (form == "switch"){
    amp.param <- params[[amp.param.i]]  # amp.new
  } else if (form == "decay"){
    amp.param <- params[[amp.param.i]]  # amp.decay
  } else if (form == "stepfree"){
    nparams <- 2
    indx.add <- nparams - 1
    amp.param.i.end <- amp.param.i + indx.add
    amp.param <- params[amp.param.i:amp.param.i.end]
  } else if (form == "decayfree"){
    nparams <- 3
    indx.add <- nparams - 1
    amp.param.i.end <- amp.param.i + indx.add
    amp.param <- params[amp.param.i:amp.param.i.end]
  } else if (form == "decaymin"){
    nparams <- 2
    indx.add <- nparams - 1
    amp.param.i.end <- amp.param.i + indx.add
    amp.param <- params[amp.param.i:amp.param.i.end]
  } else {
    stop(paste("Form not yet coded:", form))
  }
  if (length(amp.param) > 1){
    # concatenate by comma
    amp.param <- as.character(paste0(signif(amp.param, 3), collapse=","))
  }
  # # get variance estimates: TODO
  # params.var <- diag(fit.novar.optim$hessian)  # standard error is sqrt of this
  # # var estimates: circ
  # if (include.intercept){
  #   intercept.var <- params.var[[1]]; cos.part.var <- params.var[[2]]; sin.part.var <- params.var[[3]]
  # } else {
  #   cos.part.var <- params.var[[1]]; sin.part.var <- params.var[[2]]
  #   intercept.var <- NA
  # }
  # # var estimates: amplitude parameter
  # # TODO
  dat.out <- data.frame(intercept = intercept,
                        amp = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$amp,
                        phase = CosSineToAmpPhase(cos.part, sin.part, T.period = 24)$phase,
                        amp.param = amp.param,
                        form = form,
                        tswitch = tswitch,
                        logL = logL,
                        bic = bic,
                        convergence = fit.novar.optim$convergence,
                        msg = fit.novar.optim$message,
                        stringsAsFactors = FALSE)
  return(dat.out)
}

logL.Rhythmic.FreeAmp <- function(th, x, times, AmpFunc, T.period, tswitch, form, include.intercept=TRUE){
  # log likelihood for free rhythmic amplitude
  # 
  # thetas are variable depending on AmpFunc, but the base
  # Parameters depend on include.intercept
  #
  # include.intercept=TRUE:
  # th[[1]] <- mu
  # th[[2]] <- cospart
  # th[[3]] <- sinpart
  # th[4 to end] <- parameters to put into AmpFunc
  # 
  # include.intercept=FALSE (mu=0):
  # th[[1]] <- cospart
  # th[[2]] <- sinpart
  # th[3 to end] <- parameters to put into AmpFunc
  #
  
  # parse parameters
  if (include.intercept){
    yint <- th[[1]]
    cospart <- th[[2]]
    sinpart <- th[[3]]
    amp.params <- th[-c(1,2,3)]
  } else {
    yint <- 0
    cospart <- th[[1]]
    sinpart <- th[[2]]
    amp.params <- th[-c(1,2)]
  }
  mu <- Rhythmic.FreeAmp(th, times, AmpFunc, T.period, tswitch, form, include.intercept, amp.phase=FALSE)
  # biased estimator of variance
  sig2 <- sum((x - mu) ^ 2) / length(x)
  if (length(x) != length(mu)) stop("x and mu must be equal lengths")
  return(-sum(dnorm(x, mean=mu, sd=sqrt(sig2), log=T)))
}

Rhythmic.FreeAmp <- function(th, times, AmpFunc, T.period = 24, tswitch=27, form="switch", include.intercept = TRUE,
                             amp.phase = FALSE){
  # amp.phase: FALSE, assume cospart and sinpart
  # amp.phase: TRUE, assume amp and phase
  # parse parameters
  if (include.intercept){
    yint <- th[[1]]
    if (!amp.phase){
      cospart <- th[[2]]
      sinpart <- th[[3]]
    } else {
      amppart <- th[[2]]
      phasepart <- th[[3]]
    }
    amp.params <- th[-c(1,2,3)]
  } else {
    yint <- 0
    if (!amp.phase){
      cospart <- th[[1]]
      sinpart <- th[[2]]
    } else {
      amppart <- th[[1]]
      phasepart <- th[[2]]
    }
    amp.params <- th[-c(1,2)]
  }
  w <- 2 * pi / T.period
  
  # calculate Amplitude given AmpFunc which is function of time
  amp <- AmpFunc(times, amp.params, time.switch = tswitch, form=form)
  
  if (!amp.phase){
    mu <- yint + amp * (cospart * cos(w * times) + sinpart * sin(w * times))
  } else {
    mu <- yint + amp * (amppart * cos(w * times - w * phasepart))
  }
  return(mu)
}

AmpFunc <- function(times, amp.params, time.switch = 27, form = "switch"){
  # times: vector of times
  # amp.params: parameters, depends on "form"
  # time.switch: time at which AmpFunc takes into effect, otherwise A=1
  # form: type of amplitude function to be used. Possible forms are:
  #   "switch": A=1 if t <= t1. Otherwise A=A.new
  #   "decay": A=1 if t <= t1. Otherwise A=exp(-(t - t1) / tau)
  #   TODO other forms
  A <- rep(1, length(times))  # init
  if (is.numeric(time.switch)){
    switch.i <- which(times >= time.switch)  # index for replacing with AmpFunc
  } else if (is.character(time.switch)){
    switch.i <- NA  # define it later, depending on form
  } else {
    stop(paste("Time switch must be numeric or character:", time.switch, class(time.switch)))
  }
  if (form == "switch"){
    # amp.params need to be length 1.
    # amp.params[[1]]: A.new (new amplitude after time.switch)
    if (length(amp.params) != 1){
      stop(paste("Expected one parameter for 'form=switch'", paste0(amp.params, ",")))
    }
    A.new <- amp.params[[1]]
    A[switch.i] <- A.new
  } else if (form == "decay"){
    if (length(amp.params) != 1){
      stop(paste("Expected one parameter for 'form=decay'", paste0(amp.params, ",")))
    }
    tau <- amp.params[[1]]
    Avec.new <- exp(-(times[switch.i] - time.switch) / tau)
    A[switch.i] <- Avec.new
  } else if (form == "stepfree"){
    # two amp parameters
    if (length(amp.params) != 2){
      stop(paste("Expected one parameter for 'form=stepfree'", paste0(amp.params, collapse=",")))
    }
    A.new <- amp.params[[1]]
    # need to define switch.i
    time.switch <- amp.params[[2]]
    switch.i <- which(times >= time.switch)
    A[switch.i] <- A.new 
  } else if (form == "decayfree"){
    # three amp parameters
    if (length(amp.params) != 3){
      stop(paste("Expected one parameter for 'form=decayfree'", paste0(amp.params, collapse=",")))
    }
    tau <- amp.params[[1]]
    # need to define switch.i
    time.switch <- amp.params[[2]]
    switch.i <- which(times >= time.switch)
    A.min <- amp.params[[3]]
    Avec.new <- A.min + (1 - A.min) * exp(-(times[switch.i] - time.switch) / tau)
    A[switch.i] <- Avec.new
  } else if (form == "decaymin"){
    tau <- amp.params[[1]]
    A.min <- amp.params[[2]]
    Avec.new <- A.min + (1 - A.min) * exp(-(times[switch.i] - time.switch) / tau)
    A[switch.i] <- Avec.new
  } else {
    stop(paste("Form not yet coded:", form))
  }
  return(A)
}


# Fit process S  ----------------------------------------------------------

FitProcessS <- function(dat, wake.collapsed, exprs.cname = "exprs", time.cname = "time", 
                        condensed=TRUE, pseudo = 0, min.hl = 2/3, max.hl = 24, do.lowpass = FALSE, 
                        min.hl.mrna=NA, max.hl.mrna=NA, low.pass.filter.times=NA,
                        jlambda = 0, jmaxit = 100){
  # fit data to Process S
  # return ML estimates, variance of estimates (from Hessian), and BIC
  # pseudo: additional constant to allow for max/min expression in upper/lower limits
  # do.lowpass: FALSE does not consider exponential smoothing of s(t). TRUE adds effect of mRNA HL and timestep to estimate exponential smoothing.
  # jlambda: allows penalization for fits that predict differences between ZT0 and ZT24. This lambda could be fit through cross-validation? 
  max.exprs <- max(dat[[exprs.cname]])
  min.exprs <- min(dat[[exprs.cname]])
  
  if (!do.lowpass){
    n.params <- 5
    inits.novar <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8)  # 5 total parameters
    lims.S <- GetLims.S(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl)  # No MRNA half-life
  } else {
    n.params <- 6
    inits.novar <- c(dat[[exprs.cname]][1], max.exprs, 10, min.exprs, 8, 2)  # 6 total parameters
    lims.S <- GetLims.S.LowPass(max.exprs, min.exprs, pseudo = pseudo, min.hl = min.hl, max.hl = max.hl, min.hl.mrna = min.hl.mrna, max.hl.mrna = max.hl.mrna)  # HL half-life in hrs
  }
  
  lims.lower <- lims.S$lower
  lims.upper <- lims.S$upper
  
  fit.novar.optim <- tryCatch({
    fit.novar.optim <- optim(inits.novar, fn = logL.Scollapsed, 
                             x = dat[[exprs.cname]], wake.collapsed = wake.collapsed, filter.times = dat[[time.cname]], do.lowpass=do.lowpass, low.pass.filter.times=low.pass.filter.times,
                             has.reps=TRUE, 
                             lambda = jlambda,
                             estimate.var=TRUE, hessian=TRUE, method = "L-BFGS-B",
                             lower = lims.lower, upper = lims.upper,
                             control = list(maxit = jmaxit))
    # lower = rep(0, length(inits.novar)), upper = rep(Inf, length(inits.novar)))
    fit.novar.optim
  }, error = function(e) {
    # print(paste("problem with:", dat$gene[[1]]))
    fit.novar.optim <- list()
    fit.novar.optim$par <- rep(NA, n.params)
    fit.novar.optim$value <- NA
    fit.novar.optim$convergence <- paste(dat$gene[[1]], e$message)
    fit.novar.optim$hessian <- diag(NA, n.params)
    fit.novar.optim
  })
  params <- fit.novar.optim$par
  logL <- -1 * fit.novar.optim$value  # optim minimizes -logL
  
  bic <- GetBIC.logL(logL, n = nrow(dat), k = length(params))
  
  # get param estimates
  init <- params[[1]]; U <- params[[2]]; tau.w <- params[[3]]; L <- params[[4]]; tau.s <- params[[5]]
  if (do.lowpass){
    tau.mrna <- params[[6]]
  }
  if (!condensed){
    # get variance estimates
    params.var <- diag(fit.novar.optim$hessian)  # standard error is sqrt of this
    init.var <- params.var[[1]]; U.var <- params.var[[2]]; tau.w.var <- params.var[[3]]; L.var <- params.var[[4]]; tau.s.var <- params.var[[5]]; 
    if (do.lowpass){
      tau.mrna.var <- params.var[[6]]
    }
  }
  if (!condensed){
    if (!do.lowpass){
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            logL = logL,
                            init.var = init.var,
                            U.var = U.var,
                            tau.w.var = tau.w.var,
                            L.var = L.var,
                            tau.s.var = tau.s.var,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    } else {
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            tau.mrna = tau.mrna,
                            logL = logL,
                            init.var = init.var,
                            U.var = U.var,
                            tau.w.var = tau.w.var,
                            L.var = L.var,
                            tau.s.var = tau.s.var,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    }
  } else {
    if (!do.lowpass){
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            logL = logL,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    } else {
      dat.out <- data.frame(init = init,
                            U = U,  # upperobund
                            tau.w = tau.w,  # when awake
                            L = L,  # lowerbound
                            tau.s = tau.s,  # when asleep
                            tau.mrna = tau.mrna,  # delay from transfering layers
                            logL = logL,
                            bic = bic,
                            convergence = fit.novar.optim$convergence,
                            msg = fit.novar.optim$message)
    }
  }
  return(dat.out)
}

logL.Scollapsed <- function(th, x, wake.collapsed, filter.times, do.lowpass=FALSE, low.pass.filter.times=NA, 
                            estimate.var=FALSE, has.reps=FALSE, lambda = 0){
  # Can optioally do lowpass or not. The size of parameters changes depending on lowpass (6 parameters if lowpass, otherwise 5)
  # here we minimize the negative log likelihood 
  if (!estimate.var){
    sig2 <- th[1]
    params <- th[-1]  # first parameter is variance of the signal: replicates of gene expression?
  } else {
    params <- th
  }
  mu <- S.process.collapsed(params, wake.collapsed, unique(filter.times), 
                            do.lowpass = do.lowpass, low.pass.filter.times = low.pass.filter.times) 
  if (has.reps){
    # handle replicates
    mu <- DuplicateVector(mu, filter.times)
  }  
  if (estimate.var){
    # biased estimator of variance
    sig2 <- sum((x - mu) ^ 2) / length(x)
  }
  if (length(x) != length(mu)) stop("x and mu must be equal lengths")
  loss = -sum(dnorm(x, mean=mu, sd=sqrt(sig2), log=T))  # minimize the loss function
  if (lambda != 0){
    # use lambda as penalty, assumes ZT0 and ZT24 exist
    penal.factor <- PenaltyFactor(mu, filter.times, L=2)
    loss <- PenalizeLoss(loss, lambda, penal.factor)
  }
  return(loss)
}

Rcpp::cppFunction('double GetNextS(NumericVector row, NumericVector th, double Sprev) {
            double S;
            if (row[0] == 1){
              S = th[1] - (th[1] - Sprev) * exp(-row[1] / th[2]);
            } else {
              S = th[3] + (Sprev - th[3]) * exp(-row[1] / th[4]);
            }
            return S;
}')

S.process.collapsed <- function(th, wake.collapsed, filter.times, do.lowpass = FALSE, low.pass.filter.times = NA){
  # th1 -> init value
  # th2 -> Upper asymptote
  # th3 -> Decay constant for wake
  # th4 -> Lower asymptote
  # th5 -> Decay constant for sleep
  # th6 -> If do.lowpass==TRUE, then degradation rate of mRNA
  # 
  # low.pass.filter.times -> If low.pass.filter.times != NA, then vector of evenly spaced times for filtering in order to filter S.out into
  # equal timesteps. Timestep will be calculated from low.pass.filter.times, so it needs to be evenly spaced.
  # low.pass.filter.times: filter.times must be subset of low.pass.filter.times
  # 

  S.prev <<- th[1]  # init
  S.out <- apply(wake.collapsed, 1, function(row){
    S <- GetNextS(row, th, S.prev)
    # wake <- row[1]
    # t.dur <- row[2]
    # if (wake == 1){
    #   S <- th[2] - (th[2] - S.prev) * exp(-t.dur / th[3])
    # } else if (wake == 0){
    #   S <- th[4] + (S.prev - th[4]) * exp(-t.dur / th[5])
    # } else {
    #   warning(paste0("wake should be 0 or 1: ", wake))
    #   S <- NA
    # }
    S.prev <<- S
    return(S)
  })
  if (do.lowpass){
    tau.mrna <- th[6]
    if (is.na(tau.mrna)){
      stop("Tau mRNA is NA")
    }
    gamma.m <- 1 / tau.mrna
    tstep <- unique(signif(diff(low.pass.filter.times), digits = 4))  # otherwise you don't get good uniques
    if (length(tstep) != 1){
      stop(paste("Time-step must be an evenly spaced numeric vector. Found:", paste(low.pass.filter.times, collapse=",")))
    }
    
    # filter S.out into evenly spaced times
    S.evenly <- S.out[which(wake.collapsed$time.cum %in% low.pass.filter.times)]
    # print("Debug")
    # plot(low.pass.filter.times, S.evenly, type = "l", main = "before low pass")
    if (length(S.evenly) != length(low.pass.filter.times)){
      stop("S.out did not get filtered properly into evenly spaced timepoints")
    }
    # low-pass filter S.out signal
    S.lowpass <- LowPass(S.evenly, tstep, gamma.m)
    # plot(low.pass.filt.time, S.lowpass, type = "l", main = "after low pass")
    # Integrate S.lowpass into vector of size |S.out| in order to
    # allow easy filtering by filter.times
    S.lowpass.long <- rep(NA, length(S.out))
    S.lowpass.long[which(wake.collapsed$time.cum %in% low.pass.filter.times)] <- S.lowpass
    S.out <- S.lowpass.long
  }
  if (!missing(filter.times)){
    S.out <- S.out[which(wake.collapsed$time.cum %in% filter.times)]
  }
  
  if (filter.times[[1]] == 0){
    # if starts at 0, first value is skipped. This is not the case if it starts at say 3 hr
    S.out <- c(th[1], S.out)  # add init value to output
  }
  return(S.out)
}


Sprocess <- function(wake.df, U = 1, L = 0, tau.i = 10 * 3600, tau.d = 10 * 3600, start = 5, tstep = 4, filter.time){
  wake <- wake.df$wake
  time <- wake.df$time.shift
  S <- rep(NA, length(wake))
  S[1] <- start
  for (i in seq(2, length(wake))){
    w <- wake[i]
    if (w == 1){
      # awake
      S[i] <- U - (U - S[(i - 1)]) * exp(-tstep / tau.i)
    } else if (w == 0){
      # sleep
      S[i] <- L + (S[(i - 1)] - L) * exp(-tstep / tau.d)
    } else {
      warning("w must be 0 or 1")
    }
  }
  if (!missing(filter.time)){
    S <- S[time %in% filter.time]
  }
  return(S)
}


# Linear models -----------------------------------------------------------

FitGeneric <- function(dat, use.weights=FALSE, get.bic=TRUE, condensed=FALSE){
  # Fit generic model with parameters equal to number of timepoints
  # return BIC
  # Expect columns exprs and time
  # Set up dat so it is sorted by time and make time a "character"
  dat <- dat %>% 
    arrange(time) %>% 
    ungroup() %>% 
    mutate(time = as.character(time))
  if (!use.weights){
    jfit <- lm(exprs ~ 0 + time, data = dat)
  } else {
    jfit <- lm(exprs ~ 0 + time, data = dat, weights = sqrt(Var))
  }
  time.vec <- as.numeric(gsub("time", "", names(coef(jfit))))
  if (!get.bic){
    dat.out <- data.frame(coef(jfit))
  }
  if (get.bic){
    n <- length(jfit$residuals)
    k <- length(coef(jfit))
    rss <- sum(jfit$residuals ^ 2)
    sig2 <- rss / n
    mu <- coef(jfit)  # use means at each timepoint
    logL <- sum(dnorm(dat$exprs, mean=mu, sd=sqrt(sig2), log=T))
    bic.logL <- GetBIC.logL(logL, n, k)
    if (!condensed){
      dat.out <- data.frame(t(coef(jfit)),
                            rss = rss,
                            logL = logL, 
                            bic = bic.logL)
    } else {
      dat.out <- data.frame(t(coef(jfit)),
                            rss = rss,
                            logL = logL, 
                            bic = bic.logL)
    }
  }
  return(dat.out)
}


FitRhythmic <- function(dat, T.period = 24, use.weights=FALSE, get.bic=FALSE, condensed=FALSE){
  # Fit rhythmic model to long dat. Expect dat to be per tissue per gene.
  # use lapply and ddply to vectorize this code.
  # dat needs columns: tissue time experiment exprs
  
  
  # Get parameters from complex model: used later
  GetParamsRhythModel <- function(myfit){
    model.params <- coef(myfit)
    return(list(intercept = model.params[1],
                a = model.params[2],
                b = model.params[3]))
  }
  
  #   tissue <- unique(dat$tissue)
  w = 2 * pi / T.period
  # Expect columns exprs and time
  if (!use.weights){
    rhyth.fit <- lm(exprs ~ 1 + cos(w * time) + sin(w * time), data = dat)
    flat.fit <- lm(exprs ~ 1, data = dat)
  } else {
    rhyth.fit <- lm(exprs ~ 1 + cos(w * time) + sin(w * time), data = dat, weights = sqrt(Var))
    flat.fit <- lm(exprs ~ 1, data = dat, weights = sqrt(Var))
  }
  compare.fit <- anova(flat.fit, rhyth.fit)
  pval <- compare.fit["Pr(>F)"][[1]][2]
  model.params <- GetParamsRhythModel(rhyth.fit)  # y = experimentarray + experimentrnaseq + a cos(wt) + b sin(wt)
  amp <- sqrt(model.params$a ^ 2 + model.params$b ^ 2)
  phase.rad <- atan2(model.params$b, model.params$a)
  phase.time <- (phase.rad / w) %% T.period
  
  #   print(model.params)
  #   print(amp); print(phase); print(pval)
  if (!get.bic){
    dat.out <- data.frame(cos.part = model.params$a, sin.part = model.params$b,
                          amp = amp, phase = phase.time, pval = pval, intercept = model.params$intercept)
  }
  if (get.bic){
    n <- length(rhyth.fit$residuals)
    k <- length(coef(rhyth.fit))
    rss <- sum(rhyth.fit$residuals ^ 2)
    sig2 <- rss / n
    mu <- CosSine(c(model.params$intercept, model.params$a, model.params$b), dat$time)
    logL <- sum(dnorm(dat$exprs, mean=mu, sd=sqrt(sig2), log=T))
    bic.logL <- GetBIC.logL(logL, n, k)
    if (!condensed){
      dat.out <- data.frame(cos.part = model.params$a, sin.part = model.params$b,
                            amp = amp, phase = phase.time, pval = pval, intercept = model.params$intercept,
                            rss = rss,
                            logL = logL,
                            bic = bic.logL)
    } else {
      dat.out <- data.frame(amp = amp, phase = phase.time, pval = pval, intercept = model.params$intercept,
                            logL = logL,
                            bic = bic.logL)
    }
  }
  return(dat.out)
}

FitFlat <- function(dat, use.weights=FALSE, get.bic=FALSE, condensed=FALSE){
  # Fit flat and return BIC
  # Expect columns exprs and time
  if (!use.weights){
    flat.fit <- lm(exprs ~ 1, data = dat)
  } else {
    flat.fit <- lm(exprs ~ 1, data = dat, weights = sqrt(Var))
  }
  intercept <- coef(flat.fit)[["(Intercept)"]]
  
  if (!get.bic){
    dat.out <- data.frame(intercept = intercept)
  }
  if (get.bic){
    n <- length(flat.fit$residuals)
    k <- length(coef(flat.fit))
    rss <- sum(flat.fit$residuals ^ 2)
    sig2 <- rss / n
    mu <- rep(intercept, length(dat$exprs))
    logL <- sum(dnorm(dat$exprs, mean=mu, sd=sqrt(sig2), log=T))
    bic.logL <- GetBIC.logL(logL, n, k)
    if (!condensed){
      dat.out <- data.frame(intercept = intercept,
                            rss = rss,
                            logL = logL,
                            bic = bic.logL)
    } else {
      dat.out <- data.frame(intercept = intercept,
                            rss = rss,
                            logL = logL,
                            bic = bic.logL)
    }
  }
  return(dat.out)
}

WeightedLm <- function(dat.sub, w = 2 * pi / 24){
  return(lm(formula = exprs ~ 1 + sin(w * time) + cos(w * time), data = dat.sub, weights = sqrt(Var)))
}

# Get limits --------------------------------------------------------------

GetLims.S <- function(max.exprs, min.exprs, pseudo = 0, min.hl = 0.5, max.hl = 24){
  # get limits for S process
  # TODO: replace lims lower and upper for fitting S process
  # min and max half lives (hl) in HOURS.
  
  # convert time to 1/2 to time to 1/e (easier math)
  min.g <- min.hl / log(2)
  max.g <- max.hl / log(2)
  
  lims.lower <- c(min.exprs, min.exprs - pseudo, min.g, min.exprs - pseudo, min.g)  # init, U, tau.w, L, tau.s, intercept, cos.part, sin.part
  lims.upper <- c(max.exprs, max.exprs + pseudo, max.g, max.exprs + pseudo, max.g)  # init, U, tau.w, L, tau.s, intercept, cos.part, sin.part, add pseudo arbitrarily to limit U and L
  return(list(lower = lims.lower, upper = lims.upper))
}

GetLims.S.LowPass <- function(max.exprs, min.exprs, pseudo = 0, min.hl = 0.5, max.hl = 24, min.hl.mrna = 1/4, max.hl.mrna = 48){
  # get limits for S process
  # min and max half lives (hl) in HOURS.
  
  # convert time to 1/2 to time to 1/e (easier math)
  # half-time of sleep-wake response to synthesis rate??
  min.g <- min.hl / log(2)
  max.g <- max.hl / log(2)
  # half-time of mRNA
  min.tau.mrna <- min.hl.mrna / log(2)
  max.tau.mrna <- max.hl.mrna / log(2)
  
  lims.lower <- c(min.exprs, min.exprs - pseudo, min.g, min.exprs - pseudo, min.g, min.tau.mrna)  # init, U, tau.w, L, tau.s, tau.mrna
  lims.upper <- c(max.exprs, max.exprs + pseudo, max.g, max.exprs + pseudo, max.g, max.tau.mrna)  # init, U, tau.w, L, tau.s, tau.mrna
  return(list(lower = lims.lower, upper = lims.upper))
}

GetLims.circ <- function(max.exprs, min.exprs, include.intercept=TRUE, lowerpart = -Inf){
  # intercept, cos.part, sin.part
  # TODO: replace lims lower and upper for fitting S process
  # if include.intercept: max and min exprs needs to be defined. otherwise not needed
  # 
  # Use lowerpart = -Inf
  lims.lower <- c(); lims.upper <- c()
  if (include.intercept){
    lims.lower <- min.exprs
    lims.upper <- max.exprs
  }
  lims.lower <- c(lims.lower, lowerpart, lowerpart)
  lims.upper <- c(lims.upper, Inf, Inf)
  return(list(lower = lims.lower, upper = lims.upper))
}


# Cos-Sine functions ------------------------------------------------------

CosSine <- Vectorize(function(params, time, w = 2 * pi / 24, ampphase=FALSE, include.intercept=TRUE){
  # TODO: should unit test this
  if (include.intercept){
    intercept <- params[1]
    amp.or.cos <- params[2]
    phase.or.cos <- params[3]
  } else {
    intercept <- 0
    amp.or.cos <- params[1]
    phase.or.cos <- params[2]
  }
  if (!ampphase){
    y <- intercept + amp.or.cos * cos(w  * time) + phase.or.cos * sin(w * time)
  } else {
    y <- intercept + amp.or.cos * cos(w * time - w * phase.or.cos)
  }
  return(y)
}, vectorize.args = "time")

CosSineToAmpPhase <- function(cos.part, sin.part, T.period = 24){
  # TODO: use this function in FitRhythmic
  w <- 2 * pi / T.period
  amp <- sqrt(cos.part ^ 2 + sin.part ^ 2)
  phase.rad <- atan2(sin.part, cos.part)
  phase.time <- (phase.rad / w) %% T.period
  return(list(amp = amp, phase = phase.time))
}

AmpPhaseToCosSine <- function(amp, phase.time, T.period = 24){
  # Get cos.part and sin.part from amp and phase.time
  w <- 2 * pi / T.period
  phase.rad <- phase.time * w
  cos.part <- amp * cos(phase.rad)
  sin.part <- amp * sin(phase.rad)
  return(list(cos.part = cos.part, sin.part = sin.part))
}


# Auxiliary functions -----------------------------------------------------


PenalizeLoss <- function(loss, lambda, penal.factor){
  # penalize loss function with lambda with optional exponents
  # use lambda as penalty, assumes ZT0 and ZT24 exist
  return(loss + lambda * penal.factor)
}

PenaltyFactor <- function(mu, filt.time, L = 2){
  # penalty factor: initially coded as the difference between ZT0 and ZT24 (to ensure homeostatic)
  if (all(24 %in% filt.time, 0 %in% filt.time)){
    # multiple indices, but should be all identical
    mu.0 <- unique(mu[which(filt.time == 0)]) 
    mu.24 <- unique(mu[which(filt.time == 24)])
    penal.factor <- abs(mu.0 - mu.24) ^ L  # larger difference, the more you penalize
  } else {
    warning("Expect filt.time to contain 0 and 24")
    print(unique(filt.time))
  }
  return(penal.factor)
}

DoTtest <- function(dat.sub, formula.str = "counts.norm ~trtmnt"){
  # test.out <- t.test(counts.norm ~ trtmnt, dat.sub)
  test.out <- tryCatch({
    test.out <- t.test(as.formula(formula.str), dat.sub)
    # return(test.out)
  }, error = function(e) {
    # output in form that conforms with test.out if no error
    # print("Error:")
    # print(e)
    test.out <- list(p.value = NA, estimate = c(NA, NA))  # estimate is list of 2 because we do diff() to get NA
    # return(test.out)
  })
  pval <- unname(test.out$p.value)
  sd.minus.nsd <- unname(diff(test.out$estimate))
  return(data.frame(pval=pval, sd.minus.nsd = sd.minus.nsd))
}

DuplicateVector <- function(x, p){
  # Add duplicates based on number of counts of p
  counts <- table(p)
  if (length(x) != length(counts)) stop("Length of x must equal number of unique p")
  x.long <- rep(NA, sum(counts))
  j <- 1
  for (i in seq(length(x))){
    x.long[(j:(j + counts[i] - 1))] <- x[i]
    j <- j + counts[i]
  }
  return(x.long)
}


ParseFormGetParams <- function(form){
  if (form == "switch"){
    A.new <- 0.5  # switch to half amplitude at tswitch
    amp.params <- A.new
    amp.lims.lower <- 0
    amp.lims.upper <- Inf
  } else if (form == "decay"){
    tau <- 12  # 12 hours after tswitch, amplitude goes to 1/e (~0.367)
    amp.params <- tau
    amp.lims.lower <- 1/3 / log(2)  # 1/3h half life converted to natural scale
    amp.lims.upper <- 24 / log(2)  # 24h half life converted to natural scale
  } else if (form == "stepfree"){
    # like switch, but the time at which the switch happens is not predefined
    A.new <- 0.5
    tswitch.param <- 33  # at ZT 30
    amp.params <- c(A.new, tswitch.param)
    # time switch must happen somehwere in 2nd day
    amp.lims.lower <- c(0, 24)
    amp.lims.upper <- c(3, 48)
  } else if (form == "decayfree"){
    tau <- 12  # 12 hours after tswitch, amplitude goes to 1/e (~0.367)
    tswitch.param <- 35  # time at which decay happens
    amp.min <- 0.5  # decays to a minimum?
    amp.params <- c(tau, tswitch.param, amp.min)
    amp.lims.lower <- c(1/3 / log(2), 24, 0)  # 1/3h half life converted to natural scale
    amp.lims.upper <- c(24 / log(2), 48, 3)  # 24h half life converted to natural scale
  } else if (form == "decaymin"){
    # similar to "decay" but include a minimum (does not decay to 0)
    tau <- 12
    # tau <- 0.5  # for Tef
    amp.min <- 0.5
    amp.params <- c(tau, amp.min)
    amp.lims.lower <- c(1/3 / log(2), 0)  # 1/3h half life converted to natural scale
    amp.lims.upper <- c(24 / log(2), 3)  # 24h half life converted to natural scale
  } else {
    stop(paste("form:", form, "not yet coded"))
  }
  return(list(amp.params = amp.params, amp.lims.lower = amp.lims.lower, amp.lims.upper = amp.lims.upper))
}

LowPass <- function(s, dt, gamma.m){
  # https://en.wikipedia.org/wiki/Low-pass_filter#Discrete-time_realization
  # s: synthesis function
  # dt: time step
  # gamma.m: degradation rate of m
  # 
  
  # define smoothing factor
  alpha <- (dt * gamma.m) / (dt * gamma.m + 1)
  m <- rep(NA, length(s))
  # m[[1]] <- alpha * s[[1]] / gamma.m
  m[[1]] <- s[[1]]
  for (i in seq(2, length(s))){
    m[[i]] <- m[i - 1] + alpha * (s[[i]] / gamma.m - m[i - 1])
  }
  return(m)
}

# BIC utilities -----------------------------------------------------------


GetBIC <- function(RSS, n, k){
  return(n * log(RSS / n) + k * log(n))
}

GetBIC.logL <- function(logL, n, k){
  # https://en.wikipedia.org/wiki/Bayesian_information_criterion
  return(-2 * logL + k * log(n))
}

