#library(ggplot2)
#library(ggrepel)

#source("/home/yeung/projects/sleep_deprivation/scripts/functions/PlotFunctions.auxiliary.R")

# if (add.sd.period){
#   m <- AddSDPeriod(m)
#   m <- m + geom_rect(xmin = 24, xmax = 30, ymin = -Inf, ymax = -2.2, fill = "tomato", alpha = 0.005)
# }

AddSDPeriod <- function(m, xmin = 24, xmax = 30, y.frac = 1/7, jfill = "tomato", jalpha = 0.5){
  # add rectangle in plot to signify SD period. Y height is a fraction of total height, specified by y.frac
  #  https://gist.github.com/tomhopper/9076152  # get yaxis range from m and take y.frac from bottom
  # Update 2018-07-27: https://gist.github.com/tomhopper/9076152 uses panel_params
  y.min <- -Inf
  yrange <- GetYMajorSource(m)
  # yrange <- ggplot_build(m)$layout$panel_params[[1]]$y.range
  # if (is.null(yrange)){
  #   warning("Ylabs is NULL perhaps your hacky script (ggplot_build(m)$layout$panel_params[[1]]$y.major_source) has failed")
  # }
  y.max <- diff(yrange) * y.frac + min(yrange)
  print("Adding rectangle")
  # print(yrange)
  m <- m + annotate("rect", fill = jfill, alpha = jalpha, xmin = xmin, xmax = xmax, ymin = -Inf, ymax = y.max)
  return(m)
}

OrderDecreasing <- function(dat, jfactor, jval){
  # Reorder factors by decreasing value of jval
  # used in conjunction with barplots makes life easier
  dat[[jfactor]] <- factor(dat[[jfactor]], levels = dat[[jfactor]][order(dat[[jval]], decreasing = TRUE)])
  return(dat)
}

PlotPcs <- function(dat.nsd.mat, jxlim, jylim, jcenter=TRUE, jscale=FALSE, pc1=1, pc2=2, jcex=1.25, jcex.axis = 1.5, jcex.lab = 1.5){
  dat.nsd.mat <- as.matrix(dat.nsd.mat)
  # center rows
  dat.nsd.mat.cent <- t(scale(t(dat.nsd.mat), center = TRUE, scale = FALSE))
  
  # do PCA
  dat.pca <- prcomp(dat.nsd.mat.cent, center = TRUE, scale. = FALSE)
  
  labs <- colnames(dat.nsd.mat.cent)
  phases <- as.numeric(sapply(labs, function(l) strsplit(l, "_")[[1]][[1]]))
  phases <- sapply(phases, function(p){
    if (p > 72) return(p - 72)
    if (p > 48) return(p - 48)
    if (p > 24) return(p - 24)
    return(p)
  })
  phases.rad <- phases * 2 * pi / 24
  cols <- hsv(PhaseToHsv(phases.rad, min.phase = 0, max.phase = 2 * pi), s = 1, v = 1)
  # conds <- sapply(labs, function(l) strsplit(l, "_")[[1]][[3]])
  # pchs <- sapply(conds, function(cond){
  #   if (cond == "NSD") return(".")
  #   if (cond == "SD") return("*")
  # })
  
  pc1 <- 1
  pc2 <- 2
  pc.var <- dat.pca$sdev ^ 2 / sum(dat.pca$sdev ^ 2)
  pc1.var <- round(pc.var[pc1] * 100)
  pc2.var <- round(pc.var[pc2] * 100)
  
  colnames(dat.nsd.mat.cent) <- paste0("ZT", colnames(dat.nsd.mat.cent))
  
  if (missing(jxlim)){
    jxlim <- range(dat.pca$rotation[, pc1]) * 1.25
  }
  if (missing(jylim)){
    jylim <- range(dat.pca$rotation[, pc2]) * 1.25
  }
  wordcloud::textplot(dat.pca$rotation[, pc1], dat.pca$rotation[, pc2], 
                      words = colnames(dat.nsd.mat.cent), col = cols, 
                      cex = jcex, 
                      # pch = pchs, 
                      xlab = paste0("PC", pc1, " (", pc1.var, "%)"), 
                      ylab = paste0("PC", pc2, " (", pc2.var, "%)"), 
                      xlim = jxlim, 
                      ylim = jylim,
                      cex.axis = jcex.axis, cex.lab = jcex.lab)
}

fmt_dcimals <- function(decimals=0){
  # http://stackoverflow.com/questions/10035786/ggplot2-y-axis-label-decimal-precision
  # for axis labels
  # return a function responpsible for formatting the 
  # axis labels with a given number of decimals 
  function(x) as.character(round(x,decimals))
}

PlotPolarHistogram <- function(fits.sub, countstring="Count", ymax = "auto", n.breaks = 4){
  # make bins
  fits.sub$bin <- sapply(fits.sub$phase, function(p){
    if (p > 0.5){
      phase.bin <- round(p, digits = 0)
    } else {
      phase.bin <- 24
    }
    return(phase.bin)
  })
  # summarize bins
  fits.sub.forhist <- fits.sub %>%
    group_by(bin) %>%
    summarise(Count = length(bin))
  fits.sub.forhist$bin <- fits.sub.forhist$bin - 0.5
 
  if (ymax == "auto"){
    round.nearest <- 50
    ymax <- ceiling(max(fits.sub.forhist$Count) / round.nearest) * round.nearest  # round nearest integer as ymax
  } 
  
  jbreaks <- seq(0, ymax, length.out = n.breaks)
  m2 <- ggplot(fits.sub.forhist, aes_string(x = "bin", y = countstring)) +
    geom_bar(stat = "identity", width = 1) +
    scale_x_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) +
    scale_y_continuous(breaks = jbreaks) + 
    expand_limits(y = 0) +
    expand_limits(x = c(0, 24)) +
    xlab("Phase (ZT)") +
    ylab("Count") +
    theme_bw() +
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 6)) +
    coord_polar(theta = "x") + 
    geom_hline(yintercept = jbreaks, colour = "grey50", size = 0.2, linetype = "dashed") +
    geom_vline(xintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
    theme(panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
          panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank()) + 
    scale_y_continuous(labels = fmt_dcimals(1))
  # expand_limits(x = 0)  # we want center to be 0
  return(m2)
}

PlotGeneList <- function(genes, dat.long.shift, fits, wake.collapsed, time.vec = NULL){
  if (is.null(time.vec)){
    time.vec <- unique(jsub$time)
  }
  for (jgene in genes){
    jsub <- subset(dat.long.shift, gene == jgene)
    fits.sub <- subset(fits, gene == jgene)
    params <- unlist(fits.sub$fit.sleep[[1]][1:5])
    # add BIC to end of params
    params <- c(params, signif(GetBICWeights("bic.sleep", fits.sub), digits = 2))
    params.lab <- c("i", "U", "tau.w", "L", "tau.s", "model weight")
    params.title <- paste0(paste(params.lab, signif(params, 2), sep = "="), collapse = ", ")
    params.lab2 <- c("mu", "amp", "phase", "model weight")
    params.sine <- unlist(fits.sub$fit.circadian[[1]][c(4, 1, 2)])
    # add BIC to end of params
    params.sine <- c(params.sine, signif(GetBICWeights("bic.circadian", fits.sub), digits = 2))
    params.title2 <- paste0(paste(params.lab2, signif(params.sine, 2), sep = "="), collapse = ", ")
    
    plot(jsub$time, jsub$exprs, main = paste0(jgene, "\n", params.title, "\n", params.title2), xlab="Time", ylab="Exprs", pch = "*")
    # x <- time.vec
    # y <- S.process.collapsed(params, wake.collapsed, time.vec)
    # print(length(x))
    # print(length(y))
    lines(time.vec, S.process.collapsed(params, wake.collapsed, time.vec), lty = 1, col = "blue")  # 1: solid lines
    lines(time.vec, CosSine(params.sine, time = time.vec, ampphase = TRUE), lty = 2, col = "red")  # 2: dashed lines
    abline(v = c(24, 30), lty = "dotted")
    y.legpos <- min(jsub$exprs) + abs(diff(range(jsub$exprs))) * 0.1
    x.legpos <- max(jsub$time) - abs(diff(range(jsub$time))) * 0.25
    legend(x.legpos, y.legpos, legend = c("Process S", "Circadian"), lty = c(1, 2), col = c("blue", "red"))
  }  
}

PlotFit <- function(jgene, dat.long.shift, fits, wake.collapsed, dat.eeg = NULL, time.vec = NULL, label.transcript=FALSE, color.sd.time=FALSE, labsize = 20, low.pass.filter.times=NA, add.sd.period=TRUE){
  source("/home/yeung/projects/sleep_deprivation/scripts/functions/PlotFunctions.auxiliary.R")
  # Plot genelist but use ggplot2 in order to incorporate dat.eeg data
  if (is.null(time.vec)){
    time.vec <- unique(dat.long.shift$time)
  }
  for (g in jgene){
    jsub <- subset(dat.long.shift, gene == g)
    if (nrow(jsub) == 0){
      warning("Empty dataframe")
      return(NULL)
    }
    fits.sub <- subset(fits, gene == g)

    if (!label.transcript){
      glabel <- g
    } else {
      tx <- unique(jsub$transcript)
      if (length(tx) == 1){
        glabel <- paste(g, tx, sep = ": ")
      } else {
        warning("Single transcript expected")
      }
    }
    
    # init dat.pred and jtitlte
    dat.pred <- data.frame()  # append by rbind
    jtitle <- paste0(glabel)  # append by adding \n
    
    # REFACTOR: get params.title and params for S process
    if(any(grepl("fit.sleep", colnames(fits)))){
      col.i <- which(grepl("fit.sleep", colnames(fits)))
      cname <- colnames(fits)[col.i] 
      
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
      # S.pred <- S.process.collapsed(params.sleep, wake.collapsed, time.vec, do.lowpass = FALSE, low.pass.filter.times = low.pass.filter.times)
      # print(S.pred)
      # dat.pred.sleep <- data.frame(time = time.vec, exprs = S.pred, model = "Process S")
      dat.pred.sleep <- data.frame(time = time.vec, exprs = S.pred, model = "Sleep-wake (SW)")
      dat.pred <- rbind(dat.pred, dat.pred.sleep)
      jtitle <- paste0(jtitle, "\n", params.title.sleep)
    }
    
    if(any(grepl("fit.circadian", colnames(fits)))){
      col.i <- which(grepl("fit.circadian", colnames(fits)))
      cname <- colnames(fits)[col.i] 
      
      params.sine <- GetParams.circadian(fits.sub, cname)
      params.title.sine <- GetTitle.circadian(params.sine)
      
      # make dataframes for Process S and Circadian
      C.pred <- CosSine(params.sine, time = time.vec, ampphase=TRUE)
      # dat.pred.circ <- data.frame(time = time.vec, exprs = C.pred, model = "Circadian")
      dat.pred.circ <- data.frame(time = time.vec, exprs = C.pred, model = "Sinusoidal")
      dat.pred <- rbind(dat.pred, dat.pred.circ)
      jtitle <- paste0(jtitle, "\n", params.title.sine)
    }
    
    # # Auto detect LowPass model and include it
    # if ("fit.sleep.lowpass" %in% colnames(fits)){
    #   params.Slp <- GetParams.lowpass(fits.sub)
    #   params.title.Slp <- GetTitle.lowpass(params.Slp)
    #   
    #   Slp.pred <-  S.process.collapsed(params.Slp, wake.collapsed, time.vec, 
    #                                    low.pass.tau.mrna = TRUE, low.pass.filter.times = low.pass.filter.times)
    #   mix.basename <- "S Low Pass"
    #   dat.pred.Slp <- data.frame(time = time.vec, exprs = Slp.pred, model = mix.basename)
    #   dat.pred <- rbind(dat.pred, dat.pred.Slp)
    #   jtitle <- paste0(jtitle, "\n", params.title.Slp)
    # }
    
    # Auto detect mix model and include it
    if(any(grepl("fit.mix$", colnames(fits)))){
      col.i <- which(grepl("fit.mix$", colnames(fits)))
      cname <- colnames(fits)[col.i] 
      
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
      # mix.basename <- "S + Circ"
      mix.basename <- "Combined"
      dat.pred.mix <- data.frame(time = time.vec, exprs = Mix.pred, model = mix.basename)
      dat.pred <- rbind(dat.pred, dat.pred.mix)
      jtitle <- paste0(jtitle, "\n", params.title.mix)
    }
    
    # Auto detect decaying amplitude 
    if (any(grepl("fit.ampfree", colnames(fits)))){
      # get all possible colnames containig ampfree, including other suffices
      # such as ampfree.step, ... 
      ampfree.cnames <- colnames(fits)[grepl("fit.ampfree", colnames(fits))]
      for (ampfree.cname in ampfree.cnames){
        params.ampfree <- GetParams.ampfree(fits.sub, ampfree.cname)
        params.title.ampfree <- GetTitle.ampfree(fits.sub, params.ampfree, ampfree.cname)
        
        form <- GetForm(fits.sub, ampfree.cname)
        tswitch <- GetTswitch(fits.sub, ampfree.cname)
        # make pretty
        # ampfree.basename <- strsplit(ampfree.cname, split = "fit.")[[1]][[2]]
        # ampfree.basename <- SimpleCap(ampfree.basename)
        # ampfree.basename <- "Circ Free"
        ampfree.basename <- "Sine with amp change (SA)"
       
        ampfree.pred <- Rhythmic.FreeAmp(params.ampfree, time.vec, AmpFunc, 
                                         T.period = 24, tswitch = tswitch, form = form, include.intercept = TRUE,
                                         amp.phase = TRUE)
        dat.pred.ampfree <- data.frame(time = time.vec, exprs = ampfree.pred, model = ampfree.basename)
        dat.pred <- rbind(dat.pred, dat.pred.ampfree)
        
        jtitle <- paste0(jtitle, "\n", params.title.ampfree) 
      }
    }
    
    # Auto detect S + decaying amplitude
    if(any(grepl("fit.mixedaf", colnames(fits)))){
      col.i <- which(grepl("fit.mixedaf", colnames(fits)))
      cname <- colnames(fits)[col.i] 
      
      params.mixedaf <- GetParams.mixedaf(fits.sub, cname)
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
      # mixedaf.basename <- strsplit("fit.mixed", split = "fit.")[[1]][[2]]
      # mixedaf.basename <- "Mixed2"
      mixedaf.basename <- "Combined with amp change (CA)"
      mixedaf.pred <- WeightedSFreeAmp(params.mixedaf, wake.collapsed, filt.time, tswitch, form, ampphase = TRUE, 
                                       include.intercept = include.intercept, 
                                       do.lowpass = jdo.lowpass, low.pass.filter.times = jlow.pass.filter.times)
      dat.pred.mixedaf <- data.frame(time = time.vec, exprs = mixedaf.pred, model = mixedaf.basename)
      dat.pred <- rbind(dat.pred, dat.pred.mixedaf)
      jtitle <- paste0(jtitle, "\n", params.title.mixedaf) 
    }
    
    # PLOT
    jxlim <- range(jsub$time)
    # print(time.vec)
    # jylim <- range(jsub$exprs)
    # print(jylim)
    m <- ggplot() + geom_point(data = jsub, aes(x = time, y = exprs), alpha = 0.5) + geom_line(data = dat.pred, aes(x = time, y = exprs, colour = model, linetype = model)) + 
      ggtitle(jtitle) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                           legend.position = c(0.925, 0.1)) + 
      scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) +  
      # scale_y_continuous(limits = jylim) + 
      geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black") + 
      # geom_rect(aes(xmin=24, xmax=30, ymin=-Inf, ymax=Inf), alpha = 0.1) +  
      AddDarkPhases(dark.start = c(12, 36, 60), alpha = 0.05) + 
      ylab("log2 mRNA Abundance")
    if (color.sd.time){
      m <- m + AddDarkPhases(dark.start = c(24), dark.end = c(30), alpha = 0.1)
    }
    # add period of SD in plot
    if (add.sd.period){
      m <- AddSDPeriod(m)
    }
    if (is.null(dat.eeg)){
      return(m)
    } else {
      ylabs <- GetYMajorSource(m)
      # ylabs <- ggplot_build(m)$panel$ranges[[1]]$y.major_source  # pre ggplot 2.2.0
      # ylabs <- ggplot_build(m)$layout$panel_params[[1]]$y.major_source  # after ggplot2 2.2.0
      # if (is.null(ylabs)){
      #   warning("Ylabs is NULL perhaps your hacky script (ggplot_build(m)$layout$panel_params[[1]]$y.major_source) has failed")
      # }
      ndecs <- max(sapply(ylabs, NDecimals), na.rm = TRUE)
      m.eeg <- PlotEeg(dat.eeg, jxlim, n.decs = ndecs) + 
        AddDarkPhases(dark.start = c(12, 36, 60), alpha = 0.05)
      if (color.sd.time){
        m.eeg <- m.eeg + AddDarkPhases(dark.start = c(24), dark.end = c(30), alpha = 0.1)
      }
      # matrix(c(1,2,3,3), nrow=, byrow=TRUE) 
      jlay <- matrix(c(rep(1, 4), 2), ncol = 1)
      multiplot(m + xlab("") + scale_x_continuous(breaks = NULL) + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize), title=element_text(size=0.4*labsize)), 
                m.eeg + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize)), layout = jlay) 
    }    
  }
}

PlotPoints.withEEG <- function(jsub, dat.eeg.plot, jtitle, labsize = 20, colour.pts = TRUE, has.batch = TRUE){
  if (missing(jtitle)){
    jtitle <- unique(as.character(jsub$gene))
  }
  jxlim <- range(jsub$time)
  m.gene <- PlotPoints(jsub, expand.zero = FALSE, colour.pts = colour.pts, dotsize = 3.5, has.batch = has.batch, add.sd.period = TRUE) + 
    ylab("log2 mRNA abundance") +
    AddDarkPhases(alpha = 0.1) +
    ggtitle(jtitle) + 
    theme(legend.position= "top") 
  m.eeg <- PlotEeg(dat.eeg.plot, jxlim = jxlim)
  jlay <- matrix(c(rep(1, 4), 2), ncol = 1)
  multiplot(m.gene + xlab("") + scale_x_continuous(breaks = NULL) + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize), title=element_text(size=0.4*labsize)) , 
            m.eeg + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize)), layout = jlay) 
}

PlotPoints <- function(jsub, add.mean = TRUE, expand.zero = FALSE, 
                       colour.pts = FALSE, dotsize = 3, zt24tozt0=FALSE, has.batch = FALSE,
                       shift.day.7.to.4 = FALSE, group.tx = FALSE, add.sd.period = TRUE){
  if (shift.day.7.to.4){
    # shift 192 and 198 by 4 days
    jsub$time[which(jsub$time >= 192)] <- jsub$time[which(jsub$time >= 192)] - 24 * 4
  }
  # plots points, adds line as mean
  #
  if (!exists("collapse.integers", mode = "function")){
    source("/home/yeung/projects/sleep_deprivation/scripts/functions/IntegerFunctions.R")
  }
  if (zt24tozt0){
    jsub$time <- as.numeric(gsub(24, 0, jsub$time))
  }
  decimal.places <- 1
  jxlim <- range(jsub$time)
  
  # DECIDE COLORING AND SHAPE OF POINTS
  if (!group.tx){
    if (!colour.pts){
      m <- ggplot() + geom_point(data = jsub, aes(x = time, y = exprs), size = dotsize)
    } else {
      # order jsub samps by order index 
      jsub <- jsub %>%
        group_by(gene, time) %>%
        # mutate(samp = sort.int(samp, index.return=TRUE)$ix)
        mutate(samp = collapse.integers(samp))  # IntegerFunctions.R
      if (!has.batch){
        if (colour.pts){
          m <- ggplot() + geom_point(data = jsub, aes(x = time, y = exprs, colour = as.factor(ZT)), size = dotsize) + 
            scale_colour_discrete(guide = guide_legend(title = "ZT"))
        } else {
          m <- ggplot() + geom_point(data = jsub, aes(x = time, y = exprs), size = dotsize) + 
            scale_colour_discrete(guide = guide_legend(title = "ZT"))
        }
      } else {
        if (colour.pts){
          m <- ggplot() + 
            geom_point(data = jsub, aes(x = time, y = exprs, colour = as.factor(ZT), shape = as.factor(batch)), size = dotsize) + 
            scale_colour_discrete(guide = guide_legend(title = "ZT")) + 
            scale_shape_discrete(guide = guide_legend(title = "Batch"))
        } else {
          m <- ggplot() + 
            geom_point(data = jsub, aes(x = time, y = exprs, shape = as.factor(batch)), size = dotsize) + 
            scale_colour_discrete(guide = guide_legend(title = "ZT")) + 
            scale_shape_discrete(guide = guide_legend(title = "Batch"))
        }
        }
      }
  } else {
    m <- ggplot() + 
      geom_point(data = jsub, aes(x = time, y = exprs, colour = transcript, shape = transcript), size = dotsize)
  }
  
  # DECIDE LINE TYPES 
  if (add.mean){
    if (!"transcript" %in% colnames(jsub)){
      jsub.mean <- jsub %>%
        group_by(time) %>%
        summarise(exprs = mean(exprs))
    } else {
      jsub.mean <- jsub %>%
        group_by(time, transcript) %>%
        summarise(exprs = mean(exprs))
    }
    if (!group.tx){
      m <- m + geom_line(data = jsub.mean, aes(x = time, y = exprs))
    } else {
      m <- m + geom_line(data = jsub.mean, aes(x = time, y = exprs, group = transcript, linetype = transcript))
    }
  }
  # make pretty
  jxlim <- range(jsub$time)
  m <- m + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    # scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) +
    geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black")
  if (expand.zero){
    jxlim[1] <- 0
    m <- m + scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6), expand = c(0, 0)) + 
      scale_y_continuous(labels = fmt_dcimals(1))
  } else {
    m <- m + scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6), expand = c(0, 0)) + 
      scale_y_continuous(labels = fmt_dcimals(1))
  }
  # add red rectangle to represent 
  if (add.sd.period){
    m <- AddSDPeriod(m)
  }
  return(m)
}

PlotPoints.atacseq <- function(dat.atacseq, chromostartendgene.str, exprs.cname = "counts.norm", convert.log2=TRUE, zt24tozt0=TRUE, add.sd.period = FALSE){
  # chromostartendgene.str: "chr18  32438373  32438873 Bin1" <- auto split into chromo, start, end
  # filter dat.atacseq, plot the dots
  # rename counts to exprs
  
  if (!exists("collapse.integers", mode = "function")){
    source("scripts/functions/IntegerFunctions.R")
  }
  
  jchromo <- ExtractChromoStartEndGene(chromostartendgene.str)$chromo
  jstart <- ExtractChromoStartEndGene(chromostartendgene.str)$start
  jend <- ExtractChromoStartEndGene(chromostartendgene.str)$end
  jgene <- ExtractChromoStartEndGene(chromostartendgene.str)$gene
  dat.sub <- subset(dat.atacseq, chromo == jchromo & start == jstart & end == jend & gene == jgene)
  dat.sub <- dplyr::rename(dat.sub, exprs = counts.norm)
  if (convert.log2){
    eps <- 1e-3
    dat.sub$exprs <- log2(dat.sub$exprs + eps)
    ylab. <- "Log2 Normalized Counts"
  } else {
    ylab. <- "Normalized Counts"
  }
  dat.sub$samp <- sapply(dat.sub$sample, function(s) strsplit(s, "_")[[1]][[3]])
  dat.sub <- dat.sub %>%
    group_by(gene, time) %>%
   #  mutate(samp = sort.int(samp, index.return=TRUE)$ix)
    mutate(samp = collapse.integers(samp))  # IntegerFunctions.R
  if (zt24tozt0){
    dat.sub$time <- as.numeric(gsub(24, 0, dat.sub$time))
  }
  m <- PlotPoints(dat.sub, expand.zero = TRUE, colour.pts = TRUE, dotsize = 3.5, add.sd.period = add.sd.period) + ylab(ylab.)
  return(m)
}

AddDarkPhases <- function(dark.start, dark.end, alpha = 0.05){
  # mark dark phases
  if (missing(dark.start)){
    dark.start <- c(12, 36, 60)
  }
  if (missing(dark.end)){
    dark.end <- dark.start + 12
  }
  zt.marks <- data.frame(dark.start = dark.start, dark.end = dark.end, ymin = -Inf, ymax = Inf)
  return(geom_rect(data = zt.marks, aes(x = NULL, y = NULL, xmin = dark.start, xmax = dark.end, ymin = ymin, ymax = ymax), alpha = alpha))
}

PlotEeg <- function(dat.eeg, jxlim = c(-24, 78), n.decs = 1){
  if (!"mouse" %in% colnames(dat.eeg)){
    jwidth <- (4/3600) * (86401 / nrow(dat.eeg))
  } else {
    jmouse <- dat.eeg$mouse[[1]]
    jwidth <- (4/3600) * (86401 / nrow(subset(dat.eeg, mouse == jmouse)))
  }
  m.eeg <- ggplot(dat.eeg, aes(x = time.shift, y = w.smooth)) + geom_bar(stat = "identity", width = jwidth) + geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black") + 
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab("Time [h]") + ylab("Wake") + 
    # ggtitle("Minutes awake / 5 min interval. SD between 24 and 30") + 
    ggtitle("") + 
    scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) + 
    scale_y_continuous(limits = c(0, 5.1), breaks = c(0.0, 5.0), labels = format(c(0, 5), nsmall = n.decs))
  return(m.eeg)
}

fmt_dcimals <- function(decimals=0){
  # http://stackoverflow.com/questions/10035786/ggplot2-y-axis-label-decimal-precision
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

PlotExprsVsWakefulness <- function(jsub, jtitle, xstr = "Wakefulness_12", ncol = 2){
  if (missing(jtitle)){
    jtitle <- as.character(unique(jsub$gene))
  }
  jsub.samp <- jsub %>%
    group_by(trtmnt, time) %>%
    filter(exprs == max(exprs)) %>%
    mutate(log2shift = log2exprs + 0.05)
  if (xstr == "all"){
    plotslst <- lapply(t.ints, function(t.int){
      xstr <- paste0("Wakefulness_", t.int) 
      m <- ggplot(jsub, aes_string(x = xstr, y = "log2exprs", colour = "trtmnt")) + geom_point(size = 5) + ggtitle(jgene) + 
        # geom_smooth(data = jsub, method = "lm", formula = "log2exprs ~ Wakefulness_12")) + 
        geom_text(data = jsub.samp, aes_string(x = xstr, y = paste0("log2shift"), label = "time")) + theme_bw(6) + theme(aspect.ratio = 1)
    })
    # f <- function(m) multiplot(m, cols = 2)
    do.call(multiplot, c(plotslst, cols = ncol))
    # multiplot(plotslst[[1]], plotslst[[2]], plotslst[[3]], plotslst[[4]], cols = 2)
  } else {
    ggplot(jsub, aes_string(x = xstr, y = "log2exprs")) + geom_point(size = 5, aes(colour = trtmnt)) + ggtitle(jtitle) + 
      geom_text(data = jsub.samp, aes_string(x = xstr, y = "log2shift", label = "time")) + theme_bw() + geom_smooth(method = "lm")
  }
}

PlotAmpPhase <- function(dat,
                         jtitle = "", xlab = "Amplitude", ylab = "ZT", 
                         ampscale = 2, constant.amp = FALSE, dot.col = "gray85", jsize = 22,
                         take.top = 20,
                         gene.label = c()){
  # Expect amp and phase in data frame column names.
  # label as gene names
  dat$amp <- dat$amp * ampscale
  
  amp.cutoff <- sort(dat$amp, decreasing = TRUE)[20]
  
  dat$gene.orig <- dat$gene
  # unlabel amplitudes below cutoff
  dat$gene <- mapply(function(gene, amp){
    if (amp > amp.cutoff | gene %in% gene.label){
      return(gene)
    } else {
      return("")
    }
  }, as.character(dat$gene), dat$amp)
  
  amp.max <- ceiling(max(dat$amp) * 2) / 2
  if (amp.max <= 1){
    amp.step <- 0.5
  } else {
    amp.step <- 1
  }
  
  
  if (class(dot.col) == "character" & length(dot.col) == 1){
    m <- ggplot(data = dat, aes(x = amp, y = phase, label = gene)) +
      geom_point_rast(size = 0.5, colour = dot.col)
  } else if (class(dot.col) == "hash"){
    # hash maps gene to color
    dat$jcol <- sapply(as.character(dat$gene.orig), function(g) AddFromHash(g, dot.col))  # from HashFunctions.R
    m <- ggplot(data = dat, aes(x = amp, y = phase, label = gene, color = jcol)) +
      geom_point_rast(size = 0.5) + 
      scale_color_identity()
  } else if (class(dot.col) == "character" & length(dot.col) > 1){
    # by cluster colored from 1 to 10
    m <- ggplot(data = dat, aes(x = amp, y = phase, label = gene, color = cluster)) +
      geom_point_rast(size = 0.5) + 
      scale_color_manual(values = dot.col)
  } else {
    warning("dot.col must be hash or character")
    return(NULL)
  }
  
  m <- m + 
    coord_polar(theta = "y") + 
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(jtitle) +
    scale_y_continuous(limits = c(0, 24), breaks = seq(6, 24, 6)) + 
    scale_x_continuous(limits = c(0, amp.max), breaks = seq(0, amp.max, length.out = 2)) + 
    theme_bw(jsize) + 
    geom_vline(xintercept = seq(0, amp.max, length.out = 2), colour = "grey50", size = 0.2, linetype = "dashed") +
    geom_hline(yintercept = seq(6, 24, by = 6), colour = "grey50", size = 0.2, linetype = "solid") +
    theme(
          # panel.grid.major = element_line(size = 0.5, colour = "grey"), panel.grid.minor = element_blank(), 
          # panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="bottom",
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.key = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          line = element_blank())
  
  # add text
  dat.txt <- subset(dat, gene != "")
  # print(dat)
  # print(dat.txt)
  if (constant.amp != FALSE){
    # m <- m + geom_text(data = dat.txt, aes(x = amp, y = phase), vjust = 0, size = constant.amp)
    m <- m + geom_text_repel(data = dat.txt, aes(x = amp, y = phase, label = gene), size = constant.amp)
    # m <- m + geom_text_repel(data = dat.txt, aes(x = amp, y = phase))
  } else {
    # m <- m + geom_text(data = dat.txt, aes(x = amp, y = phase, size = amp), vjust = 0)
    m <- m + geom_text_repel(data = dat.txt, aes(x = amp, y = phase, size = amp, label = gene))
  }
  return(m)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  library(grid)  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

PlotGeneAcrossTimeSamps <- function(dat.sub, time_str = "time", exprs_str = "exprs", var.colname = "Var", show.plot = TRUE){
  # Plot time course across time, assumes has multiple samples per timepoint
  dat.sub$ymax <- dat.sub[[exprs_str]] + sqrt(dat.sub[[var.colname]])
  dat.sub$ymin <- dat.sub[[exprs_str]] - sqrt(dat.sub[[var.colname]])
  
  m <- ggplot(dat.sub, aes_string(x = time_str, y = exprs_str)) + geom_point()
  # add  var
  limits <- aes(ymax = ymax, ymin = ymin)
  m <- m + geom_errorbar(limits, width=0.25)
  if (show.plot){
    print(m)
  }
}

PlotGeneAcrossTime <- function(dat.sub, tmin = -24, tmax = 54,
                               jylab = "exprs",
                               jxlab = "Time [h]",
                               exprs.colname = "exprs", 
                               var.colname = "exprs.var.log2",
                               dotted.dbl = TRUE,
                               do.t.test = FALSE,
                               show.zero = FALSE,
                               jsize = 6,
                               jtitle){
  # Plot gene across time, should be dat.long filtered by gene
  # if var.colname = NA then do not show error bars
  # if do t-test, then use perform test
  
  if (dotted.dbl){
    # rename NSD 24 and beyond as "Double Plotted"
    dbl_str = "Double"
    dat.sub[[dbl_str]] <- as.factor(mapply(function(trtmnt, time){
      if (trtmnt == "NSD"){
        if (time > 24){
          return("Double")
        } else if (time == 24){
          return("Common")
        }
      }
      return("")
    }, dat.sub$trtmnt, dat.sub$time))
  }
  
  if (missing(jtitle)){
    jtitle <- dat.sub$gene[[1]]
  }
  if (is.na(var.colname)){
    var.colname <- NA
  } else {
    dat.sub$ymax <- dat.sub[[exprs.colname]] + sqrt(dat.sub[[var.colname]])
    dat.sub$ymin <- dat.sub[[exprs.colname]] - sqrt(dat.sub[[var.colname]])
  }
  time_str <- "time"
  trt_str <- "trtmnt"
  m <- ggplot(dat.sub, aes_string(x = time_str, y = exprs.colname, colour = trt_str, group = trt_str, shape = trt_str)) + 
    geom_point(size = 2)
  # add error bars
  if (!is.na(var.colname)){
    limits <- aes(ymax = ymax, ymin = ymin)
    m <- m + geom_errorbar(limits, width=0.25)
  }
  # add line on NSD, and make SD larger points
  if (!dotted.dbl){
    m <- m + geom_line(data = subset(dat.sub, trtmnt == "NSD"),
                       aes_string(x = time_str, y = exprs.colname, colour = trt_str, group = trt_str, shape = trt_str))
  } else {
    # Common allows solid and dotted lines to "connect"
    m <- m + geom_line(data = subset(dat.sub, trtmnt == "NSD" & (Double != "Double" | Double == "Common")),
                       aes_string(x = time_str, y = exprs.colname, colour = trt_str, group = trt_str, shape = trt_str), linetype = 1)
    m <- m + geom_line(data = subset(dat.sub, trtmnt == "NSD" & (Double == "Double" | Double == "Common")),
                       aes_string(x = time_str, y = exprs.colname, colour = trt_str, group = trt_str, shape = trt_str), linetype = 2)
  }
  m <- m + geom_point(data = subset(dat.sub, trtmnt == "SD"), 
                      aes_string(x = time_str, 
                                 y = exprs.colname, 
                                 colour = trt_str, 
                                 group = trt_str, 
                                 shape = trt_str),
                      size = jsize)   
  m <- m + ylab(jylab) + xlab(jxlab) + theme_bw(24) + scale_x_continuous(breaks=seq(tmin, tmax, 12)) + ggtitle(jtitle)
  if (show.zero){
    m <- m + expand_limits(y = 0)
  }
  return(m)
}

PlotGeneAcrossTimeLong <- function(dat, 
                                   jylab = "exprs",
                                   jxlab = "Time [h]",
                                   exprs.colname = "exprs", 
                                   var.colname = "exprs.var.log2",
                                   jtitle="",
                                   jlegend.pos = "top",
                                   jxlim = c(-6, 78)){
  time.str <- "time"
  trtmnt.str <- "trtmnt"
  # Plot gene across time, like PlotGeneAcrossTime but for entire time from 0 to 72. 0 to 24 is NSD, 24 to 72 is SD.
  dat$ymax <- dat[[exprs.colname]] + sqrt(dat[[var.colname]])
  dat$ymin <- dat[[exprs.colname]] - sqrt(dat[[var.colname]])
  m <- ggplot(dat, aes_string(x = time.str, y = exprs.colname)) + 
    geom_point(aes_string(x = time.str, y = exprs.colname, colour = trtmnt.str), size = 4) + geom_line(linetype = "dashed")
  # add axes and title
  m <- m  + ggtitle(jtitle) + xlab(jxlab) + ylab(jylab) + scale_x_continuous(breaks = seq(jxlim[1], jxlim[2], 6), limits = jxlim)
  # add vertical lines
  m <- m + geom_vline(xintercept = seq(0, 78, 24), linetype = "dotted") + theme_bw(24) + theme(legend.position = jlegend.pos)
  # add error bars
  limits <- aes(ymax = ymax, ymin = ymin)
  m <- m + geom_errorbar(limits, width=0.25)
  return(m)
}

GetExprsDist <- function(dat.sub, lnorm){
  common.times <- intersect(subset(dat.sub, trtmnt == "NSD")$time, subset(dat.sub, trtmnt == "SD")$time)
  exprs.dist <- subset(dat.sub, time %in% common.times) %>%
    group_by(gene, time) %>%
    summarise(exprs.dist = diff(exprs) ^ lnorm)
  return(exprs.dist)
}

PlotByDistance <- function(dat.sub, lnorm = 2, jxlab = "time", jylab = "Exprs diff (SD - NSD, log2)"){
  # 1, 2, Inf
  exprs.dist <- GetExprsDist(dat.sub, lnorm)
  gene <- exprs.dist$gene[[1]]
  m <- ggplot(exprs.dist, aes(x = time, y = exprs.dist)) + 
    geom_point() + 
    geom_line() + 
    theme_bw(24) + 
    xlab(jxlab) + 
    ylab(jylab) + 
    ggtitle(paste("Gene:", gene, "Lnorm:", lnorm)) 
  return(m)
}  


PlotSCGene <- function(fits, g, wake.collapsed, filt.time, low.pass.filt.time, dat.long.shift){
  jdo.lowpass <- TRUE
  params.sleep.gene <- GetParams.sleep(subset(fits, gene == g), "fit.sleep", append.bic = FALSE)
  # get same cnames for mix model
  sleep.cnames <- c("init", "U", "tau.w", "L", "tau.s", "tau.mrna")
  params.sleep.gene.frommix <- subset(fits, gene == g)$fit.mix[[1]]
  params.sleep.gene.frommix <- as.numeric(params.sleep.gene.frommix[names(params.sleep.gene.frommix) %in% sleep.cnames])
  # change to numeric but keep
  names(params.sleep.gene.frommix) <- sleep.cnames


  circ.cnames <- c("intercept", "amp", "phase")
  params.circ.gene.frommix <- subset(fits, gene == g)$fit.mix[[1]]
  params.circ.gene.frommix <- as.numeric(params.circ.gene.frommix[names(params.circ.gene.frommix) %in% circ.cnames])
  names(params.circ.gene.frommix) <- circ.cnames
  # set intercept to init
  params.circ.gene.frommix["intercept"] <- params.sleep.gene.frommix["init"]

  params.mix.gene <- GetParams.mixed(subset(fits, gene == g), remove.na=TRUE, mix.name = "fit.mix")

  S.pred <- S.process.collapsed(params.sleep.gene.frommix, wake.collapsed, filter.times = filt.time, do.lowpass = jdo.lowpass, low.pass.filter.times = low.pass.filt.time)
  C.pred <- CosSine(params.circ.gene.frommix, time = filt.time, ampphase=TRUE)
  Mix.pred <- WeightedSCircadian(params.mix.gene, wake.collapsed, filt.time, ampphase = TRUE,
                                 include.mix.weight = FALSE, include.intercept = FALSE,
                                 do.lowpass = jdo.lowpass, low.pass.filter.times = low.pass.filt.time)

  dat.actual <- dcast(subset(dat.long.shift, gene == g), formula = "gene ~ time + SampID", value.var = "exprs")
  dat.actual <- subset(dat.actual, select = -gene)
  y.actual <- as.numeric(dat.actual)
  tvec <- as.numeric(sapply(colnames(dat.actual), function(x) strsplit(x, "_")[[1]][[1]]))
  plot(filt.time, C.pred, ylim = range(S.pred, C.pred, y.actual), type = 'l', col = 'blue', main = paste("gene", g), lty = 3, xlab = "ZT [h]", ylab = "log2 mRNA abundance",xaxt="n")
  axis(1, at=seq(0, 78, 6), labels=seq(0, 78, 6))
  lines(filt.time, S.pred, type = 'l', col = 'red', lty = 3)
  lines(filt.time, Mix.pred, type = 'l', col = 'black', lwd = 2.5)
  legend("bottomright", inset=.05,
         c("C", "S"), col = c("blue", "red"), lty = 3, horiz=TRUE)
  points(tvec, y.actual, pch = 16)
}
