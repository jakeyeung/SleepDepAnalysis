#source("R/functions/PlotFunctions.auxiliary.R")

LoadMaraOutput <- function(act.dir="data/sleep_deprivation_gene_exprs_all"){
  act.f <- file.path(act.dir, "Activities")
  zscores.f <- file.path(act.dir, "Zscores")
  cnames.f <- file.path(act.dir, "Colnames")
  se.f <- file.path(act.dir, "StandardError")
  
  act.mat <- read.table(act.f, header = FALSE, sep = "\t")
  zscores <- read.table(zscores.f, header = FALSE, sep = "\t", col.names = c("motif", "zscore"))
  se.mat <- read.table(se.f, header = FALSE, sep = "\t")
  
  zscores <- zscores[order(zscores$zscore, decreasing = TRUE), ]
  # get colnames
  cnames <- readLines(cnames.f)
  # add "gene" to first cname
  cnames <- c("gene", cnames)
 
  act.long <- MatToLong(act.mat, cnames = cnames, jid.vars = "gene", jvar.name = "sample", jval.name = "exprs")
  se.long <- MatToLong(se.mat, cnames = cnames, jid.vars = "gene", jvar.name = "sample", jval.name = "sem")
  act.long <- inner_join(act.long, se.long)
  return(list(act.long=act.long, zscores=zscores))
}

MatToLong <- function(dat.mat, cnames, jid.vars = "gene", jvar.name = "sample", jval.name = "exprs"){
  # dat.mat -> long format, renames colnames(dat.mat) to cnames
  colnames(dat.mat) <- cnames
  
  dat.long <- melt(dat.mat, id.vars = jid.vars, variable.name = jvar.name, value.name = jval.name)
  dat.long$time <- sapply(dat.long$sample, LongSampToTime)
  dat.long$samp <- sapply(dat.long$sample, LongSampToRep)
  return(dat.long)
}

LongSampToTime <- function(s){
  # ZT_00_1 -> 0
  as.numeric(strsplit(as.character(s), "_")[[1]][[2]])
}

LongSampToRep <- function(s){
  # ZT_00_1 -> 1
  as.numeric(strsplit(as.character(s), "_")[[1]][[3]])
}

PlotMara <- function(jsub, add.sd.period = TRUE){
  # jsub is summarised across samples (expect exprs and sem)
  decimal.places <- 1
  
  jxlim <- range(jsub$time)
  
  m <- ggplot(data = jsub, aes(x = time, y = exprs, ymin = exprs - sem, ymax = exprs + sem)) + 
    geom_line() + 
    geom_errorbar(alpha = 0.25)
    
  m <- m + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    # scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6)) +
    geom_vline(xintercept = seq(0, 78, 24), linetype = "dashed", colour = "black")
  m <- m + scale_x_continuous(limits = jxlim, breaks = seq(jxlim[1], jxlim[2], 6), expand = c(0, 0)) + 
    scale_y_continuous(labels = fmt_dcimals(1))
  if (add.sd.period){
    m <- AddSDPeriod(m)
  }
  return(m)
}

PlotMara.withEEG <- function(jsub, dat.eeg.plot, jtitle, labsize = 20, ysize=20){
  if (missing(jtitle)){
    jtitle <- unique(as.character(jsub$gene))
  }
  jxlim <- range(jsub$time)
  m.gene <- PlotMara(jsub)  +  
    ylab("TF activity [A.U.]") +
    AddDarkPhases(alpha = 0.1) +
    ggtitle(jtitle) + 
    theme(legend.position= "top")
  
  # calculate ndecs required for PlotEeg
  ylabs <- GetYMajorSource(m.gene)
  # ylabs <- ggplot_build(m.gene)$layout$panel_params[[1]]$y.major_source
  # if (is.null(ylabs)){
  #   warning("Ylabs is NULL perhaps your hacky script (ggplot_build(m)$layout$panel_params[[1]]$y.major_source) has failed")
  # }
  ndecs <- max(sapply(ylabs, NDecimals), na.rm = TRUE)  # PlotFunctions.auxiliary.R
  if (min(ylabs) < 0){
    ndecs <- ndecs + 1  # account for negative sign
  }
  m.eeg <- PlotEeg(dat.eeg.plot, jxlim = jxlim, n.decs = ndecs)
  
  jlay <- matrix(c(rep(1, 4), 2), ncol = 1)
  # print("My labels")
  # print(c(paste0("-0.", paste(rep(0, ndecs), collapse="")), paste0("-5.", paste(rep(0, ndecs), collapse=""))))
  jticks <- c(paste0("-0.", paste(rep(0, ndecs), collapse="")), paste0("-5.", paste(rep(0, ndecs), collapse="")))
  multiplot(m.gene + xlab("") + scale_x_continuous(breaks = NULL) + 
              theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize), title=element_text(size=0.4*labsize)), 
            m.eeg + theme(axis.text=element_text(size=labsize), axis.title=element_text(size=labsize), axis.text.y=element_text(size=ysize)), layout = jlay)
}
