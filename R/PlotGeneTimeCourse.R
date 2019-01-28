# Jake Yeung
# Date of Creation: 2019-01-28
# File: ~/projects/SleepDepAnalysis/R/PlotGeneTimeCourse.R
# Plot genes over time with EEG data and models

source("R/functions/LoadData.R")
source("R/functions/FitFunctions_Downstream.R")

LoadPrimetimeObjs()


# Plot some genes ---------------------------------------------------------

g <- "Cry1"
print(PlotBestFit(subset(dat.long.shift, gene == g), subset(fits, gene == g), filt.time, g, wake.collapsed, low.pass.filt.time, dat.pred = NA, dat.eeg = dat.eeg.plot))



# Plot some motif activities ----------------------------------------------



library(ggplot2)
library(ggrepel)
library(dplyr)

source("R/functions/MaraDownstream.R")
source("R/functions/PlotFunctions.R")
source("R/functions/PlotFunctions.auxiliary.R")
source("R/functions/FitFunctions.R")
source("R/functions/LoadData.R")

# Load EEG data ---------------------------------------------------------------

do.linear <- FALSE
log.eps <- 1
do.test <- TRUE
if (do.test){
  print(paste("Testing mode:", do.test))
}
jtswitch <- 33
jform <- "switch"

use.merged <- TRUE
if (!use.merged){
  max.time <- 72
} else {
  max.time <- 78
}
zt0tozt24 <- FALSE

# HALF LIFE LIMITS
min.hl <- 1/3
max.hl <- 24
min.hl.mrna <- 1/3
max.hl.mrna <- 24

lambda <- 1000
jmaxiter <- 3000

# LoadPrimetimeObjs()


# Plot motifs -------------------------------------------------------------

plot.merged.models <- TRUE
top.n.motifs <- 45
weight.cutoff <- 0
merge.models <- TRUE
use.atacseq.peaks <- TRUE  # 15kb from promoter
merge.all.sleep <- TRUE
normalize.sitecounts <- TRUE
center.sitecounts <- FALSE  # center the sitecounts across columns??


col.i <- !grepl("^fit.", colnames(fits))
bic.cols <- colnames(fits)[grepl("^bic.", colnames(fits))]
fits.bic <- melt(fits[, col.i], id.vars = c("gene", "model"), measure.vars = bic.cols, variable.name = "model.names", value.name = "bic") %>%
  group_by(gene) %>%
  mutate(weight = exp(-0.5 * bic) / sum(exp(-0.5 * bic)))
genes.all <- as.character(unique(fits$gene))
jmodels <- unique(fits.bic$model)
genes.lst <- lapply(jmodels, function(m){
  genes <- as.character(unique(subset(fits.bic, model == m & weight > weight.cutoff)$gene))
})
names(genes.lst) <- jmodels

genes.lst[["genes_all"]] <- genes.all

if (merge.models){
  if (!merge.all.sleep){
    # copied from run_mara_on_models..combine_models script
    genes.lst <- MergeLists(genes.lst, "ampfree.step", "mixedaf")
    genes.lst <- MergeLists(genes.lst, "circadian", "mix")
    genes.lst <- RemoveLists(genes.lst, "ampfree.step")
    genes.lst <- RemoveLists(genes.lst, "mixedaf")
    genes.lst <- RemoveLists(genes.lst, "circadian")
    genes.lst <- RemoveLists(genes.lst, "mix")
  } else {
    # copied from run_mara.R script
    genes.lst <- MergeLists(genes.lst, "mix", "mixedaf")
    genes.lst <- MergeLists(genes.lst, "mix_mixedaf", "sleep")
    genes.lst <- MergeLists(genes.lst, "mix_mixedaf_sleep", "ampfree.step")
    genes.lst <- MergeLists(genes.lst, "mix_mixedaf_sleep_ampfree.step", "circadian")
    genes.lst <- RemoveLists(genes.lst, "mix")
    genes.lst <- RemoveLists(genes.lst, "mixedaf")
    genes.lst <- RemoveLists(genes.lst, "sleep")
    genes.lst <- RemoveLists(genes.lst, "ampfree.step")
    genes.lst <- RemoveLists(genes.lst, "circadian")
  }
}

# match TFs
tfs <- GetTFs(split.commas = FALSE, get.mat.only = TRUE)

# do.mean <- FALSE
# do.center <- TRUE
# gene.lab <- "genes_sleep"

# outmain <- paste0("/home/yeung/data/sleep_deprivation/mara_outputs", marasuffix)
# dir.create(outmain)

if (plot.merged.models){
  gene.labs <- c("mix_mixedaf_sleep_ampfree.step_circadian")
} else {
  gene.labs <- names(genes.lst)[which(names(genes.lst) != "flat")]
}

gene.lab <- c("mix_mixedaf_sleep_ampfree.step_circadian")

print(gene.lab)
# outdir <- file.path(outmain, paste0(gene.lab, ".mean.", do.mean, ".center.", do.center))

# suffix <- paste0("mean.", do.mean,  ".centered.", do.center,  marasuffix, "/sleep_deprivation_gene_exprs_all")
# act.dir <- file.path(outdir, suffix)
# act.dir <- ""
# assertthat::assert_that(dir.exists(act.dir))

# act.zscore.lst <- LoadMaraOutput("sleep_deprivation_gene_exprs_all")
act.zscore.lst <- LoadMaraOutput()
act.long <- act.zscore.lst$act.long
zscores <- act.zscore.lst$zscores
zscores$motif <- factor(as.character(zscores$motif), levels = zscores$motif)
# zscores$label <- mapply(function(m, z) ifelse(z > 1.68, m, NA), as.character(zscores$motif), zscores$zscore)
zscores$label <- mapply(function(m, z) ifelse(z > 1.46, m, NA), as.character(zscores$motif), zscores$zscore)


m.zscores <- ggplot(zscores, aes(x = motif, y = zscore, label = label)) + geom_point() + geom_text_repel(size = 7) +
  theme_bw(24) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  xlab("Index") +
  ylab("Zscore from Mean Activity")
print(m.zscores)

jmotif <- "SRF.p3"

# which model is SRF in?
fit.motif <- subset(act.long, gene == jmotif) %>%  # TODO for all motifs later and write to table in separate script
  group_by(gene) %>%
  do(fit.sleep = FitProcessS(., wake.collapsed, exprs.cname = "exprs", time.cname = "time",
                             condensed=TRUE, pseudo = 0, min.hl = min.hl, max.hl = max.hl,
                             do.lowpass=TRUE, min.hl.mrna = min.hl.mrna, max.hl.mrna = max.hl.mrna,
                             low.pass.filter.times = low.pass.filt.time, jlambda = lambda, jmaxit = jmaxiter),
     fit.flat = FitFlat(., use.weights = FALSE, get.bic = TRUE, condensed = TRUE),
     fit.circadian = FitRhythmic(., T.period = 24, use.weights=FALSE, get.bic=TRUE, condensed=TRUE),
     fit.ampfree.step = FitRhythmic.FreeAmp(., AmpFunc = AmpFunc, T.period = 24, tswitch = jtswitch, form = jform, include.intercept = TRUE, jmaxit = jmaxiter),
     fit.mixedaf = FitWeightedSFreeAmp(., AmpFunc = AmpFunc, wake.collapsed = wake.collapsed,
                                       exprs.cname = "exprs", time.cname = "time", T.period = 24,
                                       tswitch=jtswitch, form=jform, pseudo = 0, min.hl = min.hl,
                                       max.hl = max.hl, include.intercept = FALSE, do.lowpass=TRUE,
                                       min.hl.mrna = min.hl.mrna, max.hl.mrna = max.hl.mrna,
                                       low.pass.filter.times = low.pass.filt.time, jlambda = lambda, jmaxit = jmaxiter),
     fit.mix = FitWeightedSCircadian(., wake.collapsed, condensed = TRUE, pseudo = 0,
                                     min.hl = min.hl, max.hl = max.hl, include.mix.weight = FALSE,
                                     include.intercept = FALSE, do.lowpass=TRUE, min.hl.mrna = min.hl.mrna,
                                     max.hl.mrna = max.hl.mrna, low.pass.filter.times = low.pass.filt.time, jlambda = lambda, jmaxit = jmaxiter))
fit.motif <- AddBICToDat(fit.motif)
fit.motif$model <- SelectBestModel(fit.motif, colnames(fit.motif))  # SRF is sleep gene

act.means <- act.long %>%
  group_by(gene, time) %>%
  summarise(exprs = mean(exprs), sem = sqrt(sum(sem ^ 2)))

# indir <- "plots/SRF_Figure_7"
params.sleep <- SummarizeParameters(fits, "sleep")
# jgenes <- c("Egr2", "Arc", "Fos", "Egr1", "Nr4a1", "Plekhg4", "Serinc2", "Junb", "Car12", "Npas4", "Fosl2", "Egr3", "Rasl11a")
jgenes <- params.sleep$gene[1:15]
# Plekhg4, Npas4 is negative
# dir.create(indir)
# pdf(file.path(indir, paste0("SRF_Figure_7_plots-", Sys.Date(), ".pdf")), useDingbats = FALSE)
m.zscores
PlotMara.withEEG(subset(act.means, gene == jmotif), dat.eeg.plot, jtitle = jmotif)
# plot hits
for (g in jgenes){
  print(PlotBestFit(subset(dat.long.shift, gene == g), subset(fits, gene == g), filt.time, g, wake.collapsed, low.pass.filt.time, dat.pred = NA, dat.eeg = dat.eeg.plot))
}
# dev.off()
