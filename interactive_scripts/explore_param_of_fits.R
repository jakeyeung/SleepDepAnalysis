# Jake Yeung
# Date of Creation: 2021-11-02
# File: ~/projects/SleepDepAnalysis/interactive_scripts/explore_param_of_fits.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(dplyr)

library(SleepDepAnalysis)

# change this to your path 
indir <- "/Users/jakeyeung/projects/SleepDepAnalysis"

# load the object from the repo
inf.fit <- file.path(indir, "data/fits_Rcpp_sleuthfilt.maxiter.3000.model.LogLBicCorrected.model.sleep_mix_mixedaf.lambda.1000.dolinear.FALSE.minhl.0.33.RData")
load(inf.fit, v=T)

fits.all.fixed <- ungroup(fits.all.fixed)
# head(fits.all.fixed)

fits.sub <- subset(fits.all.fixed, model == "sleep")

fits.sub.summary <- SummarizeParameters(fits = fits.sub, jmodel = "sleep")

med.tau.w <- signif(median(fits.sub.summary$tau.w), digits = 3)
med.tau.s <- signif(median(fits.sub.summary$tau.s), digits = 3)
print(paste("Median of time constant for wake is:", med.tau.w))
print(paste("Median of time constant for sleep is:", med.tau.s))

outpdf <- paste0("/Users/jakeyeung/data/sleepdep_outputs/tau_values_histograms.", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)
ggplot(fits.sub.summary, aes(x = tau.w)) + 
  geom_histogram(bins = 100) + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = med.tau.w) + 
  ggtitle("Time Constant for Wake Amongst Model S Genes") + 
  xlab("Time Constant [h]") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(fits.sub.summary, aes(x = tau.w)) + 
  geom_histogram(bins = 100) + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = med.tau.w) + 
  ggtitle("Time Constant for Wake Amongst Model S Genes") + 
  xlab("Time Constant [h]") + 
  coord_cartesian(xlim = c(0.75, 24)) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


ggplot(fits.sub.summary, aes(x = tau.s)) + 
  geom_histogram(bins = 100) + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = med.tau.s) + 
  ggtitle("Time Constant for Wake Amongst Model S Genes") + 
  xlab("Time Constant [h]") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(fits.sub.summary, aes(x = tau.s)) + 
  geom_histogram(bins = 100) + 
  theme_bw() + 
  scale_x_log10() + 
  geom_vline(xintercept = med.tau.s) + 
  ggtitle("Time Constant for Wake Amongst Model S Genes") + 
  coord_cartesian(xlim = c(0.75, 24)) +
  xlab("Time Constant [h]") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# Write file to output ----------------------------------------------------

outtxt <- "/Users/jakeyeung/projects/SleepDepAnalysis/data/tables/parameter_estimates_from_fit/parameter_estimates_from_fit_sleep_genes.txt"
jgenes <- fits.sub.summary$gene
fits.sub.summary.pretty <- data.frame(fits.sub.summary)
rownames(fits.sub.summary.pretty) <- jgenes
fits.sub.summary.pretty$gene <- NULL

write.table(fits.sub.summary.pretty, file = outtxt, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


