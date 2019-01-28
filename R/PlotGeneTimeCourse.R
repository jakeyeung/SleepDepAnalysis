# Jake Yeung
# Date of Creation: 2019-01-28
# File: ~/projects/SleepDepAnalysis/R/PlotGeneTimeCourse.R
# Plot genes over time with EEG data and models

# source("R/functions/LoadData.R")
# source("R/functions/FitFunctions_Downstream.R")
#
# LoadPrimetimeObjs()
#
# # Plot some genes ---------------------------------------------------------
#
# glist <- c("Cry1", "Hif3a", "Rasl11a")
#
# pdf("/tmp/plots.pdf", useDingbats = FALSE)
# for (g in glist){
#   print(PlotBestFit(subset(dat.long.shift, gene == g), subset(fits, gene == g), filt.time, g, wake.collapsed, low.pass.filt.time, dat.pred = NA, dat.eeg = dat.eeg.plot))
# }
# dev.off()
