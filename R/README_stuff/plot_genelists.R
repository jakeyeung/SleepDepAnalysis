# 2016-04-04


# Source ------------------------------------------------------------------

#source("scripts/functions/PlotFunctions.R")

load("Robjs/dat.long.htseq.redo.Robj")

# Load --------------------------------------------------------------------



inf <- "/home/yeung/data/sleep_deprivation/gene_lists/78_genes.cut.txt"
gene.list <- read.table(inf, stringsAsFactors = FALSE)

# Make lower case on non-first letter
gene.list <- sapply(gene.list, function(jgene){
  paste(toupper(substr(jgene, 1, 1)), tolower(substr(jgene, 2, nchar(jgene))), sep="")
})

pdf("plots/gene_lists/qc_78_genes.pdf")
for (g in gene.list){
  jsub <- subset(dat.long, gene == g)
  if (nrow(jsub) > 0){
    print(PlotGeneAcrossTime(subset(dat.long, gene == g)))
  } else {
    warning(paste("Skipping", g))
  }
}
dev.off()
warnings()
