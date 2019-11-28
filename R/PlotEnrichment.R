# Jake Yeung
# Date of Creation: 2018-09-05
# File: ~/projects/sleep_deprivation/R/primetime_checks/fishers_test_on_gene_lists.R
# Do fisher's test on gene lists provided by Franken


PlotEnrichment <- function(){
  library(ggplot2)
  library(ggrepel)
  library(dplyr)

  library(xlsx)

  source("R/functions/MaraDownstream.R")
  source("R/functions/PlotFunctions.R")
  source("R/functions/PlotFunctions.auxiliary.R")
  source("R/functions/FitFunctions.R")
  source("R/functions/LoadData.R")
  source("R/functions/MdlNames.R")
  source("R/functions/GetSwitchHash.R")

  # Affymetrix IDs found at
  # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL1261
  # GeneChip Mouse Expression Set 430
  # Downloaded to: /home/yeung/projects/sleep_deprivation/tables/affymetrix_table/GPL1261-56135.txt


  # Constants ---------------------------------------------------------------

  qval.cutoff <- 0.05
  # outdir <- "/home/yeung/projects/sleep_deprivation/plots/enrichment_analysis"
  outdir <- "/tmp/enrichment_analysis"
  dir.create(outdir)

  remove.flat.genes <- TRUE

  # Functions ---------------------------------------------------------------

  AssignGene <- function(g, genes.hash, switch.hash){
    # convert known bad genes like Gm129 -> Ciart
    if (g == ""){
      return("NA")
    }
    g <- ifelse(!is.null(switch.hash[[g]]), switch.hash[[g]], g)

    # Check if in genes.hash, otherwise parse out the /// text
    val <- genes.hash[[g]]
    val <- ifelse(is.null(val), FALSE, val)
    if (val == TRUE){
      return(g)
    } else {
      # if in form "Phgdh /// LOC666422 /// LOC666875 /// LOC669985 /// LOC671102 /// LOC673015 /// LOC675010" loop through ///
      gvec <- strsplit(g, " /// ")[[1]]
      # loop through gvec
      val <- FALSE
      for (G in gvec){
        valG <- genes.hash[[G]]
        if (!is.null(valG)){
          return(G)
        }
      }
    }
    # nothing matches return NA
    return(NA)
  }


  CountsMatToLong <- function(counts.mat){
    # rename models for plotting
    counts.mat.wide <- data.frame(FG = counts.mat[1, ],
                                  BG = counts.mat[2, ],
                                  model = colnames(counts.mat))
    counts.mat.long <- melt(counts.mat.wide, value.name = "counts", variable.name = "set", id.vars = "model") %>%
      group_by(set) %>%
      mutate(counts.norm = counts / sum(counts))
    counts.mat.long$model <- ConvertMdlNames(counts.mat.long$model, mdl.hash = MdlHash(), jlevels = MdlVals())
    return(counts.mat.long)
  }

  CheckExists <- function(jgene, capitalize=TRUE){
    if (capitalize){
      jgene <- tolower(jgene)
      jgene <- paste(toupper(substr(jgene, 1, 1)), substr(jgene, 2, nchar(jgene)), sep="")
    }
    print(jgene)
    print(subset(fits, grepl(jgene, gene)))
  }

  PrintNAGenes <- function(genes.fg){
    print(genes.fg[unique(names(genes.fg[is.na(genes.fg)]))])
  }

  # Load things -------------------------------------------------------------

  LoadPrimetimeObjs()

  if (remove.flat.genes){
    fits <- subset(fits, model != "flat")
  }

  genes.hash <- hash(fits$gene, rep(TRUE, nrow(fits)))
  # switch gene names that I know are incorrect
  switch.hash <- GetSwitchHash()


  # Load Franken tables -----------------------------------------------------

  tables.dir <- "data/tables/tables_from_franken.2018-09-04"
  bg.tbl.path <- "data/tables/affymetrix_table/GPL1261-56135.txt"
  dat.bg <- read.csv2(bg.tbl.path, header = TRUE, stringsAsFactors = FALSE, comment.char = "#", sep = "\t")
  dat.bg$gene <- sapply(dat.bg$Gene.Symbol, AssignGene, genes.hash, switch.hash)

  dat.maret_rhyth <- read.xlsx2(file.path(tables.dir, "Maret_Table_3_2032_rhythmic_under_baseline.xlsx"), sheetIndex = 1)
  dat.maret_circadian <- read.xlsx2(file.path(tables.dir, "Maret_Table_4_391_still_circadian_after_4_SDs.xlsx"), sheetIndex = 1)
  dat.maret_any <- read.xlsx2(file.path(tables.dir, "Maret_Table_5_343_affected_by_SD_independent_of_time-of-day.xlsx"), sheetIndex = 1)
  dat.mongrain <- read.table(file.path(tables.dir, "78genes_Mongrain2010.txt"), header = TRUE, stringsAsFactors = FALSE)

  # switch gene names that I know are incorrect
  gene.hash <- hash("Gm129", "Ciart")

  # rhythmic in baseline of cortex

  # Cortex comparisons ------------------------------------------------------


  dat.rhyth <- read.table("data/tables/tables_for_charlotte.2018-08-16/fit.rhyth.base.annot.txt", header = TRUE)

  dat.rhyth$qval <- p.adjust(dat.rhyth$pval)

  dat.rhyth <- left_join(dat.rhyth, subset(fits, select = c(gene, model)))

  dat.rhyth.filt <- subset(dat.rhyth, qval <= qval.cutoff)

  g <- "Egr2"
  print(PlotPoints(subset(dat.long.shift, gene == g)) + ggtitle(g) + xlab("Time [h]") + ylab("log2(Counts)"))

  # load hogenesch brain tissues
  braindir <- "data/tables/rhythmic_genes_brain"
  dat.brain <- lapply(list.files(path = braindir, pattern = "*.txt", full.names = TRUE), function(f){
    return(read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  })
  dat.brain <- bind_rows(dat.brain)

  # find overlap
  dat.brain.ol <- dat.brain %>%
    group_by(gene) %>%
    filter(max(pval) <= qval.cutoff) %>%
    filter(pval == max(pval))

  # Do fisher's exact test --------------------------------------------------

  # should sample from affymetrix genes to make sure you don't get significance...



  # genes.fg.raw <- as.character(dat.maret_circadian$Gene.Symbol)
  # genes.fg.raw <- as.character(dat.maret_rhyth$Gene.Symbol)
  #fgname <- "Maret_Circadian"

  jlist <- list("Maret_Circadian" = as.character(dat.maret_circadian$Gene.Symbol),
                "Maret_RhythBaseline" = as.character(dat.maret_rhyth$Gene.Symbol),
                "Maret_Any" = as.character(dat.maret_any$Gene.Symbol),
                "Mongrain78" = as.character(dat.mongrain$Gene_symbol))

  jlist[[paste0("RhythBaselineCortex.qval.", qval.cutoff)]] <- as.character(dat.rhyth.filt$gene)
  jlist[[paste0("RhythHogeneschBrain.pval.", qval.cutoff)]] <- as.character(dat.brain.ol$gene)

  pdf(file.path(outdir, paste0("enrichment_plots.qval.", qval.cutoff, ".removeflat.", remove.flat.genes, ".", Sys.Date(), ".pdf")))
  for (fgname in names(jlist)){
    genes.fg.raw <- jlist[[fgname]]
    # check genes.fg is in universe
    genes.fg <- sapply(genes.fg.raw, AssignGene, genes.hash, switch.hash)

    genes.unknown <- genes.fg[which(genes.fg == "NA" | is.na(genes.fg))]
    genes.fg <- genes.fg[which(genes.fg != "NA")]  # those genes did not make it!
    genes.bg <- fits$gene

    fits.fg <- table(subset(fits, gene %in% genes.fg)$model)
    fits.bg <- table(subset(fits, gene %in% genes.bg)$model)

    genes.fg.unknown <- data.frame(gene.key = names(genes.unknown),
                                   gene.val = genes.unknown)

    counts.mat <- as.matrix(bind_rows(fits.fg, fits.bg))
    # NA <- 0
    counts.mat[is.na(counts.mat)] <- 0

    out <- chisq.test(counts.mat)

    counts.mat.long <- CountsMatToLong(counts.mat)
    m <- ggplot(counts.mat.long, aes(x = model, y = counts.norm, fill = set)) + geom_bar(stat = "identity", position = position_dodge()) +
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("Model") + ylab("Fraction of genes") +
      ggtitle(paste(fgname, "\n# FG genes:", sum(fits.fg), paste("\n# BG genes:", paste(fits.bg, collapse=",")), "\nChi-Sqr Pval:", signif(out$p.value, 2)))
    print(m)
    colnames(counts.mat) <- ConvertMdlNames(colnames(counts.mat))
    write.table(counts.mat, file = file.path(outdir, paste0(fgname, ".removeflat.", remove.flat.genes, ".", Sys.Date(), ".txt")), quote = FALSE, sep = "\t", row.names=FALSE)
    write.table(genes.fg.unknown, file = file.path(outdir, paste0(fgname, ".removeflat.", remove.flat.genes, ".", Sys.Date(), ".MissingGenes.txt")), row.names = FALSE, quote = FALSE, sep = "\t")
  }
  dev.off()



  # Get average p-value if I sample from Affymetrix many times --------------

  set.seed(0)
  N <- 2000
  ntrials <- 1000
  out.vec <- rep(NA, ntrials)
  out.lst <- lapply(seq(ntrials), function(i){
    # print(i)
    fits.bg <- table(subset(fits, gene %in% genes.bg)$model)
    genes.fg <- sample(dat.bg$gene, size = N, replace = FALSE)
    fits.fg <- table(subset(fits, gene %in% genes.fg)$model)
    counts.mat <- as.matrix(rbind(fits.fg, fits.bg))
    out <- chisq.test(counts.mat)
    return(out$p.value)
  })

  pdf(file.path(outdir, paste0("background_pval_distrib.", "removeflat", remove.flat.genes, Sys.Date(), ".pdf")))
  plot(hist(unlist(out.lst)), col = 'lightblue', xlab = "P-value", main = "ChiSqr p-value distrib from sampling Affymetrix genes")
  dev.off()



}
