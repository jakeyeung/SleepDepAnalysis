GetTFs <- function(split.commas = TRUE, get.motifs = FALSE, get.mat.only = FALSE){
  # Vector containing gene names (may be comma separated), get gene list
  #
  # Args:
  # tf_vector: vector containing gene names (amy be comma separated)
  # with promoter (obtained from get_TFs_from_associations.py script)
  
  # define dirs
  # data.dir <- "/home/yeung/projects/tissue-specificity/data"
  data.dir <- "data"
  tf.fname <- "motifs_and_TFs.list"
  tf.path <- file.path(data.dir, tf.fname)
  
  tf.mat <- read.table(tf.path, header=FALSE, row.names = 1, sep='\t')
  # change gene names that are not well named
  # tf.mat$V2[grepl("Zfp161", tf.mat$V2)] <- "Zbtb14"
  # tf.mat$V2[grepl("Tcfap2b", tf.mat$V2)] <- "Tfap2b"
  tf.mat$V2 <- as.character(tf.mat$V2)
  tf.mat["TFAP2B.p2", ] <- "Tfap2b"
  tf.mat["ZFP161.p2", ] <- "Zbtb14"
  tf.mat["ZNF238.p2", ] <- "Zbtb18"
  tf.mat["ZNF238.p2", ] <- "Zbtb18"
  tf.mat["TFCP2.p2", ] <- "Tfcp2"
  tf.mat["TFAP4.p2", ] <- "Tfap4"
  tf.mat["bHLH_family.p2", ] <- gsub("Tcfe3", "Tfe3", tf.mat["bHLH_family.p2", ])  # replace Tcfe3 to Tfe3
  
  if (get.mat.only){
    return(tf.mat)
  }
  
  if (get.motifs == FALSE){
    tf_vector <- as.vector(tf.mat[, 1])
    if (split.commas){
      genes.list <- strsplit(tf_vector, split=",")
      # flatten list, return uniques
      genes.list <- unique(unlist(genes.list))
    } else {
      genes.list <- tf_vector
      return(genes.list)
    }
  } else {
    return(rownames(tf.mat))
  }
  return(genes.list)
}

GetGenesFromMotifs <- function(jmotif, tfs){
  # Return gene(s) that match a motif
  genes <- strsplit(tfs[jmotif, ], ",")[[1]]
  if (all(is.na(genes))){
    # try grep
    jgrep <- strsplit(jmotif, "\\.")[[1]][[1]]
    genes <- strsplit(tfs[grepl(jgrep, rownames(tfs)), ], ",")[[1]]
  }
  return(genes)
}

GetMotifFromGene <- function(gene, tfs){
  return(rownames(tfs)[grepl(gene, tfs[, 1])])
}
