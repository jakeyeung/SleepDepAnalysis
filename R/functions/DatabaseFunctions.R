LoadDatabase <- function(sqlite3.path){
  library(dplyr)
  library(methods)
  # Load table 
  # more info: https://cran.r-project.org/web/packages/dplyr/vignettes/databases.html
  print(paste("Loading sqlite3 database:", sqlite3.path))
  atacseq.db <- src_sqlite(sqlite3.path, create=F)
  tblname <- src_tbls(atacseq.db)[[1]]
  atacseq.tbl <- tbl(atacseq.db, sql(paste0("SELECT * FROM ", tblname)))
  return(atacseq.tbl)
}


GetGenesFromDb <- function(motevo.tbl, genes){
  source("~/projects/tissue-specificity/scripts/functions/ListFunctions.R")
  print("Getting genes from database")
  start <- Sys.time()
  N.sub.lst <- expandingList()
  for (jgene in genes){
    N.long.filt.query <- filter(motevo.tbl, gene == jgene)  # peaks are not indexed, so dont take them
    N.sub.tmp <- collect(N.long.filt.query, n = Inf)
    N.sub.lst$add(N.sub.tmp)
  }
  N.long.filt <- N.sub.lst$as.list()
  N.long.filt <- bind_rows(N.long.filt)
  print(Sys.time() - start)
  return(N.long.filt)
  # rm(N.sub.tmp, N.sub.lst)  # worth it? 
}

