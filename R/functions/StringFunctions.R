ExtractChromoStartEndGene <- function(chromostartendgene.str){
  chromostartendgene <- strsplit(chromostartendgene.str, "\\s+")[[1]]
  jchromo <- chromostartendgene[[1]]
  jstart <- chromostartendgene[[2]]
  jend <- chromostartendgene[[3]]
  jgene <- chromostartendgene[[4]]
  jdist <- chromostartendgene[[5]]
  return(list(chromo = jchromo, start = jstart, end = jend, gene = jgene, dist = jdist))
}

split_str_by_index <- function(target, index) {
  index <- sort(index)
  substr(rep(target, length(index) + 1),
         start = c(1, index),
         stop = c(index -1, nchar(target)))
}

#Taken from https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
interleave <- function(v1,v2)
{
  ord1 <- 2*(1:length(v1))-1
  ord2 <- 2*(1:length(v2))
  c(v1,v2)[order(c(ord1,ord2))]
}

insert_str <- function(target, insert, index) {
  # http://stackoverflow.com/questions/13863599/insert-a-character-at-a-specific-location-in-a-string
  insert <- insert[order(index)]
  index <- sort(index)
  paste(interleave(split_str_by_index(target, index), insert), collapse="")
}