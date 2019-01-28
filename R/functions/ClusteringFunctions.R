DistFromClstr <- function(i, data, clstr){
  # stolen from https://stackoverflow.com/questions/10075122/ordering-clustered-points-using-kmeans-and-r
  #Extract cluster and center
  rows.keep <- clstr$cluster == i
  dt <- data[rows.keep, ]
  ct <- clstr$centers[i,]
  #Calculate distances
  dist <- apply((dt[, ] - ct)^2, 1, sum)
  return(dist)
}


RearrangeCluster <- function(clstr, clstr.hash, clstr.order = NA){
  clstr.ordered <- sapply(as.character(clstr), function(x) clstr.hash[[x]])
  if (!is.na(clstr.order)){
    clstr.ordered <- factor(as.character(clstr.ordered), levels = clstr.order)
  }
  return(clstr.ordered)
}