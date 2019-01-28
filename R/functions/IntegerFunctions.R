collapse.integers <- function(x){
  # Collapse so there are no gaps
  # c(1,3,5,7) -> c(1,2,3,4)
  # c(1,3,3,5) -> c(1,2,2,3)
  # c(9,8,5,5,3,2,1) -> c(6,5,4,4,3,2,1)
  x.collapsed <- rep(NA, length(x))
  indxs <- sort.int(x, index.return=TRUE, method="shell")$ix
  j <- 0
  for (i in seq(length(indxs))){
    indx <- indxs[i]
    x.now <- x[indx]
    if (i > 1){
      x.prev <- x[indxs[i - 1]]
    } else {
      x.prev <- Inf
    }
    # check if previous index is same as current index
    if (x.now == x.prev){
      x.collapsed[indx] <- j
    } else {
      j <-  j + 1
      x.collapsed[indx] <- j
    }
  }
  return(x.collapsed)
}