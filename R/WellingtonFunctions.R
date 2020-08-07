# Jake Yeung
# Tue Jul  4 09:13:00 2017


TimeFromSamp <- function(s){
  # "ZT03_NSD" -> 3
  split1 <- strsplit(s, split = "_")[[1]][[1]]
  split2 <- strsplit(split1, split = "ZT")[[1]][[2]]
  return(as.numeric(split2))
}

SleepFromSamp <- function(s){
  # "ZT03_NSD" -> NSD
  split1 <- strsplit(s, split = "_")[[1]][[2]]
  return(split1)
}