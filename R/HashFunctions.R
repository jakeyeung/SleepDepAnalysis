# Jake Yeung
# HashFunctions.R
#  
# 2018-01-15

Vectorize(AddFromHash <- function(x, hsh, none = NA){
  if(x == "") return(none)
  if(!is.null(hsh[[as.character(x)]])){
    return(hsh[[as.character(x)]])
  } else {
    return(none)
  }
})
