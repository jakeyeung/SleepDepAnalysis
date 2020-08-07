# Jake Yeung
# Date of Creation: 2018-09-06
# File: ~/projects/sleep_deprivation/scripts/functions/ChangeMdlNames.R
# Change model names as a function so we can change things during plotting in a global way

MdlKeys <- function(){
  keys <- c("flat", "sleep", "circadian", "ampfree.step", "mix", "mixedaf")
  return(keys)
}

MdlVals <- function(){
  # vals <- c("F", "S", "C", "A", "S+C", "S+A")
  vals <- c("F", "S", "C_A", "A", "S+C", "S+C_A")
  return(vals)
}

MdlHash <- function(){
  # mdl.hash <- hash(list("flat" = "F",
  #                       "sleep" = "S",
  #                       "circadian" = "C",
  #                       "ampfree.step" = "A",
  #                       "mix" = "S+C",
  #                       "mixedaf" = "S+A"))
  mdl.hash <- hash(MdlKeys(),
		   MdlVals())
  return(mdl.hash)
}

ConvertMdlNames <- function(mdl.old, mdl.hash, jlevels){
  if (missing(mdl.hash)){
    mdl.hash <- MdlHash()
  }
  if (missing(jlevels)){
    jlevels <- MdlVals()
  }
  mdl.old <- as.character(mdl.old)
  mdl.new <- sapply(mdl.old, function(x){
   xnew <- mdl.hash[[x]]
   if (is.null(xnew)){
     warning(paste("Cannot convert", x, "using hash", mdl.hash, "defaulting NA"))
     xnew <- NA
   }
  return(xnew)
  })
  # set as factor with levels 
  mdl.new <- factor(mdl.new, levels = jlevels)
  return(mdl.new)
}
