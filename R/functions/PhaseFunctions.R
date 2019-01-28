# copied from NcondsFunctions.R
DiffPhase <- function(phase1, phase2, period=24){
  # Get phase1 and phase2 difference, modulo period
  jdiff <- abs(diff(c(phase1, phase2)))
  return(min(period - jdiff, jdiff))
}

MaxDiffPhaseVector <- function(phases, period=24, get.i = FALSE){
  phasediff.max <- 0
  if (get.i) ijpair <- c(0, 0)
  for (phase.i in phases){
    for (phase.j in phases){
      jdiff <- DiffPhase(phase.j, phase.i, period)
      if (jdiff > phasediff.max){
        phasediff.max <- jdiff
        if (get.i) ijpair <- c(which(phases == phase.i), which(phases == phase.j))
      }
    }
  }
  if (!get.i){
    return(phasediff.max)
  } else {
    return(ijpair)
  }
}

FilterDatMaxPhase <- function(dat.sub){
  # Filter dat to keep only the two entries of max phase
  return(dat.sub[MaxDiffPhaseVector(dat.sub$phase, get.i = TRUE), ])
}