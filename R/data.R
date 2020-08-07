#' Long format of data containing genes over time "dat.long.cleaned"
#'
#' @format A data frame with 797754 rows and 7 columns
#' \describe{
#'   \item{gene}{Gene name}
#'   \item{SampID}{Sample ID}
#'   \item{exprs}{Log2 + 1 transformed and batch-corrected expression}
#'   \item{trtmnt}{Treatment NSD or SD}
#'   \item{time}{0 to 78 hours}
#'   \item{batch}{Either 1; 2; or 1,2 (technical replicate merged)}
#'   \item{samp}{biological replicates}
#' }
"dat.long.cleaned.techmerged_zt24assigned"

#' 
#'
#' @format Dataframe 91801 rows, 4 columns
#' \describe{
#'   \item{time.shift}{Time shifted so time starts from -24, to give prediction at t=0}
#'   \item{wake}{Wake or not wake, scored manually, averaged across mice}
#'   \item{N}{Number of mice}
#'   \item{w.smooth}{Smoothed EEG signal from 0 to 5. 5 means mice was awake for last 5 minutes.}
#' }
"wake.df.method.mode"

#'
#' @format Fits output
#' \describe{
#'   \item{}{}
#' }
"fits_Rcpp_bugfix.maxiter.2000.LogLBicCorrected.model.sleep_mix_mixedaf_ampfreestep.lambda.1000.dolinear.FALSE.minhl.0.33"

#'
#' @format Fits output
#' \describe{
#'   \item{}{}
#' }
"fits_Rcpp_sleuthfilt.maxiter.3000.model.LogLBicCorrected.model.sleep_mix_mixedaf.lambda.1000.dolinear.FALSE.minhl.0.33.RData"

#'
#' @format Object containing input for kmeans analysis
#' \describe{
#'   \item{}{}
#' }
"kmeans_objs.Robj"

