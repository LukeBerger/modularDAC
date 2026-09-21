#' Select soft-thresholding power from WGCNA::pickSoftThreshold output
#' @param sft a list, the output of WGCNA::pickSoftThreshold

#' @return an integer, the selected soft-thresholding power

#' @keywords internal
.sft_check <- function(sft) {
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- 6 # default
    message("Using the following power: ", beta)
  } else {
    message("Optimal power selected: ", beta)
  }
  return(beta)
}
