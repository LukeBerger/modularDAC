# other

#' Impute missing values in a data matrix using half the row minimum
#' @param dat a numeric matrix with p features (rows) and n samples (columns)

#' @return the data matrix with NAs filled


#' @keywords internal
.halfmin_impute <- function(dat) {
  halfmin <- matrix(matrixStats::rowMins(dat, na.rm = TRUE) / 2,
                    nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat)
  )
  # the imputed value is constant along each row
  stopifnot(all.equal(
    matrixStats::rowMeans2(halfmin),
    matrixStats::rowMins(dat, na.rm = TRUE) / 2
  ))
  dat[is.na(dat)] <- halfmin[is.na(dat)]
  return(dat)
}
