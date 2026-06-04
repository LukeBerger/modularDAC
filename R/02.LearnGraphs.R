####################
### Learn Graphs ###
####################

# allow for dplyr piping
utils::globalVariables(".")

# graph learning

# takes a p x n matrix of data and inputs to SILGGM function
# then generates networks based on partial correlation between nodes

#' Learn a gene co-expression graph from a data matrix using SILGGM
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param method a character, the statistical inference method; one of 'B_NW_SL', 'D-S_NW_SL', 'D-S_GL', 'GFC_SL', or 'GFC_L'
#' @param global a logical, if TRUE global inference is performed across all feature pairs
#' @param alpha a numeric or numeric vector, the pre-specified FDR significance level(s) for edge inclusion (default 0.05)
#' @param fdr.filter a logical, if TRUE the adjacency matrix is filtered to retain only edges below max.fdr
#' @param max.fdr a numeric, the FDR q-value threshold above which an edge is set to zero
#' @param pos.cut a numeric, the minimum absolute partial correlation for a positive edge to be retained
#' @param neg.cut a numeric, the minimum absolute partial correlation for a negative edge to be retained
#' @param ... additional arguments passed to SILGGM::SILGGM

#' @return an igraph object

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom dplyr  %>%

#' @export
learn_SILGGM_graph <- function(x,
                               method  = "B_NW_SL",
                               global = T,
                               alpha = 0.05,
                               fdr.filter = T,
                               max.fdr = 0.05,
                               pos.cut = 0,
                               neg.cut = 0,
                               ...){
  if (!requireNamespace("SILGGM", quietly = TRUE)) {
    stop("Package SILGGM is required. Install with: install.packages('SILGGM')", call. = FALSE)
  }

  # run silggm
  silggm.output <- suppressMessages(
    SILGGM::SILGGM(x,
                   method = method,
                   alpha = alpha,
                   global = global,
                   ...)
  )
  # get partial.cor and zScore based on method
  if(method == "D-S_NW_SL"){
    partial.cor <- .upper_tri_vec(silggm.output$partialCor)
    z.score.partial.cor <- sapply(partial.cor, .pcor_zscore, nrow(x))
  }
  if(method == "B_NW_SL"){
    partial.cor <- .upper_tri_vec(silggm.output$partialCor)
    z.score.partial.cor <- .upper_tri_vec(silggm.output$z_score_partialCor)
  }

  # estimate the adjusted p-value of each edge
  qval.pcor <- sapply(z.score.partial.cor, .pvalue) %>%
    .upper_tri_to_matrix(., variable_names = colnames(x), diagl = 1) %>%
    .matrix_p_adjust(.)

  # build partial correlation adjacency matrix
  pcor.avg <- partial.cor %>%
    round(.,2) %>%
    .upper_tri_to_matrix(., variable_names = colnames(x), diagl = 1)
  adj.mat <- apply(pcor.avg, c(1,2), .abs_pcor_filter, pos.cut, neg.cut)

  # further filtering by significance
  if(fdr.filter){
    adj.mat[qval.pcor>=max.fdr] <- 0
  }

  # diagonal to 0 before creating igraph object
  diag(adj.mat) <- 0

  # convert weighted adjacency matrix to igraph object
  g  <- igraph::graph_from_adjacency_matrix(adj.mat,
                                            mode = "undirected",
                                            weighted = TRUE)

  return(g)
}

# helpers to simpleSILGGM graph taken from RSCGGM
# so that i can run it without loading the full package

#' Reconstruct the symmetric matrix from upper triangular vector
#'
#' @param upper_tri_values a numeric vector of the upper triangle of the matrix
#' @param variable_names a character vector of row and column names for the output matrix
#' @param diagl a numeric scalar or vector of length p, the diagonal value(s) of the output matrix

#' @return a symmetric matrix

#' @import methods utils

#' @keywords internal
.upper_tri_to_matrix <- function(upper_tri_values,
                                variable_names =NULL,
                                diagl=1){
  p <- (1 + sqrt(1 + 8 * length(upper_tri_values))) / 2
  if( (length(diagl)>1) & (length(diagl) != p) ) stop("invalid dignoal!")

  mat <- matrix(0, p, p)

  if(!is.null(variable_names)){
    row.names(mat) <- variable_names
    colnames(mat) <- variable_names
  }

  # fill the diagonal
  diag(mat) <- diagl

  # fill the upper triangular part
  mat[upper.tri(mat, diag = FALSE)] <- upper_tri_values

  # fill the lower triangular part (mirror the upper triangular part)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  return(mat)
}

#' Perform adjustments of p-values on a p x p matrix
#'
#' @param mx_p a p x p matrix, (i,j) represents the p-value of the partial correlation between node i and node j
#'
#' @return a p x p matrix, (i,j) represents the adjusted p-value of the partial correlation between node i and node j

#' @import methods utils

#' @keywords internal
.matrix_p_adjust <- function( mx_p ) {
  # initialize
  mx_q <- mx_p
  # adjust upper triangle
  mx_q[upper.tri(mx_q)] <-
    stats::p.adjust(mx_p[upper.tri(mx_p)], method = "BH")
  # copy to lower triangle
  mx_q[lower.tri(mx_q)] <-
    t(mx_q)[lower.tri(mx_q)]
  # some checks
  stopifnot(isTRUE(
    all.equal(
      mx_q[upper.tri(mx_q)],
      t(mx_q)[upper.tri(mx_q)]
    )
  ))
  stopifnot(!isTRUE(all.equal(mx_p,mx_q)))
  return( mx_q)
}

#' Take the upper triangular vector of a symmetric matrix
#'
#' @param mat the input symmetric matrix
#'
#' @return a numeric vector

#' @keywords internal
.upper_tri_vec <- function(mat){
  vec <- mat[upper.tri(mat)]
  if(length(vec) != (ncol(mat)*(ncol(mat)-1))/2) stop ("invalid length!")
  return(vec)
}

#' Compute the z-score of a partial correlation, see https://github.com/cran/SILGGM/blob/master/src/SILGGMCpp.cpp
#'
#' @param pcor the partial correlation
#' @param n the number of samples
#'
#' @return the z-score of the corresponding partial correlation

#' @keywords internal
.pcor_zscore <- function(pcor,n){
  std_new <- sqrt(.pow((1-.pow(pcor,2)),2)/n)
  return(pcor/std_new)
}

#' Raise of number to its n-th power
#'
#' @param x the input number
#' @param n the order of power
#'
#' @return a number

#' @keywords internal
.pow <- function(x,n){
  return(x^n)
}

#' Compute the two-sided p-value of a z-score
#'
#' @param z_score a z-score
#'
#' @return the p-value
#' @keywords internal
.pvalue <- function(z_score){
  return(2*stats::pnorm(q=abs(z_score), lower.tail=FALSE))
}

#' Filter edges whose effect size is small
#'
#' @param x the average partial correlation
#' @param pos_cut a threshold for positive partial correlation
#' @param neg_cut a threshold for negative partial correlation

#' @return the average value or 0

#' @keywords internal
.abs_pcor_filter <- function(x,
                             pos_cut,
                             neg_cut){

  if( (x) > 0 & (abs(x) <= abs(pos_cut)) ){
    return(0)
  }

  if( (x) < 0 & (abs(x) <= abs(neg_cut)) ){
    return(0)
  }

  return(x)
}

#' Learn a gene co-expression graph from a data matrix using WGCNA
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param min.sft a numeric between 0 and 1, the minimum R-squared threshold for soft-thresholding power selection
#' @param beta an integer, the soft-thresholding power for WGCNA::adjacency; if NULL it is selected automatically via pickSoftThreshold
#' @param cor.FN a character, the correlation function to use in WGCNA::adjacency; either 'bicor' or 'cor'
#' @param powers an integer vector, candidate soft-thresholding powers evaluated by WGCNA::pickSoftThreshold
#' @param threshold a numeric, the adjacency value above which an edge is retained; if NULL the 99.9th percentile is used

#' @return an igraph object

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats quantile

#' @export
learn_WGCNA_graph <- function(x,
                              min.sft=0.85,
                              beta=NULL,
                              cor.FN=c("bicor", "cor"),
                              powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                              threshold = NULL
) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package WGCNA is required. Install with: install.packages('WGCNA')", call. = FALSE)
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package matrixStats is required. Install with: install.packages('matrixStats')", call. = FALSE)
  }
  # handle arguments
  cor.FN <- match.arg(cor.FN)

  # correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=x,
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # check selected power
    beta <- .sft_check(sft)
  }

  # construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=x,
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # binarize based on threshold
  if(is.null(threshold)){threshold <- quantile(.upper_tri_vec(adj), .999)}
  binary.adj <- (adj > threshold) * 1

  # remove looped edges
  diag(binary.adj) <- 0

  # return as igraph object
  return(igraph::graph_from_adjacency_matrix(binary.adj, mode = "undirected"))
}

#' Learn a gene co-expression graph from a data matrix using the ARACNE algorithm
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param threshold a numeric, the mutual information value above which an edge is retained
#' @param eps a numeric, the data processing inequality threshold used by ARACNE to remove indirect edges

#' @return an igraph object

#' @importFrom igraph graph_from_adjacency_matrix

#' @export
learn_ARACNE_graph <- function(x, eps=0, threshold = 0.05) {
  if (!requireNamespace("minet", quietly = TRUE)) {
    stop("Package minet is required. Install with: install.packages('minet')", call. = FALSE)
  }

  # build mutual information matrix
  mim <- minet::build.mim(dataset = x)

  # apply ARANCE alogrithm
  aracne.mat <- minet::aracne(mim, eps=eps)

  # binarize and convert to graph object
  return(igraph::graph_from_adjacency_matrix((aracne.mat > threshold ) * 1, mode = "undirected"))
}

######################
### Compare Graphs ###
######################

#' Calculate the F1 score comparing a predicted graph against a true graph
#' @param g.true an igraph object, the true reference graph
#' @param g.pred an igraph object, the predicted graph to evaluate

#' @return a list of F1 score and its components

#' @importFrom igraph V as_edgelist

#' @export
calc_F1 <- function(g.true, g.pred) {
  # check that the graphs have the same nodes
  if (!all(igraph::V(g.true)$name %in% igraph::V(g.pred)$name)) {
    stop("True and Predicted graphs must contain the same nodes")
  }

  # extract edges as character vectors in a consistent format
  edges.true <- apply(igraph::as_edgelist(g.true), 1, function(x) paste(sort(x), collapse = " -- "))
  edges.pred <- apply(igraph::as_edgelist(g.pred), 1, function(x) paste(sort(x), collapse = " -- "))

  # count True Positives, False Positives, False Negatives
  tp <- sum(edges.pred %in% edges.true)
  fp <- sum(!(edges.pred %in% edges.true))
  fn <- sum(!(edges.true %in% edges.pred))

  # calc precision and recall
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))

  # calc F1 Score
  if (precision + recall == 0) {
    f1 <- 0 # to avoid zero divison error in edge case
  } else {
    f1 <- 2 * precision * recall / (precision + recall)
  }

  # return values as list
  return(list(
    F1 = f1,
    Precision = precision,
    Recall = recall,
    TruePos = tp,
    FalsePos = fp,
    FalseNeg = fn
  ))
}

# other

#' Impute missing values in a data matrix using half the row minimum
#' @param dat a numeric matrix with p features (rows) and n samples (columns)

#' @return the data matrix with NAs filled


#' @keywords internal
.halfmin_impute <- function(dat) {
  halfmin <- matrix(matrixStats::rowMins(dat, na.rm = TRUE) / 2,
                    nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat)
  )
  # just checking
  stopifnot(all.equal(
    matrixStats::rowMeans2(halfmin),
    matrixStats::rowMins(dat, na.rm = TRUE) / 2
  ))
  dat[is.na(dat)] <- halfmin[is.na(dat)]
  return(dat)
}
