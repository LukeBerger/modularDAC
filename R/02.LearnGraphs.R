### Allow for dplyr piping
utils::globalVariables(".")

######################
### Graph learning ###
######################

### Takes a n x p matrix of data and inputs to silgggm function
### Then generates networks based on partial correlation between nodes

#' Learns a graph from data using SILGGM
#' @param x an  n by p matrix of numeric data
#' @param method Methods for statistical inference with 5 options: "B_NW_SL", "D-S_NW_SL", "D-S_GL", "GFC_SL" and "GFC_L"
#' @param global 	If global = TRUE, the global inference of all gene pairs is performed.
#' @param alpha A user-supplied sequence of pre-sepecified alpha levels for FDR control. The default is alpha = 0.05, 0.1 if no sequence is provided.
#' @param fdr.filter Whether to filter adjacency matrix results based on
#' @param max.fdr The threshold to define an edge based on partial cor q val
#' @param pos.cut Threhold to zero values in adjacency matrix
#' @param neg.cut Threhold to zero values in adjacency matrix
#' @param ... Additional arguments passed to SILGGM function

#' @return an igraph object

#' @importFrom SILGGM SILGGM
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

  #run silggm
  silggm.output <- suppressMessages(
    SILGGM::SILGGM(x,
                   method = method,
                   alpha = alpha,
                   global = global,
                   ...)
  )
  #get partial.cor and zScore based on method
  if(method == "D-S_NW_SL"){
    partial.cor <- .upper_tri_vec(silggm.output$partialCor)
    z.score.partial.cor <- sapply(partial.cor, .pcor_zscore, nrow(x))
  }
  if(method == "B_NW_SL"){
    partial.cor <- .upper_tri_vec(silggm.output$partialCor)
    z.score.partial.cor <- .upper_tri_vec(silggm.output$z_score_partialCor)
  }

  #estimate the adjusted p-value of each edge
  qval.pcor <- sapply(z.score.partial.cor, .pvalue) %>%
    .upper_tri_to_matrix(., variable_names = colnames(x), diagl = 1) %>%
    .matrix_p_adjust(.)

  #make network
  pcor.avg <- partial.cor %>%
    round(.,2) %>%
    .upper_tri_to_matrix(., variable_names = colnames(x), diagl = 1)
  adj.mat <- apply(pcor.avg, c(1,2), .abs_pcor_filter, pos.cut, neg.cut)

  ## further filtering by significance
  if(fdr.filter){
    adj.mat[qval.pcor>=max.fdr] <- 0
  }

  ## diagonal to 0 before creating igraph object
  diag(adj.mat) <- 0

  # to igraph object
  g  <- igraph::graph_from_adjacency_matrix(adj.mat,
                                            mode = "undirected",
                                            weighted = TRUE)

  return(g)
}


## Helpers to simpleSILGGM graph taken from RSCGGM
## So that i can run it without loading the full package


#' Reconstruct the symmetric matrix from upper triangular vector
#'
#' @param upper_tri_values a numeric vector of the upper triangle of the matrix
#' @param variable_names row&column names
#' @param diagl the diagnoal elemenet, a number of a numerci vector

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

  ## Fill the diagonal
  diag(mat) <- diagl

  ## Fill the upper triangular part
  mat[upper.tri(mat, diag = FALSE)] <- upper_tri_values

  # Fill the lower triangular part (mirror the upper triangular part)
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
  ## initialize
  mx_q <- mx_p
  ## adjust upper triangle
  mx_q[upper.tri(mx_q)] <-
    stats::p.adjust(mx_p[upper.tri(mx_p)], method = "BH")
  ## copy to lower triangle
  mx_q[lower.tri(mx_q)] <-
    t(mx_q)[lower.tri(mx_q)]
  ## some checks
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

#' Uses WGCNA to learn a graph
#' @param x a n x p matrix of features
#' @param min.sft an integer, the min sft used by the pickSoftThreshold function
#' @param beta an integer, the power value used by WGCNA::adjacency
#' @param cor.FN a character, the cor.options WGCNA::adjacency
#' @param powers an integer vector, the powers vector used by WGCNA::pickSoftThreshold
#' @param threshold value to binarize adj matrix against

#' @return a igraph object object

#' @importFrom WGCNA pickSoftThreshold adjacency bicor
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
  # Handle arguments
  cor.FN <- match.arg(cor.FN)

  # Correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=x,
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # Check selected power
    beta <- .sft_check(sft)
  }

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=x,
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # Binarize based on threshold
  if(is.null(threshold)){threshold <- quantile(.upper_tri_vec(adj), .999)}
  binary.adj <- (adj > threshold) * 1

  # Remove looped edges
  diag(binary.adj) <- 0

  # Return as Igraph
  return(igraph::graph_from_adjacency_matrix(binary.adj, mode = "undirected"))
}



####################
###Compare Graphs###
####################

#' Calculates F1 score from a true and predicted graph
#' @param g.true an igraph object, the true graph for the purpose of F1 calculation
#' @param g.pred an igraph object, the predicted graph for the purpose of F1 calculation

#' @return a list of F1 score and its components

#' @importFrom igraph V as_edgelist

#' @export
calc_F1 <- function(g.true, g.pred) {
  #check that the graphs have the same nodes
  if (!all(igraph::V(g.true)$name %in% igraph::V(g.pred)$name)) {
    stop("True and Predicted graphs must contain the same nodes")
  }

  #extract edges as character vectors in a consistent format
  edges.true <- apply(igraph::as_edgelist(g.true), 1, function(x) paste(sort(x), collapse = " -- "))
  edges.pred <- apply(igraph::as_edgelist(g.pred), 1, function(x) paste(sort(x), collapse = " -- "))

  #count True Positives, False Positives, False Negatives
  tp <- sum(edges.pred %in% edges.true)
  fp <- sum(!(edges.pred %in% edges.true))
  fn <- sum(!(edges.true %in% edges.pred))

  #calc precision and recall
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))

  #calc F1 Score
  if (precision + recall == 0) {
    f1 <- 0 #to avoid zero divison error in edge case
  } else {
    f1 <- 2 * precision * recall / (precision + recall)
  }

  #return values as list
  return(list(
    F1 = f1,
    Precision = precision,
    Recall = recall,
    TruePos = tp,
    FalsePos = fp,
    FalseNeg = fn
  ))
}


###############
#### Other ####
###############

#' Preforms halfmin imputation to fill NAs in matrix missing values, useful for filling NAs in matrix before graph learning
#' @param dat data matrix of numeric values

#' @return the data matrix with NAs filled

#' @importFrom matrixStats rowMins rowMeans2

#' @keywords internal
.halfmin_impute <- function(dat) {
  halfmin <- matrix(matrixStats::rowMins(dat, na.rm = TRUE) / 2,
                    nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat)
  )
  ## just checking
  stopifnot(all.equal(
    matrixStats::rowMeans2(halfmin),
    matrixStats::rowMins(dat, na.rm = TRUE) / 2
  ))
  dat[is.na(dat)] <- halfmin[is.na(dat)]
  return(dat)
}


###########
## TEST ###
###########

# # make data
# er <- make_modular_graph()
# x <- sim_graph_data(er, n.samples = 100)
#
# # test silggm learning
# s <- learn_SILGGM_graph(t(x))
# calc_F1(er, s)
