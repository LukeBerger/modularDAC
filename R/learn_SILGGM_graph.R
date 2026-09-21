# allow for dplyr piping
utils::globalVariables(".")

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

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the symmetric partial correlation matrix (or, for the global-FDR methods 'GFC_SL'/'GFC_L', the symmetric test-statistic matrix)

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom dplyr  %>%

#' @export
learn_SILGGM_graph <- function(x,
                               method  = "B_NW_SL",
                               global = TRUE,
                               alpha = 0.05,
                               fdr.filter = TRUE,
                               max.fdr = 0.05,
                               pos.cut = 0,
                               neg.cut = 0,
                               ...){
  if (!requireNamespace("SILGGM", quietly = TRUE)) {
    stop("Package SILGGM is required. Install with: install.packages('SILGGM')", call. = FALSE)
  }

  # SILGGM expects an n samples x p features matrix
  t.x <- t(x)

  # run silggm
  silggm.output <- suppressMessages(
    SILGGM::SILGGM(t.x,
                   method = method,
                   alpha = alpha,
                   global = global,
                   ...)
  )
  # build the returned weight matrix (pcor.avg) and the thresholded edge matrix
  # (adj.mat) in a method-dependent way
  if(method %in% c("GFC_SL", "GFC_L")){
    # global-FDR methods report no partial correlations: weight edges by the test
    # statistic and select them with SILGGM's own FDR-controlled global decision
    pcor.avg <- silggm.output$T_stat
    dimnames(pcor.avg) <- list(rownames(x), rownames(x))
    adj.mat <- pcor.avg
    if(fdr.filter){
      # global_decision holds one binary edge set per supplied alpha level
      lvl <- which(alpha == max.fdr)[1]
      if(is.na(lvl)) lvl <- 1L
      adj.mat <- adj.mat * silggm.output$global_decision[[lvl]]
    }
  } else {
    # partial-correlation methods: extract partial correlations and edge z-scores
    if(method == "B_NW_SL"){
      # this method reports a partial-correlation z-score directly
      partial.cor <- .upper_tri_vec(silggm.output$partialCor)
      z.score.partial.cor <- .upper_tri_vec(silggm.output$z_score_partialCor)
    } else if(method %in% c("D-S_NW_SL", "D-S_GL")){
      # these report partialCor but no partial-correlation z-score, so derive it
      partial.cor <- .upper_tri_vec(silggm.output$partialCor)
      z.score.partial.cor <- .pcor_zscore(partial.cor, ncol(x))
    } else {
      stop("method '", method, "' is not supported by learn_SILGGM_graph.")
    }

    # estimate the adjusted p-value of each edge
    qval.pcor <- .pvalue(z.score.partial.cor) %>%
      .upper_tri_to_matrix(., variable_names = rownames(x), diagl = 1) %>%
      .matrix_p_adjust(.)

    # build partial correlation adjacency matrix
    pcor.avg <- partial.cor %>%
      round(.,2) %>%
      .upper_tri_to_matrix(., variable_names = rownames(x), diagl = 1)

    # zero out edges whose effect size is below the positive/negative cutoffs
    adj.mat <- pcor.avg
    adj.mat[pcor.avg > 0 & abs(pcor.avg) <= abs(pos.cut)] <- 0
    adj.mat[pcor.avg < 0 & abs(pcor.avg) <= abs(neg.cut)] <- 0

    # further filtering by significance
    if(fdr.filter){
      adj.mat[qval.pcor>=max.fdr] <- 0
    }
  }

  # diagonal to 0 before creating igraph object
  diag(adj.mat) <- 0
  adj.mat[is.na(adj.mat)] <- 0

  # convert weighted adjacency matrix to igraph object
  g  <- igraph::graph_from_adjacency_matrix(adj.mat,
                                            mode = "undirected",
                                            weighted = TRUE)

  # return the learned graph and the partial correlation matrix
  return(
    list(
      graph   = g,
      weights = pcor.avg
    )
  )
}

# Matrix helpers for learn_SILGGM_graph, adapted from the RSCGGM package.

#' Reconstruct the symmetric matrix from upper triangular vector
#'
#' @param upper_tri_values a numeric vector of the upper triangle of the matrix
#' @param variable_names a character vector of row and column names for the output matrix
#' @param diagl a numeric scalar or vector of length p, the diagonal value(s) of the output matrix

#' @return a symmetric matrix

#' @keywords internal
.upper_tri_to_matrix <- function(upper_tri_values,
                                variable_names =NULL,
                                diagl=1){
  p <- (1 + sqrt(1 + 8 * length(upper_tri_values))) / 2
  if( (length(diagl)>1) && (length(diagl) != p) ) stop("invalid diagonal!")

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
  std_new <- sqrt((1 - pcor^2)^2 / n)
  return(pcor/std_new)
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
