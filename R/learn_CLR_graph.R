#' Learn a gene co-expression graph from a data matrix using the CLR algorithm
#'
#' CLR (Context Likelihood of Relatedness) is a mutual-information network
#' method from the same family as ARACNE and shares its first step: build the
#' pairwise mutual information matrix. Where ARACNE prunes indirect edges with
#' the data processing inequality (a hard, triplet-wise deletion), CLR rescores
#' every pair against the background distribution of its two endpoints, so each
#' edge becomes a joint z-score of how surprising its MI is relative to all other
#' MIs involving feature i and feature j. Edges are down-weighted rather than
#' deleted outright, and CLR keeps many more pairs than ARACNE does.
#'
#' The CLR matrix is still sparse: `minet::clr` clamps a negative endpoint
#' z-score to 0, so every pair whose MI falls below the background mean of both
#' its endpoints is a structural zero. A percentile spanning all node pairs
#' therefore degenerates on CLR as it does on ARACNE, which is why
#' `nonzero.only` defaults to TRUE here (see `.percentile_threshold`).
#'
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param percentile.threshold a numeric, the top percentile of edges by strength will be kept in the final graph; ranked within the non-zero CLR scores unless `nonzero.only = FALSE`
#' @param clr.threshold a numeric, a fixed CLR score threshold that overrides the percentile, edges bellow this strength will be removed
#' @param skip.diagonal an integer passed to \code{minet::clr}; 1 (the default) excludes self-MI from each feature's background distribution, 0 includes it
#' @param nonzero.only a logical; if TRUE (the default) `percentile.threshold` is a percentile of the non-zero CLR scores, the pairs CLR left in contention, so the cut drops the weakest X% of those. If FALSE the percentile spans all p(p-1)/2 node pairs, which collapses to "keep every non-zero score" (and warns) whenever more than 1 - percentile.threshold of pairs are zero. The two settings are on different scales -- the top 5% of non-zero scores is far fewer edges than the top 5% of pairs -- so a percentile tuned under one does not carry over

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the CLR score matrix

#' @references Faith JJ, Hayete B, Thaden JT, Mogno I, Wierzbowski J, Cottarel G,
#'   Kasif S, Collins JJ, Gardner TS (2007). Large-scale mapping and validation
#'   of Escherichia coli transcriptional regulation from a compendium of
#'   expression profiles. PLoS Biology 5(1):e8.

#' @importFrom igraph graph_from_adjacency_matrix

#' @export
learn_CLR_graph <- function(x, percentile.threshold = 0.95, clr.threshold = NULL,
                            skip.diagonal = 1, nonzero.only = TRUE) {
  if (!requireNamespace("minet", quietly = TRUE)) {
    stop("Package minet is required. Install with: install.packages('minet')", call. = FALSE)
  }

  # minet expects an n samples x p features matrix
  t.x <- t(x)

  # build mutual information matrix (same first step as learn_ARACNE_graph)
  mim <- minet::build.mim(dataset = t.x)

  # apply CLR: rescore each pair against the MI background of both endpoints
  clr.mat <- minet::clr(mim, skipDiagonal = skip.diagonal)
  dimnames(clr.mat) <- dimnames(mim)

  # if no fixed clr.threshold is provided, keep only the 'percentile.threshold'
  # percentile of edges. The quantile is taken over the UPPER TRIANGLE (as
  # learn_WGCNA_graph does) rather than the whole matrix: a full-matrix quantile
  # also counts the diagonal and both copies of every pair, so the kept fraction
  # is not the requested one. By default it ranks within the non-zero scores
  # rather than all node pairs, see .percentile_threshold()
  if(is.null(clr.threshold)){
    clr.threshold <- .percentile_threshold(clr.mat, percentile.threshold, "CLR score",
                                           nonzero.only = nonzero.only)
  }

  # threshold to select edges based on a minimum CLR score
  weighted.clr <- clr.mat * (clr.mat > clr.threshold)

  # remove looped edges
  diag(weighted.clr) <- 0

  # build the weighted graph and return it alongside the full CLR score matrix
  g <- igraph::graph_from_adjacency_matrix(weighted.clr,
                                           mode = "undirected",
                                           weighted = TRUE)
  return(
    list(
      graph   = g,
      weights = clr.mat
    )
  )
}
