#' Learn a gene co-expression graph from a data matrix using GENIE3
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param regulators a character or integer vector, the subset of features (by name or row index) treated as candidate regulators; if NULL all features are considered potential regulators
#' @param targets a character or integer vector, the subset of features (by name or row index) whose regulators are inferred; if NULL all features are treated as targets
#' @param tree.method a character, the tree-based ensemble method used to compute the variable importances; either 'RF' (Random Forests) or 'ET' (Extra-Trees)
#' @param K the number of candidate regulators randomly considered at each tree node; one of 'sqrt' (the square root of the number of regulators), 'all' (all regulators), or a positive integer
#' @param n.trees an integer, the number of trees grown in each feature's ensemble
#' @param symmetrize a character, how the directed GENIE3 importance matrix is collapsed into a symmetric one before building the undirected graph; one of 'max', 'mean', or 'min' of the two link directions
#' @param percentile.threshold a numeric, the top percentile of edges by importance weight will be kept in the final graph
#' @param weight.threshold a numeric, a fixed importance-weight threshold that overrides the percentile; edges below this strength will be removed
#' @param n.cores an integer, the number of cores GENIE3 uses to parallelise the computation

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the symmetric GENIE3 importance (regulatory-link) matrix

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats quantile

#' @export
learn_GENIE3_graph <- function(x,
                               regulators = NULL,
                               targets = NULL,
                               tree.method = c("RF", "ET"),
                               K = "sqrt",
                               n.trees = 1000,
                               symmetrize = c("max", "mean", "min"),
                               percentile.threshold = 0.95,
                               weight.threshold = NULL,
                               n.cores = 1){
  if (!requireNamespace("GENIE3", quietly = TRUE)) {
    stop("Package GENIE3 is required. Install with: BiocManager::install('GENIE3')", call. = FALSE)
  }

  # handle arguments
  tree.method <- match.arg(tree.method)
  symmetrize  <- match.arg(symmetrize)

  # GENIE3 already expects a p features (rows) x n samples (columns) matrix,
  # matching modularDAC's native orientation, so x is passed through untransposed
  link.matrix <- GENIE3::GENIE3(exprMatrix = x,
                                regulators = regulators,
                                targets = targets,
                                treeMethod = tree.method,
                                K = K,
                                nTrees = n.trees,
                                nCores = n.cores,
                                returnMatrix = TRUE,
                                verbose = FALSE)

  # GENIE3 returns a directed importance matrix (rows = regulators, columns =
  # targets) that may span only the requested regulators/targets. Re-embed it in
  # a full p x p feature-named matrix so every feature pair is represented.
  imp.matrix <- matrix(0, nrow = nrow(x), ncol = nrow(x),
                       dimnames = list(rownames(x), rownames(x)))
  imp.matrix[rownames(link.matrix), colnames(link.matrix)] <- link.matrix

  # collapse the two directed link weights of each pair into one undirected weight
  sym.matrix <- switch(symmetrize,
                       max  = pmax(imp.matrix, t(imp.matrix)),
                       min  = pmin(imp.matrix, t(imp.matrix)),
                       mean = (imp.matrix + t(imp.matrix)) / 2)
  dimnames(sym.matrix) <- list(rownames(x), rownames(x))

  # threshold to select edges, but keep the importance values as edge weights.
  # if no fixed weight.threshold is provided, keep only the 'percentile.threshold'
  # top percentile of edges by importance strength (mirrors learn_WGCNA_graph)
  if(is.null(weight.threshold)){
    weight.threshold <- .percentile_threshold(sym.matrix, percentile.threshold, "GENIE3 importance")
  }
  weighted.imp <- sym.matrix * (sym.matrix > weight.threshold)

  # remove looped edges
  diag(weighted.imp) <- 0

  # build the weighted graph and return it alongside the full importance matrix
  g <- igraph::graph_from_adjacency_matrix(weighted.imp,
                                           mode = "undirected",
                                           weighted = TRUE)
  return(
    list(
      graph   = g,
      weights = sym.matrix
    )
  )
}
