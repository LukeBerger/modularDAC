#' Simulate gene expression data for n.samples samples based on the covariance structure implied by a graph
#' @param g an igraph object defining the dependency structure between features
#' @param n.samples an integer, the number of samples to generate
#' @param mean.vec a numeric vector of length p, the mean expression value for each feature; defaults to a zero vector if NULL
#' @param b a numeric (> 2), the degrees-of-freedom parameter of the G-Wishart distribution used to sample the precision matrix (passed to BDgraph::rgwish)

#' @return a p x n data matrix

#' @importFrom igraph as_adjacency_matrix V

#' @keywords internal
#' @note Benchmarking helper, kept on this branch only to back the test
#'   suite. It is not part of the user-facing API; the exported version
#'   lives on the dev branch.
sim_graph_data <- function(g, n.samples, mean.vec = NULL, b = 3){
  if (!requireNamespace("BDgraph", quietly = TRUE)) {
    stop("Package BDgraph is required. Install with: install.packages('BDgraph')", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package MASS is required. Install with: install.packages('MASS')", call. = FALSE)
  }
  # if mean.vec is not set, default to a zero vector
  if(is.null(mean.vec)){
    # zero mean for all features
    mean.vec <- rep(0, length(g))
  }
  # derive the covariance matrix implied by the graph structure: sample a
  # precision matrix from the G-Wishart distribution constrained to the graph
  # (forcing zeros at non-edges), then convert it to a correlation matrix
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(g))
  diag(adj.mat) <- 0
  K <- BDgraph::rgwish(adj = adj.mat, b = b, D = diag(nrow(adj.mat)), threshold = 1e-8)
  sigma.mat <- stats::cov2cor(solve(K)) # correlation matrix implied by the graph

  # make n samples with MASS
  sim.data <- t(MASS::mvrnorm(n = n.samples, mu = mean.vec, Sigma =  sigma.mat))

  # add row and column names
  colnames(sim.data) <- paste("Sample", seq_len(ncol(sim.data)), sep = "_")
  rownames(sim.data) <- igraph::V(g)$name

  # return data
  return(sim.data)
}
