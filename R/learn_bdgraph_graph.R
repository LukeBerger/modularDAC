#' Learn a gene co-expression graph from a data matrix using BDgraph
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param method a character, the BDgraph inference method; 'ggm' for Gaussian graphical models (continuous data) or 'gcgm' for Gaussian copula graphical models (non-Gaussian / mixed data)
#' @param algorithm a character, the sampling algorithm passed to BDgraph::bdgraph; either 'bdmcmc' (birth-death MCMC) or 'rjmcmc' (reversible-jump MCMC)
#' @param iter an integer, the number of MCMC sampling iterations
#' @param burnin an integer, the number of burn-in iterations to discard (must be < iter); defaults to half of iter
#' @param g.prior a numeric between 0 and 1 (or a p x p matrix), the prior probability of edge inclusion
#' @param df.prior a numeric > 2, the degrees of freedom of the G-Wishart prior on the precision matrix
#' @param prob.cut a numeric between 0 and 1, the posterior edge-inclusion probability at or above which an edge is retained (0.5 selects the Bayesian median-probability model)
#' @param pos.cut a numeric, the minimum absolute partial correlation for a positive edge to be retained
#' @param neg.cut a numeric, the minimum absolute partial correlation for a negative edge to be retained
#' @param ... additional arguments passed to BDgraph::bdgraph

#' @return a named list with two elements: 'graph', the learned weighted igraph object whose edge weights are the partial correlations of the retained edges, and 'weights', the symmetric partial correlation matrix derived from the posterior precision matrix

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats cov2cor

#' @export
learn_bdgraph_graph <- function(x,
                                method = c("ggm", "gcgm"),
                                algorithm = c("bdmcmc", "rjmcmc"),
                                iter = 5000,
                                burnin = floor(iter / 2),
                                g.prior = 0.5,
                                df.prior = 3,
                                prob.cut = 0.5,
                                pos.cut = 0,
                                neg.cut = 0,
                                ...){
  if (!requireNamespace("BDgraph", quietly = TRUE)) {
    stop("Package BDgraph is required. Install with: install.packages('BDgraph')", call. = FALSE)
  }

  # handle arguments
  method <- match.arg(method)
  algorithm <- match.arg(algorithm)

  # BDgraph expects an n samples x p features matrix
  t.x <- t(x)

  # fit the Bayesian graphical model
  bd.fit <- suppressMessages(
    BDgraph::bdgraph(data = t.x,
                     method = method,
                     algorithm = algorithm,
                     iter = iter,
                     burnin = burnin,
                     g.prior = g.prior,
                     df.prior = df.prior,
                     ...)
  )

  # posterior edge-inclusion probabilities (upper triangular); mirror to a full
  # symmetric matrix and name by feature
  p.links <- BDgraph::plinks(bd.fit, round = 15)
  p.links[lower.tri(p.links)] <- t(p.links)[lower.tri(p.links)]
  dimnames(p.links) <- list(rownames(x), rownames(x))

  # partial correlations from the posterior mean precision matrix K_hat:
  # pcor_ij = -K_ij / sqrt(K_ii * K_jj), i.e. the negative of cov2cor(K) off the
  # diagonal, with a unit diagonal
  K <- bd.fit$K_hat
  if (is.null(K)) {
    stop("BDgraph did not return a precision matrix ('K_hat'); partial correlations are only available for method 'ggm' or 'gcgm'.", call. = FALSE)
  }
  pcor.avg <- -cov2cor(K)
  diag(pcor.avg) <- 1
  dimnames(pcor.avg) <- list(rownames(x), rownames(x))

  # keep the partial correlation as the edge weight only for edges whose
  # posterior inclusion probability clears prob.cut
  adj.mat <- pcor.avg
  adj.mat[p.links < prob.cut] <- 0

  # zero out edges whose effect size is below the positive/negative cutoffs
  adj.mat[pcor.avg > 0 & abs(pcor.avg) <= abs(pos.cut)] <- 0
  adj.mat[pcor.avg < 0 & abs(pcor.avg) <= abs(neg.cut)] <- 0

  # diagonal to 0 before creating igraph object
  diag(adj.mat) <- 0

  # convert weighted adjacency matrix to igraph object
  g <- igraph::graph_from_adjacency_matrix(adj.mat,
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
