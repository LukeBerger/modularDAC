#' Wrapper for netUtils::sample_lfr that generates a Lancichinetti-Fortunato-Radicchi (LFR) benchmark graph and labels nodes
#' @param n an integer, the number of nodes in the graph
#' @param tau1 a numeric, the power-law exponent for the degree distribution
#' @param tau2 a numeric, the power-law exponent for the community size distribution
#' @param mu a numeric between 0 and 1, the mixing parameter controlling the fraction of inter-community edges per node
#' @param average.degree a numeric, the average node degree in the graph
#' @param max.degree an integer, the maximum node degree allowed in the graph
#' @param min.community an integer, the minimum number of nodes in each community
#' @param max.community an integer, the maximum number of nodes in each community

#' @return an igraph object

#' @importFrom igraph V

#' @keywords internal
#' @note Benchmarking helper, kept on this branch only to back the test
#'   suite. It is not part of the user-facing API; the exported version
#'   lives on the dev branch.
make_lfr <- function(n = 120,
                     tau1 = 3,
                     tau2 = 2,
                     mu = 0.08,
                     average.degree = 6,
                     max.degree = 10,
                     min.community = 30,
                     max.community = 50){
  if (!requireNamespace("netUtils", quietly = TRUE)) {
    stop("Package netUtils is required. Install with: install.packages('netUtils')", call. = FALSE)
  }
  # generate LFR benchmark graph
  g <- netUtils::sample_lfr(
    n = n,
    tau1 = tau1,
    tau2 = tau2,
    mu = mu,
    average_degree = average.degree,
    max_degree = max.degree,
    min_community = min.community,
    max_community = max.community
  )

  # assign sequential node names and copy community membership to module attribute
  igraph::V(g)$name <- paste0("Node_", igraph::V(g))
  igraph::V(g)$module <- igraph::V(g)$membership

  return(g)
}
