#' Learn a gene co-expression graph from a data matrix using WGCNA
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param min.sft a numeric between 0 and 1, the minimum R-squared threshold for soft-thresholding power selection
#' @param beta an integer, the soft-thresholding power for WGCNA::adjacency; if NULL it is selected automatically via pickSoftThreshold
#' @param cor.FN a character, the correlation function to use in WGCNA::adjacency; either 'bicor' or 'cor'
#' @param powers an integer vector, candidate soft-thresholding powers evaluated by WGCNA::pickSoftThreshold
#' @param percentile.threshold a numeric, the top percentile of edges by adjacency strength will be kept in the final graph
#' @param adj.threshold a numeric, a fixed adjacency threshold that overrides the percentile; edges below this strength will be removed

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the WGCNA co-expression adjacency matrix

#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom stats quantile

#' @export
learn_WGCNA_graph <- function(x,
                              min.sft=0.85,
                              beta=NULL,
                              cor.FN=c("bicor", "cor"),
                              powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                              percentile.threshold = 0.95,
                              adj.threshold = NULL
) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package WGCNA is required. Install with: install.packages('WGCNA')", call. = FALSE)
  }
  # handle arguments
  cor.FN <- match.arg(cor.FN)

  # correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # WGCNA expects an n samples x p features matrix
  t.x <- t(x)

  # pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=t.x,
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # check selected power
    beta <- .sft_check(sft)
  }

  # construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=t.x,
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # threshold to select edges, but keep the adjacency values as edge weights.
  # if no fixed adj.threshold is provided, keep only the 'percentile.threshold'
  # top percentile of edges by adjacency strength (mirrors learn_ARACNE_graph)
  if(is.null(adj.threshold)){
    adj.threshold <- .percentile_threshold(adj, percentile.threshold, "WGCNA adjacency")
  }
  weighted.adj <- adj * (adj > adj.threshold)

  # remove looped edges
  diag(weighted.adj) <- 0

  # build the weighted graph and return it alongside the full adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(weighted.adj,
                                           mode = "undirected",
                                           weighted = TRUE)
  return(
    list(
      graph   = g,
      weights = adj
    )
  )
}
