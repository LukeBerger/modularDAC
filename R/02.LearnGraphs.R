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

# helpers to simpleSILGGM graph taken from RSCGGM
# so that i can run it without loading the full package

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

# Shared percentile cut for every learn_*_graph that exposes percentile.threshold
# (ARACNE, WGCNA, CLR, GENIE3), so the four cannot drift apart again.
#
# The quantile is taken over the UPPER TRIANGLE. A whole-matrix quantile -- which
# learn_ARACNE_graph used until this was fixed -- also counts the diagonal and
# both copies of every pair, so the fraction of EDGES kept is not the fraction
# requested.
#
# Sparse weight matrices need care beyond that. ARACNE's DPI step zeroes most
# entries, so when more than `percentile.threshold` of pairs are structurally
# zero the quantile lands ON the zero mass: the cut degenerates to "keep every
# non-zero entry" and different percentiles silently return the SAME graph (this
# is why V.2.0's ARACNE top-2.5% and top-5% lines were identical). That is a real
# property of percentile-thresholding a sparse matrix, not something a quantile
# can fix, so it is surfaced as a warning rather than hidden. Use an explicit
# threshold, or a rank-based edge budget, when it fires.
.percentile_threshold <- function(mat, percentile.threshold, what = "edges") {
  v <- .upper_tri_vec(mat)
  thr <- stats::quantile(v, percentile.threshold, names = FALSE)
  if (thr <= 0) {
    n.pos <- sum(v > 0)
    warning(sprintf(
      "percentile.threshold = %g falls on the zero mass of the %s matrix (%.1f%% of pairs are zero); the cut keeps all %d non-zero %s regardless of the percentile requested.",
      percentile.threshold, what, 100 * mean(v <= 0), n.pos, what), call. = FALSE)
  }
  thr
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

#' Learn a gene co-expression graph from a data matrix using the CLR algorithm
#'
#' CLR (Context Likelihood of Relatedness) is a mutual-information network method
#' from the same family as ARACNE, and shares its first step: build the pairwise
#' mutual information matrix. Where ARACNE prunes indirect edges with the data
#' processing inequality (a hard, triplet-wise deletion), CLR instead RESCORES
#' every pair against the background distribution of its two endpoints -- each
#' edge is converted to a joint z-score of how surprising its MI is relative to
#' all other MIs involving feature i and feature j. Edges are therefore
#' down-weighted rather than deleted, so unlike ARACNE the returned weight matrix
#' has no structural zeros and a percentile threshold cannot degenerate.
#'
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param percentile.threshold a numeric, the top percentile of edges by strength will be kept in the final graph
#' @param clr.threshold a numeric, a fixed CLR score threshold that overrides the percentile, edges bellow this strength will be removed
#' @param skip.diagonal an integer passed to \code{minet::clr}; 1 (the default) excludes self-MI from each feature's background distribution, 0 includes it

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the CLR score matrix

#' @importFrom igraph graph_from_adjacency_matrix

#' @export
learn_CLR_graph <- function(x, percentile.threshold = 0.95, clr.threshold = NULL,
                            skip.diagonal = 1) {
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
  # is not the requested one.
  if(is.null(clr.threshold)){
    clr.threshold <- .percentile_threshold(clr.mat, percentile.threshold, "CLR score")
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

#' Learn a gene co-expression graph from a data matrix using the ARACNE algorithm
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param percentile.threshold a numeric, the top percentile of edges by strength will be kept in the final graph
#' @param mim.threshold a numeric, the a fixed mutual information threshold that overrides the percentile, edges bellow this strength will be removed
#' @param eps a numeric, the data processing inequality threshold used by ARACNE to remove indirect edges

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the mutual information (ARACNE) matrix

#' @importFrom igraph graph_from_adjacency_matrix

#' @export
learn_ARACNE_graph <- function(x, eps=0, percentile.threshold = 0.95, mim.threshold = NULL) {
  if (!requireNamespace("minet", quietly = TRUE)) {
    stop("Package minet is required. Install with: install.packages('minet')", call. = FALSE)
  }

  # minet expects an n samples x p features matrix
  t.x <- t(x)

  # build mutual information matrix
  mim <- minet::build.mim(dataset = t.x)

  # apply ARACNE algorithm
  aracne.mat <- minet::aracne(mim, eps=eps)

  # if no fixed mim.threshold is provided, keep only the 'percentile.threshold' percentile of edges
  if(is.null(mim.threshold)){
    mim.threshold <- .percentile_threshold(aracne.mat, percentile.threshold, "ARACNE mutual information")
  }

  # threshold to select edges based on a minimum mutal information threshold
  weighted.mim <- aracne.mat * (aracne.mat > mim.threshold)

  # remove looped edges
  diag(weighted.mim) <- 0

  # build the weighted graph and return it alongside the full mutual information matrix
  g <- igraph::graph_from_adjacency_matrix(weighted.mim,
                                           mode = "undirected",
                                           weighted = TRUE)
  return(
    list(
      graph   = g,
      weights = aracne.mat
    )
  )
}

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
