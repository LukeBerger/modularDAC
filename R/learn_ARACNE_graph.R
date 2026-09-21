#' Learn a gene co-expression graph from a data matrix using the ARACNE algorithm
#'
#' ARACNE builds the pairwise mutual information matrix, then prunes indirect
#' edges with the data processing inequality (DPI): for every triplet it deletes
#' the weakest of the three MI values outright. Most of the matrix is therefore a
#' structural zero by the time it is thresholded, which is what `nonzero.only`
#' handles -- a percentile spanning all node pairs lands on that zero mass and
#' stops responding to the percentile (see `.percentile_threshold`).
#'
#' Four edge-selection rules are available, in precedence order: `mim.threshold`
#' (an absolute MI cutoff), `p.threshold` (the raw, unadjusted p-value of each MI
#' under the independence null), `percentile.threshold` (keep the top slice by MI
#' strength) and `max.fdr` (the same null test on the Benjamini-Hochberg adjusted
#' p-value). `max.fdr = 0.05` is the default and sits last in the precedence
#' order, so setting any other rule disables it; setting all four to NULL is an
#' error.
#'
#' `max.fdr` is the default because it is the rule Margolin et al. describe, and
#' the only truth-free one that needs no calibration against graph size: a
#' percentile fixes a fraction of node pairs (~p^2) while the number of true
#' edges scales with ~p, so a single percentile cannot be right at two sizes,
#' whereas the MI null adapts to the sample size on its own.
#'
#' Margolin et al. threshold the raw MI matrix and then run DPI over what
#' survives, where this function runs DPI over the complete MI graph and
#' thresholds afterwards. For the significance rules, and for `mim.threshold`,
#' the two orders give the same graph, so no reordering is offered. The p-value
#' is a strictly decreasing function of MI at fixed n, so any rule built on
#' p-values (raw, or BH, which takes a prefix of the sorted p-values) is exactly
#' `MI >= c` for some c. DPI prunes (i,j) only via a witness k with
#' `min(MI(i,k), MI(j,k)) > MI(i,j) + eps`, so if (i,j) clears c then both
#' witness edges clear c as well and survive the mask: masking can never remove a
#' witness that would have mattered.
#'
#' Note what a significance threshold does and does not do. It rules out pairs
#' whose MI is explicable as sampling noise. It does not rule out indirect edges:
#' a pair joined only by a path through other features is genuinely dependent, so
#' its true MI is non-zero and it becomes significant at large enough n. DPI,
#' applied afterwards, is what removes those, and being a three-node rule it
#' cannot see paths of length three or more.
#'
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param percentile.threshold a numeric, or NULL (the default, leaving this rule off); the top percentile of edges by strength will be kept in the final graph, ranked within the pairs that survive DPI unless `nonzero.only = FALSE`. Supplying it overrides the default `max.fdr` rule. Note that it is not calibrated across graph sizes, and that a value tuned under one `nonzero.only` setting does not carry over to the other
#' @param mim.threshold a numeric, the a fixed mutual information threshold that overrides the percentile, edges bellow this strength will be removed
#' @param eps a numeric, the data processing inequality threshold used by ARACNE to remove indirect edges
#' @param nonzero.only a logical; if TRUE (the default) `percentile.threshold` is a percentile of the non-zero entries of the ARACNE matrix, the pairs still in contention after DPI pruning, so the cut drops the weakest X% of surviving edges. If FALSE the percentile spans all p(p-1)/2 node pairs: because DPI zeroes most entries, that cut collapses to "keep every non-zero entry" (and warns) whenever more than 1 - percentile.threshold of pairs are zero, making different percentiles return the same graph. The two settings are on different scales -- the top 5% of survivors is far fewer edges than the top 5% of pairs -- so a percentile tuned under one does not carry over

#' @param max.fdr a numeric between 0 and 1, or NULL; the default edge rule (0.05). An edge is kept when the Benjamini-Hochberg adjusted p-value of its mutual information, tested against the independence null by \code{.mi_pvalue}, falls below it. Mirrors `learn_SILGGM_graph`'s `max.fdr`. It sits last in the precedence order, so it applies only when no other rule is set; pass `max.fdr = NULL` together with another rule, or simply set that other rule, to turn it off
#' @param p.threshold a numeric between 0 and 1, or NULL (the default); as `max.fdr` but testing the raw, unadjusted per-edge p-value, the form ARACNE's authors use. With p(p-1)/2 tests this needs a very small value to control false positives (the Bonferroni level is `0.05 / (p * (p - 1) / 2)`). Supplying it overrides `max.fdr`

#' @return a named list with two elements: 'graph', the learned weighted igraph object, and 'weights', the mutual information (ARACNE) matrix

#' @references Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G,
#'   Dalla Favera R, Califano A (2006). ARACNE: an algorithm for the
#'   reconstruction of gene regulatory networks in a mammalian cellular context.
#'   BMC Bioinformatics 7(Suppl 1):S7.

#' @importFrom igraph graph_from_adjacency_matrix

#' @export
learn_ARACNE_graph <- function(x, eps=0, percentile.threshold = NULL, mim.threshold = NULL,
                               nonzero.only = TRUE, max.fdr = 0.05, p.threshold = NULL) {
  if (!requireNamespace("minet", quietly = TRUE)) {
    stop("Package minet is required. Install with: install.packages('minet')", call. = FALSE)
  }

  # Resolve the edge-selection rule. max.fdr is the default and sits last in
  # precedence: the other three are NULL unless the caller sets them, so setting
  # any one of them overrides the default.
  rule <- if (!is.null(mim.threshold))            "absolute"
          else if (!is.null(p.threshold))         "pvalue"
          else if (!is.null(percentile.threshold)) "percentile"
          else if (!is.null(max.fdr))             "fdr"
          else stop("no edge-selection rule: supply one of max.fdr, p.threshold, ",
                    "percentile.threshold or mim.threshold.", call. = FALSE)

  # minet expects an n samples x p features matrix
  t.x <- t(x)

  # build mutual information matrix
  mim <- minet::build.mim(dataset = t.x)

  # apply ARACNE algorithm
  aracne.mat <- minet::aracne(mim, eps = eps)

  # Select edges. The significance rules test the raw MI, never the DPI output:
  # DPI zeroes entries as a pruning decision, so a pruned pair would have a
  # p-value of 1 however dependent it is, and testing the pruned matrix would
  # re-apply DPI's verdict instead of asking about independence. Intersecting
  # that mask with the DPI output gives the same graph as re-running DPI on the
  # masked MI; see the note on ordering in this function's documentation.
  keep <- switch(rule,
    absolute = aracne.mat > mim.threshold,
    fdr      = .matrix_p_adjust(.mi_pvalue(mim, ncol(x))) < max.fdr,
    pvalue   = .mi_pvalue(mim, ncol(x)) < p.threshold,
    # keep only the 'percentile.threshold' percentile of edges -- by default
    # ranked within the DPI survivors rather than all node pairs, see
    # .percentile_threshold()
    percentile = aracne.mat > .percentile_threshold(aracne.mat, percentile.threshold,
                                                    "ARACNE mutual information",
                                                    nonzero.only = nonzero.only)
  )
  weighted.mim <- aracne.mat * keep

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

#' Two-sided p-value for a mutual information estimate under the independence null
#'
#' `minet::build.mim()`'s default estimator returns `MI = -0.5 * log(1 - rho^2)`
#' for the Spearman rank correlation `rho`, a monotone function of a rank
#' statistic. Under the null of independent, continuously distributed features
#' the null is therefore distribution-free and closed form: `rho` has variance
#' `1 / (n - 1)`, so `2 * (n - 1) * MI ~ chi-squared(1)`. This is the analytic
#' form of the permutation threshold described by Margolin et al. (2006), and it
#' agrees with a permutation null to within a few percent across the tail.
#'
#' The result is tied to the Spearman estimator `learn_ARACNE_graph` uses. A
#' discretised or kernel MI estimator would need a permutation null instead,
#' since its bias at fixed bin width is not captured by this reference
#' distribution.
#'
#' @param mim a p x p mutual information matrix from `minet::build.mim()`
#' @param n.samples an integer, the number of samples the MI was estimated from
#'
#' @return a p x p matrix of p-values for the null hypothesis that the pair is independent
#'
#' @references Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G,
#'   Dalla Favera R, Califano A (2006). ARACNE: an algorithm for the
#'   reconstruction of gene regulatory networks in a mammalian cellular context.
#'   BMC Bioinformatics 7(Suppl 1):S7.

#' @importFrom stats pchisq
#' @keywords internal
.mi_pvalue <- function(mim, n.samples) {
  if (n.samples < 3L) {
    stop("at least 3 samples are needed to test a mutual information estimate.", call. = FALSE)
  }
  pv <- stats::pchisq(2 * (n.samples - 1) * mim, df = 1, lower.tail = FALSE)
  dim(pv) <- dim(mim)
  dimnames(pv) <- dimnames(mim)
  pv
}
