# Shared percentile cut for the learners that expose percentile.threshold
# (ARACNE, WGCNA, CLR, GENIE3).
# The quantile is taken over the upper triangle, so the fraction of node pairs
# kept is the fraction requested; a whole-matrix quantile would also count the
# diagonal and both copies of every pair.
# `nonzero.only` selects which population the percentile ranks within, which
# matters for any learner whose weight matrix carries structural zeros:
#   FALSE -- all p(p-1)/2 pairs, making the percentile a density target ("keep
#     the top X% of node pairs"). Appropriate when every pair is a live
#     candidate, as for WGCNA's soft-power adjacency and GENIE3's importances.
#     On a sparse matrix it degenerates: once more than 1 - percentile.threshold
#     of pairs are zero the quantile lands on the zero mass, the cut collapses
#     to "keep every non-zero entry", and different percentiles return the same
#     graph. That case warns.
#   TRUE -- the non-zero pairs only, making the percentile a relative-strength
#     filter ("drop the weakest X% of the pairs still in contention").
#     Appropriate for ARACNE and CLR, whose zeros are a deletion the algorithm
#     has already decided on (DPI pruning, or CLR clamping a negative background
#     z-score to 0) rather than a measurement, so those pairs are no longer
#     candidates. Ranking within the survivors keeps the cut monotone in
#     percentile.threshold at any sparsity.
# The two modes are on different scales -- the top 5% of survivors is a far
# smaller edge count than the top 5% of pairs -- so a percentile tuned under one
# does not carry over to the other. Neither is comparable across graph sizes: a
# percentile fixes a fraction of node pairs (~p^2) while the number of true
# edges scales with ~p.
.percentile_threshold <- function(mat, percentile.threshold, what = "edges",
                                  nonzero.only = FALSE) {
  v <- .upper_tri_vec(mat)

  # rank within the pairs still in contention: there is no zero mass for the
  # quantile to land on, so the cut stays monotone in percentile.threshold
  if (nonzero.only) {
    pos <- v[v > 0]
    if (length(pos) == 0L) {
      warning(sprintf(
        "every pair of the %s matrix is zero, so no edges can be kept.",
        what), call. = FALSE)
      return(0)
    }
    return(stats::quantile(pos, percentile.threshold, names = FALSE))
  }

  thr <- stats::quantile(v, percentile.threshold, names = FALSE)
  if (thr <= 0) {
    n.pos <- sum(v > 0)
    warning(sprintf(
      "percentile.threshold = %g falls on the zero mass of the %s matrix (%.1f%% of pairs are zero); the cut keeps all %d non-zero %s regardless of the percentile requested. Pass nonzero.only = TRUE to rank within the non-zero entries instead.",
      percentile.threshold, what, 100 * mean(v <= 0), n.pos, what), call. = FALSE)
  }
  thr
}
