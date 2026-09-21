# Tests for the shared .percentile_threshold() edge cut in
# R/percentile_threshold.R

test_that(".percentile_threshold() ignores the diagonal and the lower triangle", {
  m <- matrix(0, 4, 4)
  m[upper.tri(m)] <- 1:6
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 100                            # a whole-matrix quantile would see these
  expect_equal(.percentile_threshold(m, 0.5), stats::quantile(1:6, 0.5, names = FALSE))
})

test_that(".percentile_threshold(nonzero.only = TRUE) excludes the zero mass", {
  m <- matrix(0, 10, 10)
  m[upper.tri(m)] <- c(rep(0, 40), 1:5)     # 45 pairs, 5 of them non-zero
  m[lower.tri(m)] <- t(m)[lower.tri(m)]

  # all pairs: 89% are zero, so any percentile below that lands on the zero mass
  expect_warning(thr.all <- .percentile_threshold(m, 0.5), "falls on the zero mass")
  expect_equal(thr.all, 0)

  # non-zero only: ranks within the 5 survivors, so the cut still responds
  thr.nz <- .percentile_threshold(m, 0.5, nonzero.only = TRUE)
  expect_equal(thr.nz, stats::quantile(1:5, 0.5, names = FALSE))
  expect_gt(thr.nz, 0)
})

test_that(".percentile_threshold(nonzero.only = TRUE) is monotone at any sparsity", {
  set.seed(42)
  m <- matrix(0, 30, 30)
  v <- stats::runif(435)
  v[v < 0.85] <- 0                          # ~85% structural zeros
  m[upper.tri(m)] <- v
  m[lower.tri(m)] <- t(m)[lower.tri(m)]

  ut <- m[upper.tri(m)]
  kept <- vapply(c(0.1, 0.5, 0.9), function(p) {
    sum(ut > .percentile_threshold(m, p, nonzero.only = TRUE))
  }, numeric(1))
  expect_equal(kept, sort(kept, decreasing = TRUE))
  expect_equal(anyDuplicated(kept), 0L)
})

test_that(".percentile_threshold(nonzero.only = TRUE) warns and keeps nothing on an all-zero matrix", {
  m <- matrix(0, 5, 5)
  expect_warning(thr <- .percentile_threshold(m, 0.95, nonzero.only = TRUE),
                 "every pair")
  expect_equal(thr, 0)
  expect_equal(sum(m[upper.tri(m)] > thr), 0L)   # an empty graph, not an error
})
