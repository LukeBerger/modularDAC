# Tests for learn_CLR_graph() in R/learn_CLR_graph.R

test_that("learn_CLR_graph() returns a weighted igraph alongside the CLR score matrix", {
  skip_if_not_installed("minet")
  out <- quiet(learn_CLR_graph(.x))        # mutual information over features (p x n input)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))

  clr <- out$weights
  expect_true(is.matrix(clr))
  expect_equal(dim(clr), c(nrow(.x), nrow(.x)))
  # feature names must survive so the weights matrix can be thresholded
  # externally (the benchmark's edge-budget rule works off dimnames)
  expect_equal(rownames(clr), rownames(.x))
})

test_that("learn_CLR_graph() keeps fewer edges at a stricter percentile", {
  skip_if_not_installed("minet")
  loose  <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.99))$graph
  expect_lt(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_CLR_graph() ranks the percentile within the non-zero scores", {
  skip_if_not_installed("minet")
  # minet::clr clamps a negative endpoint z-score to 0, so the CLR matrix is
  # ~70% structural zeros and the default cut ranks within the rest
  w   <- quiet(learn_CLR_graph(.x))$weights
  ut  <- w[upper.tri(w)]
  pos <- ut[ut > 0]
  expect_gt(mean(ut <= 0), 0)               # the zero mass is real, not hypothetical

  g <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.90))$graph
  expect_equal(igraph::gsize(g), sum(pos > stats::quantile(pos, 0.90)))
})

test_that("learn_CLR_graph(nonzero.only = FALSE) tracks the requested fraction of PAIRS", {
  skip_if_not_installed("minet")
  # the pre-0.0.0.9015 behaviour: the upper-triangle quantile spans every pair,
  # so the kept count tracks the requested fraction of pairs -- unlike a
  # whole-matrix quantile, which double-counts
  w  <- quiet(learn_CLR_graph(.x))$weights
  ut <- w[upper.tri(w)]
  p  <- 0.99
  skip_if(mean(ut <= 0) >= p, "matrix too sparse for a non-degenerate all-pairs cut")

  g <- quiet(learn_CLR_graph(.x, percentile.threshold = p, nonzero.only = FALSE))$graph
  n.pairs <- nrow(.x) * (nrow(.x) - 1) / 2
  expect_equal(igraph::gsize(g), sum(ut > stats::quantile(ut, p)))
  expect_lt(abs(igraph::gsize(g) - (1 - p) * n.pairs), (1 - p) * n.pairs)
})

test_that("learn_CLR_graph() honours an explicit clr.threshold over the percentile", {
  skip_if_not_installed("minet")
  w <- quiet(learn_CLR_graph(.x))$weights
  cut <- as.numeric(stats::quantile(w[upper.tri(w)], 0.80))
  g <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.99, clr.threshold = cut))$graph
  # the explicit cut (80th pct) must win over the percentile argument (99th)
  expect_equal(igraph::gsize(g), sum(w[upper.tri(w)] > cut))
})
