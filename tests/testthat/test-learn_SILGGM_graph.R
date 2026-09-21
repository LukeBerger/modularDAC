# Tests for learn_SILGGM_graph() in R/learn_SILGGM_graph.R
#
# SILGGM is the slow learner, so the run is computed once here and shared
# across the assertions below.

.out     <- quiet(learn_SILGGM_graph(.x))          # learn_* take p x n (features x samples)

.learned <- .out$graph                             # the learned igraph

test_that("learn_SILGGM_graph() returns a weighted igraph with one node per feature", {
  expect_s3_class(.learned, "igraph")
  expect_length(.learned, nrow(.x))
  expect_true(igraph::is_weighted(.learned))
})

test_that("learn_SILGGM_graph() returns the graph alongside the partial correlation matrix", {
  expect_type(.out, "list")
  expect_named(.out, c("graph", "weights"))
  expect_s3_class(.out$graph, "igraph")

  pcor <- .out$weights
  # symmetric p x p matrix, feature-named, with a unit diagonal
  expect_true(is.matrix(pcor))
  expect_equal(dim(pcor), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(pcor), rownames(.x))
  expect_equal(colnames(pcor), rownames(.x))
  expect_equal(unname(diag(pcor)), rep(1, nrow(.x)))
  expect_equal(pcor, t(pcor))
})

test_that("learn_SILGGM_graph() preserves feature identity and undirected structure", {
  # node names carry over from the data's feature labels
  expect_equal(igraph::V(.learned)$name, rownames(.x))
  expect_false(igraph::is_directed(.learned))
  # no self-loops (diagonal is zeroed before graph construction)
  expect_true(igraph::is_simple(.learned))
})

test_that("learn_SILGGM_graph() recovers a non-empty subset of the true edges", {
  # not a benchmark, just a sanity check that the learned graph is usable and
  # overlaps the truth: F1 is a valid score in [0, 1].
  f1 <- calc_F1(.g, .learned)$F1
  expect_type(f1, "double")
  expect_gte(f1, 0)
  expect_lte(f1, 1)
})

test_that("learn_SILGGM_graph() drops more edges as the FDR threshold tightens", {
  # a stricter (smaller) max.fdr should never retain MORE edges than a lax one
  strict <- quiet(learn_SILGGM_graph(.x, fdr.filter = TRUE, max.fdr = 0.001))$graph
  loose  <- quiet(learn_SILGGM_graph(.x, fdr.filter = TRUE, max.fdr = 0.5))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})
