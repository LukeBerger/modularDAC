# Tests for learn_GENIE3_graph() in R/learn_GENIE3_graph.R

test_that("learn_GENIE3_graph() returns a weighted igraph alongside the importance matrix", {
  skip_if_not_installed("GENIE3")
  out <- quiet(learn_GENIE3_graph(.x, n.trees = 50))   # p x n (features x samples)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))                 # one node per feature
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))          # no self-loops (diagonal zeroed)
  expect_equal(igraph::V(g)$name, rownames(.x))

  imp <- out$weights
  # symmetric p x p importance matrix, feature-named
  expect_true(is.matrix(imp))
  expect_equal(dim(imp), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(imp), rownames(.x))
  expect_equal(colnames(imp), rownames(.x))
  expect_equal(imp, t(imp))
})

test_that("learn_GENIE3_graph() retains fewer edges as the fixed weight threshold rises", {
  skip_if_not_installed("GENIE3")
  loose  <- quiet(learn_GENIE3_graph(.x, n.trees = 50, weight.threshold = 0.001))$graph
  strict <- quiet(learn_GENIE3_graph(.x, n.trees = 50, weight.threshold = 0.05))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_GENIE3_graph() retains fewer edges as the kept percentile rises", {
  skip_if_not_installed("GENIE3")
  # a higher percentile.threshold keeps a smaller top slice of edges
  loose  <- quiet(learn_GENIE3_graph(.x, n.trees = 50, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_GENIE3_graph(.x, n.trees = 50, percentile.threshold = 0.99))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})
