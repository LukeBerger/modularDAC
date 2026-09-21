# Tests for learn_WGCNA_graph() in R/learn_WGCNA_graph.R

test_that("learn_WGCNA_graph() returns a weighted igraph alongside the adjacency matrix", {
  out <- quiet(learn_WGCNA_graph(.x))      # p x n (features x samples)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))                # one node per feature
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))         # no self-loops (diagonal zeroed)
  expect_equal(igraph::V(g)$name, rownames(.x))

  adj <- out$weights
  # symmetric p x p adjacency, feature-named
  expect_true(is.matrix(adj))
  expect_equal(dim(adj), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(adj), rownames(.x))
  expect_equal(colnames(adj), rownames(.x))
  expect_equal(adj, t(adj))
})

test_that("learn_WGCNA_graph() retains fewer edges as the fixed adjacency threshold rises", {
  loose  <- quiet(learn_WGCNA_graph(.x, adj.threshold = 0.01))$graph
  strict <- quiet(learn_WGCNA_graph(.x, adj.threshold = 0.5))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_WGCNA_graph() retains fewer edges as the kept percentile rises", {
  # a higher percentile.threshold keeps a smaller top slice of edges
  loose  <- quiet(learn_WGCNA_graph(.x, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_WGCNA_graph(.x, percentile.threshold = 0.99))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})
