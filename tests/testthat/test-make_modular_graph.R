# Tests for make_modular_graph() in R/make_modular_graph.R
#
# The generator is stochastic, so every test seeds the RNG.

test_that("make_modular_graph() returns a simple igraph with the requested size", {
  set.seed(1)
  g <- make_modular_graph()

  expect_s3_class(g, "igraph")
  # default n.nodes = 120 split evenly across n.mods = 3 -> 120 nodes total
  expect_length(g, 120)
  # simplify() should leave no self-loops or multiple edges
  expect_true(igraph::is_simple(g))
})

test_that("make_modular_graph() labels every node with a module and a name", {
  set.seed(1)
  g <- make_modular_graph()

  # module attribute: one per node, exactly n.mods distinct values
  expect_length(igraph::V(g)$module, 120)
  expect_equal(length(unique(igraph::V(g)$module)), 3)

  # name attribute: one unique "Node_*" label per node
  expect_length(igraph::V(g)$name, 120)
  expect_true(all(grepl("^Node_", igraph::V(g)$name)))
  expect_equal(anyDuplicated(igraph::V(g)$name), 0)
})

test_that("make_modular_graph() honours n.mods and n.nodes", {
  set.seed(1)
  g <- make_modular_graph(n.mods = 4, n.nodes = 80)

  # nodes-per-module = round(80 / 4) = 20, so total is 80
  expect_length(g, 80)
  expect_equal(length(unique(igraph::V(g)$module)), 4)
  expect_equal(as.vector(table(igraph::V(g)$module)), rep(20, 4))
})

test_that("make_modular_graph() supports the scale-free generator", {
  set.seed(1)
  g <- make_modular_graph(g.type = "sf")

  expect_s3_class(g, "igraph")
  expect_length(g, 120)
  expect_equal(length(unique(igraph::V(g)$module)), 3)
})

test_that("make_modular_graph() leaves no isolated nodes when no.uncon = TRUE", {
  set.seed(1)
  g <- make_modular_graph(no.uncon = TRUE)

  expect_gte(min(igraph::degree(g)), 1)
})

test_that("make_modular_graph() adds inter-module edges when link.all = TRUE", {
  set.seed(1)
  g <- make_modular_graph(link.all = TRUE)

  # an edge is "inter-module" when its endpoints carry different module labels
  el <- igraph::as_edgelist(g, names = FALSE)
  mods <- igraph::V(g)$module
  cross.module.edges <- sum(mods[el[, 1]] != mods[el[, 2]])

  expect_gt(cross.module.edges, 0)
})

test_that("make_modular_graph() rejects an invalid graph type", {
  expect_error(make_modular_graph(g.type = "not-a-type"))
})
