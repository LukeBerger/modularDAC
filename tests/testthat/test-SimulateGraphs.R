# Tests for the graph- and data-simulation functions in R/00.SimulateGraphs.R
# Scope: make_modular_graph(), sim_graph_data()
#
# These functions are stochastic, so every test seeds the RNG to keep results
# deterministic and the assertions reproducible.

# ---------------------------------------------------------------------------
# make_modular_graph()
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# sim_graph_data()
# ---------------------------------------------------------------------------

test_that("sim_graph_data() returns a clean numeric p x n matrix", {
  set.seed(1)
  g <- make_modular_graph()
  x <- sim_graph_data(g, n.samples = 10)

  expect_true(is.matrix(x))
  expect_true(is.numeric(x))
  expect_false(anyNA(x))

  # p features (graph nodes) as rows, n samples as columns
  expect_equal(nrow(x), length(g))
  expect_equal(ncol(x), 10)
})

test_that("sim_graph_data() names rows by node and columns by sample", {
  set.seed(1)
  g <- make_modular_graph()
  x <- sim_graph_data(g, n.samples = 10)

  expect_equal(rownames(x), igraph::V(g)$name)
  expect_equal(colnames(x), paste("Sample", 1:10, sep = "_"))
})

test_that("sim_graph_data() defaults to zero-mean features", {
  set.seed(1)
  g <- make_modular_graph()
  # large sample so the empirical mean is a stable estimate of the true mean
  x <- sim_graph_data(g, n.samples = 2000)

  expect_equal(rowMeans(x), rep(0, nrow(x)), tolerance = 0.2,
               ignore_attr = TRUE)
})

test_that("sim_graph_data() shifts features toward a supplied mean vector", {
  set.seed(1)
  g <- make_modular_graph()
  mean.vec <- rep(5, length(g))
  x <- sim_graph_data(g, n.samples = 2000, mean.vec = mean.vec)

  expect_equal(nrow(x), length(g))
  expect_equal(rowMeans(x), mean.vec, tolerance = 0.2, ignore_attr = TRUE)
})

# ---------------------------------------------------------------------------
# make_lfr()
# ---------------------------------------------------------------------------

test_that("make_lfr() returns an igraph of the requested size", {
  set.seed(1)
  g <- make_lfr()

  expect_s3_class(g, "igraph")
  # default n = 120
  expect_length(g, 120)

  set.seed(1)
  expect_length(make_lfr(n = 200), 200)
})

test_that("make_lfr() labels nodes with names and community membership", {
  set.seed(1)
  g <- make_lfr()

  # sequential "Node_*" names, one unique label per node
  expect_length(igraph::V(g)$name, 120)
  expect_true(all(grepl("^Node_", igraph::V(g)$name)))
  expect_equal(anyDuplicated(igraph::V(g)$name), 0)

  # module is copied from the LFR community membership
  expect_length(igraph::V(g)$module, 120)
  expect_gt(length(unique(igraph::V(g)$module)), 1)
  expect_equal(igraph::V(g)$module, igraph::V(g)$membership)
})
