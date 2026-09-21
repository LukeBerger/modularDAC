# Tests for sim_graph_data() in R/sim_graph_data.R
#
# Sampling is stochastic, so every test seeds the RNG.

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
