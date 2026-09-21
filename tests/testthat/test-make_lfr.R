# Tests for make_lfr() in R/make_lfr.R
#
# The generator is stochastic, so every test seeds the RNG.

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
