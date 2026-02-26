test_that("SILGMM Learning", {
  # Simulate True Graph and Data
  er <- make_modular_graph()
  x <- sim_graph_data(er, n.samples = 100)

  # Learn Graph
  g <- learn_SILGGM_graph(t(x))

  # Check Graph
  expect_s3_class(g, "igraph")
  expect_length(g, 120)
  expect_true(all(igraph::V(g)$name == row.names(x)))

  # Check F1
  f1 <-calc_F1(er, g)$F1
  expect_type(f1, "double")
  expect_length(f1, 1)
  expect_gte(f1, 0)
  expect_lte(f1, 1)
  expect_equal(calc_F1(er, er)$F1, 1)

})

# Commenting out when not working on this function since it slows the run
# test_that("WGCNA Learning", {
#   # Simulate True Graph and Data
#   er <- make_modular_graph()
#   x <- sim_graph_data(er, n.samples = 100)
#
#   # Learn Graph
#   g <- learn_WGCNA_graph(t(x))
#
#   # Check Graph
#   expect_s3_class(g, "igraph")
#   expect_length(g, 120)
#   expect_true(all(igraph::V(g)$name == row.names(x)))
#
#   # Check F1
#   f1 <-calc_F1(er, g)$F1
#   expect_type(f1, "double")
#   expect_length(f1, 1)
#   expect_gte(f1, 0)
#   expect_lte(f1, 1)
#
# })


test_that("ARACNE Learning", {
  # Simulate True Graph and Data
  er <- make_modular_graph()
  x <- sim_graph_data(er, n.samples = 100)

  # Learn Graph
  g <- learn_ARACNE_graph(t(x))

  # Check Graph
  expect_s3_class(g, "igraph")
  expect_length(g, 120)
  expect_true(all(igraph::V(g)$name == row.names(x)))

  # Check F1
  f1 <-calc_F1(er, g)$F1
  expect_type(f1, "double")
  expect_length(f1, 1)
  expect_gte(f1, 0)
  expect_lte(f1, 1)

})
