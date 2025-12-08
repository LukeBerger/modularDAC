test_that("Graph Simulation", {
  # Create graphs
  er <- make_modular_graph()
  lfr <- make_lfr()
  hybrid <- make_submodular_graph()


  # Test class
  expect_s3_class(er, "igraph")
  expect_s3_class(lfr, "igraph")
  expect_s3_class(hybrid, "igraph")


  # Test default length
  expect_length(er, 120)
  expect_length(lfr, 120)
  expect_length(hybrid, 360)


  # Test Modularity
  expect_length(igraph::V(er)$module, 120)
  expect_length(unique(igraph::V(er)$module), 3)

  expect_length(igraph::V(lfr)$module, 120)
  expect_length(unique(igraph::V(lfr)$module), 3)

  expect_length(igraph::V(hybrid)$subgraph, 360)
  expect_length(unique(igraph::V(hybrid)$subgraph), 3)
})

test_that("Plotting", {
  vn <- modular_plot(make_modular_graph())
  expect_s3_class(vn, "visNetwork")
  expect_s3_class(vn, "htmlwidget")
})

test_that("Data Simulation", {
  # Simualte Data
  x <- sim_graph_data(make_modular_graph(), n.samples = 10)

  # Check that it is a numeric matrix, with no NAs
  expect_true(is.matrix(x) && is.numeric(x))
  expect_equal(sum(is.na(x)), 0)

  # Check dimensions of 120 x 10
  expect_equal(ncol(x), 10)
  expect_equal(nrow(x), 120)
})

