test_that("Fast DaC Run with Defaults", {
  # Graphs and modules
  er <- make_modular_graph()
  x <- sim_graph_data(er, n.samples = 100)
  t <- .true_modules(er)
  f <- eigen_fuzzy_modules(x, t, 90)

  # Run fast dac
  dac <- divide_and_conquer(x, f)
  g <- dac$final.graph

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

test_that("Fast DaC Run with BDgraph", {
  # Graphs and modules
  er <- make_modular_graph()
  x <- sim_graph_data(er, n.samples = 100)
  t <- .true_modules(er)
  f <- eigen_fuzzy_modules(x, t, 90)

  # Run fast dac
  dac <- divide_and_conquer(x, f,
                            graph.learning.func = .bd_wrapper,
                            arg.wrapping.func = .bd_arg_wrapper,
                            out.parsing.func = .default_output_parser,
                            packages.to.each = c("igraph", "BDgraph"),
                            export.to.each = c(".bd_wrapper")
  )
  g <- dac$final.graph

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

