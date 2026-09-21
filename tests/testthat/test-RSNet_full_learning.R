# Tests for the RSNet wrapper set in R/RSNet_full_learning.R
#
# .RSNet_arg_wrapper / .RSNet_full_learning / .RSNet_output_parser are
# exercised through divide_and_conquer(), which is how they are used.

# ---------------------------------------------------------------------------
# divide_and_conquer(): alternative (RSNet) graph-learning helpers
# ---------------------------------------------------------------------------
test_that("divide_and_conquer() works with the RSNet helper functions", {
  skip_if_not_installed("RSNet")

  # RSNet parallelises with parallel::mclapply; pass n_cores = 1 so it runs
  # serially (mc.cores > 1 is unsupported on Windows). n_cores flows through
  # the ... to .RSNet_arg_wrapper and on to RSNet::ensemble_ggm.
  dac <- quiet(divide_and_conquer(
    .x, .fuzzy,
    graph.learning.func = .RSNet_full_learning,
    arg.wrapping.func   = .RSNet_arg_wrapper,
    out.parsing.func    = .RSNet_output_parser,
    packages.to.each    = c("igraph", "RSNet"),
    export.to.each      = c(".RSNet_full_learning"),
    n_cores = 1
  ))

  # same container contract as the default SILGGM path
  expect_named(dac, c("module.subgraphs", "graph", "weights", "other.outputs"))
  expect_length(dac$module.subgraphs, length(.fuzzy@index.list))
  expect_true(all(vapply(dac$module.subgraphs, igraph::is_igraph, logical(1))))

  # stitched graph covers the full feature set
  expect_s3_class(dac$graph, "igraph")
  expect_length(dac$graph, nrow(.x))
  expect_setequal(igraph::V(dac$graph)$name, rownames(.x))

  # unlike the default parser, the RSNet parser returns the extra consensus
  # outputs alongside the learned graphs
  expect_false(is.null(dac$other.outputs))
  expect_length(dac$other.outputs, length(.fuzzy@index.list))
})
