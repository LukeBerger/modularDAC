# Tests for the divide-and-conquer driver in R/03.Multithreaded_DaC.R
# Scope: divide_and_conquer()
#
# divide_and_conquer() splits the data into overlapping modules, learns a graph
# per module (SILGGM by default), and stitches them back together. The default
# learner is verbose via cat(), so quiet() captures that output and the full run
# is computed once and shared across assertions.

# silence cat()/message()/warning() noise while returning the value
quiet <- function(expr) {
  utils::capture.output(val <- suppressWarnings(suppressMessages(expr)))
  val
}

set.seed(1)
.g <- make_modular_graph()
.x <- sim_graph_data(.g, n.samples = 100)   # p x n
.fuzzy <- true_fuzzy(true_modules(.g), .g)  # overlapping modules to divide along
.dac <- quiet(divide_and_conquer(.x, .fuzzy))

# ---------------------------------------------------------------------------
# divide_and_conquer(): successful run with defaults
# ---------------------------------------------------------------------------

test_that("divide_and_conquer() returns subgraphs, a final graph, and other outputs", {
  expect_type(.dac, "list")
  expect_named(.dac, c("modular.subgraphs", "final.graph", "other.outputs"))

  # one learned igraph per input module
  expect_length(.dac$modular.subgraphs, length(.fuzzy@index.list))
  expect_true(all(vapply(.dac$modular.subgraphs, igraph::is_igraph, logical(1))))

  # default output parser carries no extra outputs
  expect_null(.dac$other.outputs)
})

test_that("divide_and_conquer() stitches modules back into the full feature set", {
  fg <- .dac$final.graph
  expect_s3_class(fg, "igraph")
  # every feature from the input data appears exactly once in the final graph
  expect_length(fg, nrow(.x))
  expect_setequal(igraph::V(fg)$name, rownames(.x))
})

# ---------------------------------------------------------------------------
# divide_and_conquer(): input validation
# ---------------------------------------------------------------------------

test_that("divide_and_conquer() rejects non-matrix data", {
  expect_error(divide_and_conquer(as.data.frame(.x), .fuzzy),
               regexp = "numeric matrix")
})

test_that("divide_and_conquer() requires every row to be covered by a module", {
  # drop the high-index features from every module so some rows are uncovered
  uncovered <- .fuzzy
  uncovered@index.list <- lapply(.fuzzy@index.list, function(v) v[v <= 100])
  expect_error(suppressMessages(divide_and_conquer(.x, uncovered)),
               regexp = "covered")
})

test_that("divide_and_conquer() rejects a non-function learner", {
  expect_error(divide_and_conquer(.x, .fuzzy, graph.learning.func = 42),
               regexp = "function")
})

test_that("divide_and_conquer() warns when a module has no overlap with the others", {
  # build two fully disjoint modules covering all rows but sharing no nodes
  p <- nrow(.x)
  disjoint <- .fuzzy
  disjoint@index.list <- list(seq_len(p / 2), (p / 2 + 1):p)
  disjoint@name.list  <- list(rownames(.x)[seq_len(p / 2)],
                              rownames(.x)[(p / 2 + 1):p])

  # capture_warnings() collects every warning (the run also emits an unrelated
  # "%dopar% sequentially" warning), so we assert that at least one reports the
  # missing overlap rather than matching a single warning.
  warns <- testthat::capture_warnings(
    utils::capture.output(suppressMessages(divide_and_conquer(.x, disjoint)))
  )
  expect_match(warns, "no overlaps", all = FALSE)
})

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
  expect_named(dac, c("modular.subgraphs", "final.graph", "other.outputs"))
  expect_length(dac$modular.subgraphs, length(.fuzzy@index.list))
  expect_true(all(vapply(dac$modular.subgraphs, igraph::is_igraph, logical(1))))

  # stitched graph covers the full feature set
  expect_s3_class(dac$final.graph, "igraph")
  expect_length(dac$final.graph, nrow(.x))
  expect_setequal(igraph::V(dac$final.graph)$name, rownames(.x))

  # unlike the default parser, the RSNet parser returns the extra consensus
  # outputs alongside the learned graphs
  expect_false(is.null(dac$other.outputs))
  expect_length(dac$other.outputs, length(.fuzzy@index.list))
})
