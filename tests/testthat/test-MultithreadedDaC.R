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

  # the default (SILGGM) parser now returns one partial correlation matrix per module
  expect_type(.dac$other.outputs, "list")
  expect_length(.dac$other.outputs, length(.fuzzy@index.list))
  expect_true(all(vapply(.dac$other.outputs, is.matrix, logical(1))))
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
  # each half owns itself (keeps the ownership consistent with the new index list)
  disjoint@core.list  <- disjoint@index.list

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

# ---------------------------------------------------------------------------
# .connect_subgraphs(): weight reconciliation
# ---------------------------------------------------------------------------

# small helper: build a weighted, undirected igraph from a symmetric matrix
.mk_wg <- function(m) {
  igraph::graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE)
}

# two sub-graphs over the same three nodes {a, b, c}. Only the row names of x
# matter to .connect_subgraphs(), so its values are arbitrary.
.xc <- matrix(0, nrow = 3, ncol = 2, dimnames = list(c("a", "b", "c"), NULL))

.A1 <- matrix(0, 3, 3, dimnames = list(c("a","b","c"), c("a","b","c")))
.A1["a","b"] <- .A1["b","a"] <- 0.5
.A1["b","c"] <- .A1["c","b"] <- 0.8
.A1["a","c"] <- .A1["c","a"] <- 0.3   # present only in sub-graph 1

.A2 <- matrix(0, 3, 3, dimnames = list(c("a","b","c"), c("a","b","c")))
.A2["a","b"] <- .A2["b","a"] <- 0.4
.A2["b","c"] <- .A2["c","b"] <- 0.2
# no a-c edge in sub-graph 2

.g1 <- .mk_wg(.A1)
.g2 <- .mk_wg(.A2)

wmat <- function(g) {
  m <- as.matrix(igraph::as_adjacency_matrix(g, attr = "weight"))
  m[c("a","b","c"), c("a","b","c")]
}

test_that(".connect_subgraphs() min keeps the smallest-magnitude weight and drops non-consensus edges", {
  g <- .connect_subgraphs(.xc, list(.g1, .g2), "min")
  expect_s3_class(g, "igraph")
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))

  W <- wmat(g)
  expect_equal(W["a","b"], 0.4)   # min(|0.5|, |0.4|)
  expect_equal(W["b","c"], 0.2)   # min(|0.8|, |0.2|)
  expect_equal(W["a","c"], 0)     # possible in both, present in one -> dropped
})

test_that(".connect_subgraphs() mean averages signed weights and keeps soft edges", {
  g <- .connect_subgraphs(.xc, list(.g1, .g2), "mean")
  W <- wmat(g)
  expect_equal(W["a","b"], 0.45)            # (0.5 + 0.4) / 2
  expect_equal(W["b","c"], 0.5)             # (0.8 + 0.2) / 2
  expect_equal(W["a","c"], 0.3 / 2)         # (0.3 + 0) / 2, survives softly
})

test_that(".connect_subgraphs() defaults to min", {
  expect_equal(wmat(.connect_subgraphs(.xc, list(.g1, .g2))),
               wmat(.connect_subgraphs(.xc, list(.g1, .g2), "min")))
})

test_that(".connect_subgraphs() warns on sign disagreement and keeps the signed min-magnitude value", {
  A2neg <- .A2
  A2neg["a","b"] <- A2neg["b","a"] <- -0.4   # disagrees in sign with g1 (+0.5)
  g2neg <- .mk_wg(A2neg)

  expect_warning(
    g <- .connect_subgraphs(.xc, list(.g1, g2neg), "min"),
    "different signs"
  )
  W <- suppressWarnings(wmat(.connect_subgraphs(.xc, list(.g1, g2neg), "min")))
  expect_equal(W["a","b"], -0.4)   # smaller magnitude wins, sign retained
})

test_that(".connect_subgraphs() rejects an unknown weight.summary", {
  expect_error(.connect_subgraphs(.xc, list(.g1, .g2), "median"))
})

# ---------------------------------------------------------------------------
# divide_and_conquer(): weight.summary plumbing
# ---------------------------------------------------------------------------

test_that("divide_and_conquer() returns a weighted final graph and respects weight.summary", {
  expect_true(igraph::is_weighted(.dac$final.graph))   # .dac uses the default (min)

  dac.mean <- quiet(divide_and_conquer(.x, .fuzzy, weight.summary = "mean"))
  # mean keeps every edge present in any possible sub-graph, so it is a superset
  # of the consensus (min) edge set
  expect_gte(igraph::gsize(dac.mean$final.graph), igraph::gsize(.dac$final.graph))
})
