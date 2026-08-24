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

test_that("divide_and_conquer() returns subgraphs, a graph, weights, and other outputs", {
  expect_type(.dac, "list")
  expect_named(.dac, c("module.subgraphs", "graph", "weights", "other.outputs"))

  # one learned igraph per input module
  expect_length(.dac$module.subgraphs, length(.fuzzy@index.list))
  expect_true(all(vapply(.dac$module.subgraphs, igraph::is_igraph, logical(1))))

  # the default (SILGGM) parser now returns one partial correlation matrix per module
  expect_type(.dac$other.outputs, "list")
  expect_length(.dac$other.outputs, length(.fuzzy@index.list))
  expect_true(all(vapply(.dac$other.outputs, is.matrix, logical(1))))
})

test_that("divide_and_conquer() stitches modules back into the full feature set", {
  fg <- .dac$graph
  expect_s3_class(fg, "igraph")
  # every feature from the input data appears exactly once in the final graph
  expect_length(fg, nrow(.x))
  expect_setequal(igraph::V(fg)$name, rownames(.x))
})

test_that("divide_and_conquer() returns a combined feature-by-feature weight matrix", {
  W <- .dac$weights
  # a square, symmetric, feature-named matrix over the full feature set (like the
  # single-shot learners' 'weights' output)
  expect_true(is.matrix(W))
  expect_equal(dim(W), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(W), rownames(.x))
  expect_equal(colnames(W), rownames(.x))
  expect_equal(W, t(W))

  # output.weights = FALSE skips the matrix entirely
  dac.now <- quiet(divide_and_conquer(.x, .fuzzy, output.weights = FALSE))
  expect_null(dac.now$weights)
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
# .connect_weights(): combining per-module weight matrices
# ---------------------------------------------------------------------------

test_that(".connect_weights() min keeps the smallest-magnitude weight per pair", {
  # both modules cover {a, b, c}; ownership-agnostic (core.sets = NULL)
  W <- .connect_weights(.xc, list(.A1, .A2), "min")
  expect_equal(W, t(W))                       # symmetric
  expect_equal(dim(W), c(3, 3))
  expect_equal(W["a","b"], 0.4)               # min(|0.5|, |0.4|)
  expect_equal(W["b","c"], 0.2)               # min(|0.8|, |0.2|)
  expect_equal(W["a","c"], 0)                 # 0.3 in one module, 0 in the other -> 0
  expect_equal(diag(W), c(a = 0, b = 0, c = 0))
})

test_that(".connect_weights() mean averages the signed weights over proposers", {
  W <- .connect_weights(.xc, list(.A1, .A2), "mean")
  expect_equal(W["a","b"], 0.45)              # (0.5 + 0.4) / 2
  expect_equal(W["b","c"], 0.5)               # (0.8 + 0.2) / 2
  expect_equal(W["a","c"], 0.15)              # (0.3 + 0) / 2
})

test_that(".connect_weights() honours core ownership (drops auxiliary-auxiliary pairs)", {
  # one module over {a, b, c, d}; only a and b are core, so c and d are auxiliary
  xc4 <- matrix(0, 4, 2, dimnames = list(c("a","b","c","d"), NULL))
  A <- matrix(0.7, 4, 4, dimnames = list(c("a","b","c","d"), c("a","b","c","d")))
  diag(A) <- 0

  W <- .connect_weights(xc4, list(A), "min", core.sets = list(c("a","b")))
  # pairs touching a core node survive with their weight ...
  expect_equal(W["a","c"], 0.7)
  expect_equal(W["b","d"], 0.7)
  # ... but the auxiliary-auxiliary pair c-d is discarded
  expect_equal(W["c","d"], 0)
})

# ---------------------------------------------------------------------------
# .subgraph_weight_mats(): choosing a weight source per module
# ---------------------------------------------------------------------------

test_that(".subgraph_weight_mats() prefers a parser weight matrix, then the weighted graph, else NULL", {
  # module 1: parser supplies a feature-named weight matrix -> use it directly
  # module 2: parser output is not a matrix, but the graph is weighted -> use its adjacency
  # module 3: no weight matrix and an unweighted graph -> NULL
  g.unweighted <- igraph::graph_from_data_frame(
    data.frame(from = "a", to = "b"), directed = FALSE,
    vertices = data.frame(name = c("a","b","c")))

  wm <- .subgraph_weight_mats(
    learned.graphs = list(.g1, .g2, g.unweighted),
    other.outputs  = list(.A1, "not-a-matrix", NULL),
    x              = .xc)

  expect_identical(wm[[1]], .A1)                     # source 1: parser matrix
  expect_true(is.matrix(wm[[2]]))                    # source 2: from weighted graph
  expect_equal(wm[[2]][c("a","b","c"), c("a","b","c")], .A2)
  expect_null(wm[[3]])                               # source 3: nothing usable
})

# ---------------------------------------------------------------------------
# divide_and_conquer(): output.weights fallback
# ---------------------------------------------------------------------------

# a stub learner that returns an UNWEIGHTED graph and no weight matrix, used to
# exercise the output.weights fallback (warn + NULL, but do not stop)
.stub_arg <- function(sub.x, ...) lapply(sub.x, function(xx) list(x = xx))
.stub_learn <- function(x) {
  ns <- rownames(x)
  edges <- if (length(ns) >= 2) {
    data.frame(from = ns[-length(ns)], to = ns[-1], stringsAsFactors = FALSE)
  } else {
    data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  }
  g <- igraph::graph_from_data_frame(edges, directed = FALSE,
                                     vertices = data.frame(name = ns))
  list(graph = g, weights = NULL)   # unweighted graph, no weight matrix
}
.stub_parse <- function(outs) list(
  learned.graphs = lapply(outs, function(o) o$graph),
  other.outputs  = lapply(outs, function(o) o$weights)
)

test_that("divide_and_conquer() warns and returns NULL weights when learners expose no weights", {
  dac <- NULL
  warns <- testthat::capture_warnings(
    dac <- suppressMessages(divide_and_conquer(
      .x, .fuzzy,
      graph.learning.func = .stub_learn,
      arg.wrapping.func   = .stub_arg,
      out.parsing.func    = .stub_parse,
      packages.to.each    = "igraph",
      export.to.each      = character(0)))
  )
  # the missing-weights warning is emitted (alongside foreach's %dopar% notice)
  expect_match(warns, "neither weight matrices nor weighted graphs", all = FALSE)
  # weights are NULL, but the run still completes and stitches a final graph
  expect_null(dac$weights)
  expect_s3_class(dac$graph, "igraph")
  expect_setequal(igraph::V(dac$graph)$name, rownames(.x))
})

# ---------------------------------------------------------------------------
# divide_and_conquer(): weight.summary plumbing
# ---------------------------------------------------------------------------

test_that("divide_and_conquer() returns a weighted final graph and respects weight.summary", {
  expect_true(igraph::is_weighted(.dac$graph))   # .dac uses the default (min)

  dac.mean <- quiet(divide_and_conquer(.x, .fuzzy, weight.summary = "mean"))
  # mean keeps every edge present in any possible sub-graph, so it is a superset
  # of the consensus (min) edge set
  expect_gte(igraph::gsize(dac.mean$graph), igraph::gsize(.dac$graph))
})
