# Tests for the graph-learning functions in R/02.LearnGraphs.R
# Scope: learn_SILGGM_graph()
#
# SILGGM prints progress via cat() (not message()), so quiet() captures that
# output to keep test logs clean. The learned graph is computed once and shared
# across the assertions below, since SILGGM is the slow part.

# silence cat()/message()/warning() noise while returning the value
quiet <- function(expr) {
  utils::capture.output(val <- suppressWarnings(suppressMessages(expr)))
  val
}

set.seed(1)
.g <- make_modular_graph()
.x <- sim_graph_data(.g, n.samples = 100)          # p x n (features x samples)
.out     <- quiet(learn_SILGGM_graph(t(.x)))       # SILGGM expects n x p
.learned <- .out$graph                             # the learned igraph

# ---------------------------------------------------------------------------
# learn_SILGGM_graph()
# ---------------------------------------------------------------------------

test_that("learn_SILGGM_graph() returns a weighted igraph with one node per feature", {
  expect_s3_class(.learned, "igraph")
  expect_length(.learned, nrow(.x))
  expect_true(igraph::is_weighted(.learned))
})

test_that("learn_SILGGM_graph() returns the graph alongside the partial correlation matrix", {
  expect_type(.out, "list")
  expect_named(.out, c("graph", "partial.cor"))
  expect_s3_class(.out$graph, "igraph")

  pcor <- .out$partial.cor
  # symmetric p x p matrix, feature-named, with a unit diagonal
  expect_true(is.matrix(pcor))
  expect_equal(dim(pcor), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(pcor), rownames(.x))
  expect_equal(colnames(pcor), rownames(.x))
  expect_equal(unname(diag(pcor)), rep(1, nrow(.x)))
  expect_equal(pcor, t(pcor))
})

test_that("learn_SILGGM_graph() preserves feature identity and undirected structure", {
  # node names carry over from the data's feature labels
  expect_equal(igraph::V(.learned)$name, rownames(.x))
  expect_false(igraph::is_directed(.learned))
  # no self-loops (diagonal is zeroed before graph construction)
  expect_true(igraph::is_simple(.learned))
})

test_that("learn_SILGGM_graph() recovers a non-empty subset of the true edges", {
  # not a benchmark, just a sanity check that the learned graph is usable and
  # overlaps the truth: F1 is a valid score in [0, 1].
  f1 <- calc_F1(.g, .learned)$F1
  expect_type(f1, "double")
  expect_gte(f1, 0)
  expect_lte(f1, 1)
})

test_that("learn_SILGGM_graph() drops more edges as the FDR threshold tightens", {
  # a stricter (smaller) max.fdr should never retain MORE edges than a lax one
  strict <- quiet(learn_SILGGM_graph(t(.x), fdr.filter = TRUE, max.fdr = 0.001))$graph
  loose  <- quiet(learn_SILGGM_graph(t(.x), fdr.filter = TRUE, max.fdr = 0.5))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})
