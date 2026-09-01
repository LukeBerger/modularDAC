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
.out     <- quiet(learn_SILGGM_graph(.x))          # learn_* take p x n (features x samples)
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
  expect_named(.out, c("graph", "weights"))
  expect_s3_class(.out$graph, "igraph")

  pcor <- .out$weights
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
  strict <- quiet(learn_SILGGM_graph(.x, fdr.filter = TRUE, max.fdr = 0.001))$graph
  loose  <- quiet(learn_SILGGM_graph(.x, fdr.filter = TRUE, max.fdr = 0.5))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

# ---------------------------------------------------------------------------
# learn_WGCNA_graph()
# ---------------------------------------------------------------------------

test_that("learn_WGCNA_graph() returns a weighted igraph alongside the adjacency matrix", {
  out <- quiet(learn_WGCNA_graph(.x))      # p x n (features x samples)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))                # one node per feature
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))         # no self-loops (diagonal zeroed)
  expect_equal(igraph::V(g)$name, rownames(.x))

  adj <- out$weights
  # symmetric p x p adjacency, feature-named
  expect_true(is.matrix(adj))
  expect_equal(dim(adj), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(adj), rownames(.x))
  expect_equal(colnames(adj), rownames(.x))
  expect_equal(adj, t(adj))
})

test_that("learn_WGCNA_graph() retains fewer edges as the fixed adjacency threshold rises", {
  loose  <- quiet(learn_WGCNA_graph(.x, adj.threshold = 0.01))$graph
  strict <- quiet(learn_WGCNA_graph(.x, adj.threshold = 0.5))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_WGCNA_graph() retains fewer edges as the kept percentile rises", {
  # a higher percentile.threshold keeps a smaller top slice of edges
  loose  <- quiet(learn_WGCNA_graph(.x, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_WGCNA_graph(.x, percentile.threshold = 0.99))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

# ---------------------------------------------------------------------------
# learn_ARACNE_graph()
# ---------------------------------------------------------------------------

test_that("learn_ARACNE_graph() returns a weighted igraph alongside the MI matrix", {
  skip_if_not_installed("minet")
  out <- quiet(learn_ARACNE_graph(.x))     # mutual information over features (p x n input)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))

  mim <- out$weights
  expect_true(is.matrix(mim))
  expect_equal(dim(mim), c(nrow(.x), nrow(.x)))
})

# ---------------------------------------------------------------------------
# learn_CLR_graph()
# ---------------------------------------------------------------------------

test_that("learn_CLR_graph() returns a weighted igraph alongside the CLR score matrix", {
  skip_if_not_installed("minet")
  out <- quiet(learn_CLR_graph(.x))        # mutual information over features (p x n input)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))

  clr <- out$weights
  expect_true(is.matrix(clr))
  expect_equal(dim(clr), c(nrow(.x), nrow(.x)))
  # feature names must survive so the weights matrix can be thresholded
  # externally (the benchmark's edge-budget rule works off dimnames)
  expect_equal(rownames(clr), rownames(.x))
})

test_that("learn_CLR_graph() keeps fewer edges at a stricter percentile", {
  skip_if_not_installed("minet")
  loose  <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.99))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))

  # the upper-triangle quantile means the kept edge count tracks the requested
  # fraction of PAIRS -- unlike a whole-matrix quantile, which double-counts
  n.pairs <- nrow(.x) * (nrow(.x) - 1) / 2
  expect_lt(abs(igraph::gsize(strict) - 0.01 * n.pairs), 0.01 * n.pairs)
})

test_that("learn_CLR_graph() honours an explicit clr.threshold over the percentile", {
  skip_if_not_installed("minet")
  w <- quiet(learn_CLR_graph(.x))$weights
  cut <- as.numeric(stats::quantile(w[upper.tri(w)], 0.80))
  g <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.99, clr.threshold = cut))$graph
  # the explicit cut (80th pct) must win over the percentile argument (99th)
  expect_equal(igraph::gsize(g), sum(w[upper.tri(w)] > cut))
})

# ---------------------------------------------------------------------------
# learn_GENIE3_graph()
# ---------------------------------------------------------------------------

test_that("learn_GENIE3_graph() returns a weighted igraph alongside the importance matrix", {
  skip_if_not_installed("GENIE3")
  out <- quiet(learn_GENIE3_graph(.x, n.trees = 50))   # p x n (features x samples)

  expect_type(out, "list")
  expect_named(out, c("graph", "weights"))

  g <- out$graph
  expect_s3_class(g, "igraph")
  expect_length(g, nrow(.x))                 # one node per feature
  expect_true(igraph::is_weighted(g))
  expect_false(igraph::is_directed(g))
  expect_true(igraph::is_simple(g))          # no self-loops (diagonal zeroed)
  expect_equal(igraph::V(g)$name, rownames(.x))

  imp <- out$weights
  # symmetric p x p importance matrix, feature-named
  expect_true(is.matrix(imp))
  expect_equal(dim(imp), c(nrow(.x), nrow(.x)))
  expect_equal(rownames(imp), rownames(.x))
  expect_equal(colnames(imp), rownames(.x))
  expect_equal(imp, t(imp))
})

test_that("learn_GENIE3_graph() retains fewer edges as the fixed weight threshold rises", {
  skip_if_not_installed("GENIE3")
  loose  <- quiet(learn_GENIE3_graph(.x, n.trees = 50, weight.threshold = 0.001))$graph
  strict <- quiet(learn_GENIE3_graph(.x, n.trees = 50, weight.threshold = 0.05))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_GENIE3_graph() retains fewer edges as the kept percentile rises", {
  skip_if_not_installed("GENIE3")
  # a higher percentile.threshold keeps a smaller top slice of edges
  loose  <- quiet(learn_GENIE3_graph(.x, n.trees = 50, percentile.threshold = 0.90))$graph
  strict <- quiet(learn_GENIE3_graph(.x, n.trees = 50, percentile.threshold = 0.99))$graph
  expect_lte(igraph::gsize(strict), igraph::gsize(loose))
})
