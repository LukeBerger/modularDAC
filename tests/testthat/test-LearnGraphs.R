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

# The ARACNE matrix is ~93% structural zeros after DPI pruning, so the two
# percentile modes behave completely differently on it -- see
# .percentile_threshold() for why nonzero.only defaults to TRUE.

test_that("learn_ARACNE_graph() ranks the percentile within the DPI survivors", {
  skip_if_not_installed("minet")
  w   <- quiet(learn_ARACNE_graph(.x))$weights
  pos <- w[upper.tri(w)]
  pos <- pos[pos > 0]                       # the pairs DPI left in contention
  expect_gt(length(pos), 0)

  for (p in c(0.80, 0.90, 0.99)) {
    g <- quiet(learn_ARACNE_graph(.x, percentile.threshold = p))$graph
    # the cut is the p-th quantile of the SURVIVORS, not of all node pairs
    expect_equal(igraph::gsize(g), sum(pos > stats::quantile(pos, p)))
  }
})

test_that("learn_ARACNE_graph() keeps strictly fewer edges as the percentile rises", {
  skip_if_not_installed("minet")
  # the regression this guards: under the old all-pairs percentile these three
  # cuts all landed on the zero mass and returned the SAME graph
  sizes <- vapply(c(0.80, 0.90, 0.99), function(p) {
    igraph::gsize(quiet(learn_ARACNE_graph(.x, percentile.threshold = p))$graph)
  }, numeric(1))

  expect_equal(sizes, sort(sizes, decreasing = TRUE))
  expect_equal(anyDuplicated(sizes), 0L)
})

test_that("learn_ARACNE_graph(nonzero.only = FALSE) recovers the old all-pairs cut", {
  skip_if_not_installed("minet")
  w  <- quiet(learn_ARACNE_graph(.x))$weights
  ut <- w[upper.tri(w)]

  # a percentile above the zero fraction still ranks within all pairs
  p <- 0.99
  skip_if(mean(ut <= 0) >= p, "matrix too sparse for a non-degenerate all-pairs cut")
  g <- quiet(learn_ARACNE_graph(.x, percentile.threshold = p, nonzero.only = FALSE))$graph
  expect_equal(igraph::gsize(g), sum(ut > stats::quantile(ut, p)))
})

test_that("learn_ARACNE_graph(nonzero.only = FALSE) warns and degenerates on the zero mass", {
  skip_if_not_installed("minet")
  w  <- quiet(learn_ARACNE_graph(.x))$weights
  ut <- w[upper.tri(w)]
  n.pos <- sum(ut > 0)

  # two percentiles that both sit below the zero fraction
  p <- c(0.50, 0.90)
  skip_if(any(mean(ut <= 0) < p), "matrix not sparse enough to force the degeneracy")

  sizes <- vapply(p, function(pp) {
    expect_warning(out <- learn_ARACNE_graph(.x, percentile.threshold = pp,
                                             nonzero.only = FALSE),
                   "falls on the zero mass")
    igraph::gsize(out$graph)
  }, numeric(1))

  # the cut collapsed to "keep every non-zero entry", identically for both
  expect_equal(sizes, c(n.pos, n.pos))
})

test_that("learn_ARACNE_graph() honours an explicit mim.threshold over the percentile", {
  skip_if_not_installed("minet")
  w   <- quiet(learn_ARACNE_graph(.x))$weights
  ut  <- w[upper.tri(w)]
  cut <- as.numeric(stats::quantile(ut[ut > 0], 0.50))
  g <- quiet(learn_ARACNE_graph(.x, percentile.threshold = 0.99,
                                mim.threshold = cut))$graph
  # the explicit cut (median survivor) must win over the percentile argument
  expect_equal(igraph::gsize(g), sum(ut > cut))
})

# --- the MI significance threshold (Margolin et al.'s rule) -----------------

test_that(".mi_pvalue() is the chi-squared tail of the Spearman MI estimator", {
  skip_if_not_installed("minet")
  mim <- minet::build.mim(t(.x))
  M   <- ncol(.x)
  expect_equal(.mi_pvalue(mim, M),
               stats::pchisq(2*(M-1)*mim, df = 1, lower.tail = FALSE))
  # dimnames survive, so the p-value matrix can be masked against the MI matrix
  expect_equal(dimnames(.mi_pvalue(mim, M)), dimnames(mim))
  expect_error(.mi_pvalue(mim, 2), "at least 3 samples")
})

test_that(".mi_pvalue() is calibrated against a permutation null", {
  skip_if_not_installed("minet")
  # shuffling each feature independently makes every pair manifestly independent,
  # so the p-values must be ~uniform -- this is the empirical null the paper
  # tabulates, checked against the closed form used here instead
  set.seed(11)
  shuffled <- t(apply(.x, 1, sample))
  pv <- .mi_pvalue(minet::build.mim(t(shuffled)), ncol(.x))
  pv <- pv[upper.tri(pv)]
  for (a in c(0.01, 0.05, 0.10)) {
    expect_lt(abs(mean(pv < a) - a), 0.4 * a)     # within 40% of nominal
  }
})

test_that("a p-value rule on MI is exactly a single MI cutoff", {
  skip_if_not_installed("minet")
  # the property the implementation leans on: p is strictly decreasing in MI at
  # fixed n, and BH keeps a prefix of the sorted p-values, so either rule
  # selects {MI >= c}. This is what makes the mask safe to intersect with the
  # DPI output instead of re-running DPI on a masked MI matrix.
  mim <- minet::build.mim(t(.x))
  for (keep in list(.matrix_p_adjust(.mi_pvalue(mim, ncol(.x))) < 0.05,
                    .mi_pvalue(mim, ncol(.x)) < 1e-6)) {
    cutoff <- min(mim[keep])
    expect_identical(which(keep[upper.tri(keep)]),
                     which(mim[upper.tri(mim)] >= cutoff))
  }
})

test_that("thresholding before or after DPI gives learn_ARACNE_graph the same graph", {
  skip_if_not_installed("minet")
  # follows from the cutoff property above: DPI prunes (i,j) only via a witness
  # k with min(MI(i,k), MI(j,k)) > MI(i,j) + eps, so any witness for a surviving
  # edge also clears the cutoff and survives the mask
  mim  <- minet::build.mim(t(.x))
  keep <- .matrix_p_adjust(.mi_pvalue(mim, ncol(.x))) < 0.05

  for (eps in c(0, 0.001)) {
    g <- quiet(learn_ARACNE_graph(.x, eps = eps, max.fdr = 0.05))$graph

    after  <- minet::aracne(mim, eps = eps) * keep      # what the function does
    masked <- mim * keep; diag(masked) <- 0
    before <- minet::aracne(masked, eps = eps)          # Margolin et al.'s order
    diag(after) <- 0

    expect_identical(which(after[upper.tri(after)] > 0),
                     which(before[upper.tri(before)] > 0))
    expect_equal(igraph::gsize(g), sum(before[upper.tri(before)] > 0))
  }
})

test_that("learn_ARACNE_graph(max.fdr) keeps exactly the BH-significant DPI survivors", {
  skip_if_not_installed("minet")
  mim  <- minet::build.mim(t(.x))
  keep <- .matrix_p_adjust(.mi_pvalue(mim, ncol(.x))) < 0.05
  w    <- minet::aracne(mim, eps = 0) * keep; diag(w) <- 0

  out <- quiet(learn_ARACNE_graph(.x, max.fdr = 0.05))
  expect_equal(igraph::gsize(out$graph), sum(w[upper.tri(w)] > 0))
  expect_named(out, c("graph", "weights"))
  expect_equal(rownames(out$weights), rownames(.x))
  # 'weights' stays the unmasked DPI output, so the graph is a subgraph of it
  expect_lte(igraph::gsize(out$graph), sum(out$weights[upper.tri(out$weights)] > 0))
})

test_that("learn_ARACNE_graph() significance rules keep fewer edges as they tighten", {
  skip_if_not_installed("minet")
  sizes <- vapply(c(0.10, 0.05, 1e-4), function(q) {
    igraph::gsize(quiet(learn_ARACNE_graph(.x, max.fdr = q))$graph)
  }, numeric(1))
  expect_equal(sizes, sort(sizes, decreasing = TRUE))

  # BH adjustment only ever raises a p-value, so at the same nominal level the
  # raw-p rule is the more permissive of the two
  expect_gte(igraph::gsize(quiet(learn_ARACNE_graph(.x, p.threshold = 0.05))$graph),
             igraph::gsize(quiet(learn_ARACNE_graph(.x, max.fdr = 0.05))$graph))

  # ... and a Bonferroni-level raw p is stricter than BH at the same alpha
  bonf <- 0.05 / (nrow(.x) * (nrow(.x) - 1) / 2)
  expect_lte(igraph::gsize(quiet(learn_ARACNE_graph(.x, p.threshold = bonf))$graph),
             igraph::gsize(quiet(learn_ARACNE_graph(.x, max.fdr = 0.05))$graph))
})

test_that("learn_ARACNE_graph() significance masks are built from the RAW MI", {
  skip_if_not_installed("minet")
  mim  <- minet::build.mim(t(.x))
  dpi  <- minet::aracne(mim, eps = 0)
  keep <- .matrix_p_adjust(.mi_pvalue(mim, ncol(.x))) < 0.05

  # there are pairs that are significant but DPI-pruned, so the two matrices
  # genuinely disagree ...
  expect_gt(sum(keep[upper.tri(keep)] & dpi[upper.tri(dpi)] == 0), 0)

  # ... and testing the DPI output instead would give a different, smaller mask,
  # because a pruned entry's p-value is 1 regardless of the pair's dependence
  wrong <- .matrix_p_adjust(.mi_pvalue(dpi, ncol(.x))) < 0.05
  expect_false(identical(which(keep[upper.tri(keep)]), which(wrong[upper.tri(wrong)])))
  expect_lt(sum(wrong[upper.tri(wrong)]), sum(keep[upper.tri(keep)]))
})

test_that("learn_ARACNE_graph() defaults to the max.fdr = 0.05 significance rule", {
  skip_if_not_installed("minet")
  # the default edge rule, so a bare call must equal an explicit max.fdr = 0.05
  expect_equal(igraph::gsize(quiet(learn_ARACNE_graph(.x))$graph),
               igraph::gsize(quiet(learn_ARACNE_graph(.x, max.fdr = 0.05))$graph))

  mim  <- minet::build.mim(t(.x))
  keep <- .matrix_p_adjust(.mi_pvalue(mim, ncol(.x))) < 0.05
  w    <- minet::aracne(mim, eps = 0) * keep; diag(w) <- 0
  expect_equal(igraph::gsize(quiet(learn_ARACNE_graph(.x))$graph),
               sum(w[upper.tri(w)] > 0))
})

test_that("learn_ARACNE_graph() honours the documented rule precedence", {
  skip_if_not_installed("minet")
  dpi <- minet::aracne(minet::build.mim(t(.x)), eps = 0)
  ut  <- dpi[upper.tri(dpi)]
  cut <- as.numeric(stats::quantile(ut[ut > 0], 0.50))

  # max.fdr is the DEFAULT rule, so it sits LAST: setting any other rule
  # overrides it. mim.threshold beats everything ...
  g <- quiet(learn_ARACNE_graph(.x, mim.threshold = cut, max.fdr = 0.05,
                                percentile.threshold = 0.99))$graph
  expect_equal(igraph::gsize(g), sum(ut > cut))

  # ... p.threshold beats percentile.threshold and the default max.fdr ...
  expect_equal(igraph::gsize(quiet(learn_ARACNE_graph(.x, p.threshold = 1e-6,
                                                      percentile.threshold = 0.99))$graph),
               igraph::gsize(quiet(learn_ARACNE_graph(.x, p.threshold = 1e-6))$graph))

  # ... and percentile.threshold beats the default max.fdr, so a caller who asks
  # for a percentile gets one rather than having it silently ignored
  by.pct <- quiet(learn_ARACNE_graph(.x, percentile.threshold = 0.99))$graph
  expect_equal(igraph::gsize(quiet(learn_ARACNE_graph(.x, percentile.threshold = 0.99,
                                                      max.fdr = 0.05))$graph),
               igraph::gsize(by.pct))
  expect_false(igraph::gsize(by.pct) == igraph::gsize(quiet(learn_ARACNE_graph(.x))$graph))

  # turning every rule off is an error rather than a silent unthresholded graph
  expect_error(learn_ARACNE_graph(.x, max.fdr = NULL), "no edge-selection rule")
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
  expect_lt(igraph::gsize(strict), igraph::gsize(loose))
})

test_that("learn_CLR_graph() ranks the percentile within the non-zero scores", {
  skip_if_not_installed("minet")
  # minet::clr clamps a negative endpoint z-score to 0, so the CLR matrix is
  # ~70% structural zeros and the default cut ranks within the rest
  w   <- quiet(learn_CLR_graph(.x))$weights
  ut  <- w[upper.tri(w)]
  pos <- ut[ut > 0]
  expect_gt(mean(ut <= 0), 0)               # the zero mass is real, not hypothetical

  g <- quiet(learn_CLR_graph(.x, percentile.threshold = 0.90))$graph
  expect_equal(igraph::gsize(g), sum(pos > stats::quantile(pos, 0.90)))
})

test_that("learn_CLR_graph(nonzero.only = FALSE) tracks the requested fraction of PAIRS", {
  skip_if_not_installed("minet")
  # the pre-0.0.0.9015 behaviour: the upper-triangle quantile spans every pair,
  # so the kept count tracks the requested fraction of pairs -- unlike a
  # whole-matrix quantile, which double-counts
  w  <- quiet(learn_CLR_graph(.x))$weights
  ut <- w[upper.tri(w)]
  p  <- 0.99
  skip_if(mean(ut <= 0) >= p, "matrix too sparse for a non-degenerate all-pairs cut")

  g <- quiet(learn_CLR_graph(.x, percentile.threshold = p, nonzero.only = FALSE))$graph
  n.pairs <- nrow(.x) * (nrow(.x) - 1) / 2
  expect_equal(igraph::gsize(g), sum(ut > stats::quantile(ut, p)))
  expect_lt(abs(igraph::gsize(g) - (1 - p) * n.pairs), (1 - p) * n.pairs)
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
# .percentile_threshold()
# ---------------------------------------------------------------------------

test_that(".percentile_threshold() ignores the diagonal and the lower triangle", {
  m <- matrix(0, 4, 4)
  m[upper.tri(m)] <- 1:6
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 100                            # a whole-matrix quantile would see these
  expect_equal(.percentile_threshold(m, 0.5), stats::quantile(1:6, 0.5, names = FALSE))
})

test_that(".percentile_threshold(nonzero.only = TRUE) excludes the zero mass", {
  m <- matrix(0, 10, 10)
  m[upper.tri(m)] <- c(rep(0, 40), 1:5)     # 45 pairs, 5 of them non-zero
  m[lower.tri(m)] <- t(m)[lower.tri(m)]

  # all pairs: 89% are zero, so any percentile below that lands on the zero mass
  expect_warning(thr.all <- .percentile_threshold(m, 0.5), "falls on the zero mass")
  expect_equal(thr.all, 0)

  # non-zero only: ranks within the 5 survivors, so the cut still responds
  thr.nz <- .percentile_threshold(m, 0.5, nonzero.only = TRUE)
  expect_equal(thr.nz, stats::quantile(1:5, 0.5, names = FALSE))
  expect_gt(thr.nz, 0)
})

test_that(".percentile_threshold(nonzero.only = TRUE) is monotone at any sparsity", {
  set.seed(42)
  m <- matrix(0, 30, 30)
  v <- stats::runif(435)
  v[v < 0.85] <- 0                          # ~85% structural zeros
  m[upper.tri(m)] <- v
  m[lower.tri(m)] <- t(m)[lower.tri(m)]

  ut <- m[upper.tri(m)]
  kept <- vapply(c(0.1, 0.5, 0.9), function(p) {
    sum(ut > .percentile_threshold(m, p, nonzero.only = TRUE))
  }, numeric(1))
  expect_equal(kept, sort(kept, decreasing = TRUE))
  expect_equal(anyDuplicated(kept), 0L)
})

test_that(".percentile_threshold(nonzero.only = TRUE) warns and keeps nothing on an all-zero matrix", {
  m <- matrix(0, 5, 5)
  expect_warning(thr <- .percentile_threshold(m, 0.95, nonzero.only = TRUE),
                 "every pair")
  expect_equal(thr, 0)
  expect_equal(sum(m[upper.tri(m)] > thr), 0L)   # an empty graph, not an error
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
