# Tests for learn_ARACNE_graph() and its private .mi_pvalue() helper in
# R/learn_ARACNE_graph.R

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
