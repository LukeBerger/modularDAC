# Tests for the module detection functions in R/01.DetectModules.R
# Scope: true_modules(), true_fuzzy(), find_WGCNA_mods(),
#        eigen_fuzzy_modules(), adj_fuzzy_modules()
#
# Module objects are validated with the internal .module_check(), which enforces
# the structural contract of the "module" S4 class (matching index/name lists,
# correct overlap behaviour, feature names drawn from the data, etc.).
#
# Functions are stochastic, so every fixture seeds the RNG for reproducibility.

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

# An Erdos-Renyi modular graph + data: cheap, used for the index-based
# (true_*, eigen_*, adj_*) functions that only need consistent node ordering.
make_modular_fixture <- function(seed = 1) {
  set.seed(seed)
  g <- make_modular_graph()
  list(g = g, x = sim_graph_data(g, n.samples = 100))
}

# An LFR graph + data: has the community structure WGCNA needs to actually
# recover modules. find_WGCNA_mods() is the only function that requires it.
make_lfr_fixture <- function(seed = 1) {
  set.seed(seed)
  g <- make_lfr(n = 120)
  list(g = g, x = sim_graph_data(g, n.samples = 100))
}

# ---------------------------------------------------------------------------
# true_modules()
# ---------------------------------------------------------------------------

test_that("true_modules() builds a valid non-overlapping module object", {
  fx <- make_modular_fixture()

  tm <- true_modules(fx$g)

  expect_s4_class(tm, "module")
  expect_false(tm@overlapping)
  expect_true(.module_check(fx$x, tm))
})

test_that("true_modules() mirrors the graph's ground-truth modules", {
  fx <- make_modular_fixture()

  tm <- true_modules(fx$g)

  # one entry per ground-truth module
  expect_equal(length(tm@index.list),
               length(unique(igraph::V(fx$g)$module)))
  # the membership vector covers every node
  expect_equal(length(tm@index.vector), length(fx$g))
  expect_equal(tm@index.vector, igraph::V(fx$g)$module)
  # names come straight from the graph
  expect_equal(sort(unlist(tm@name.list, use.names = FALSE)),
               sort(igraph::V(fx$g)$name))
})

# ---------------------------------------------------------------------------
# true_fuzzy()
# ---------------------------------------------------------------------------

test_that("true_fuzzy() builds a valid overlapping module object", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  tf <- true_fuzzy(tm, fx$g)

  expect_s4_class(tf, "module")
  expect_true(tf@overlapping)
  expect_true(.module_check(fx$x, tf))
})

test_that("true_fuzzy() expands each module to its graph neighbourhood", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  tf <- true_fuzzy(tm, fx$g)

  # same number of modules, but each is grown to include neighbours, so the
  # fuzzy module is a superset of the original
  expect_equal(length(tf@index.list), length(tm@index.list))
  for (i in seq_along(tm@index.list)) {
    expect_true(all(tm@index.list[[i]] %in% tf@index.list[[i]]))
  }
  expect_true(all(lengths(tf@index.list) >= lengths(tm@index.list)))
})

# ---------------------------------------------------------------------------
# find_WGCNA_mods()
# ---------------------------------------------------------------------------

test_that("find_WGCNA_mods() returns the adjacency matrix and two module objects", {
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(find_WGCNA_mods(fx$x, max.size = 60)))

  expect_type(w, "list")
  expect_named(w, c("wgcna.adj", "initial.mods", "final.mods"))

  # wgcna.adj: square p x p matrix named by feature
  expect_true(is.matrix(w$wgcna.adj))
  expect_equal(dim(w$wgcna.adj), c(nrow(fx$x), nrow(fx$x)))
  expect_equal(rownames(w$wgcna.adj), rownames(fx$x))

  # initial.mods / final.mods: valid non-overlapping module objects
  expect_s4_class(w$initial.mods, "module")
  expect_s4_class(w$final.mods, "module")
  expect_false(w$initial.mods@overlapping)
  expect_false(w$final.mods@overlapping)
  expect_true(.module_check(fx$x, w$initial.mods))
  expect_true(.module_check(fx$x, w$final.mods))
})

test_that("find_WGCNA_mods() respects the max.size constraint after trading", {
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(find_WGCNA_mods(fx$x, max.size = 60)))

  # the membership vector covers every feature
  expect_equal(length(w$final.mods@index.vector), nrow(fx$x))
  # node trading should bring every final module within max.size
  expect_lte(max(lengths(w$final.mods@index.list)), 60)
})

# ---------------------------------------------------------------------------
# eigen_fuzzy_modules()
# ---------------------------------------------------------------------------

test_that("eigen_fuzzy_modules() builds a valid overlapping module object", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  ef <- eigen_fuzzy_modules(fx$x, tm, max.size = 60, ratio = 1.5)

  expect_s4_class(ef, "module")
  expect_true(ef@overlapping)
  expect_true(.module_check(fx$x, ef))
})

test_that("eigen_fuzzy_modules() grows modules by ratio up to max.size", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  ef <- eigen_fuzzy_modules(fx$x, tm, max.size = 60, ratio = 1.5)

  # fuzzy modules are supersets of the originals, capped at max.size
  for (i in seq_along(tm@index.list)) {
    expect_true(all(tm@index.list[[i]] %in% ef@index.list[[i]]))
  }
  expect_lte(max(lengths(ef@index.list)), 60)
  expect_true(all(lengths(ef@index.list) >= lengths(tm@index.list)))
})

test_that("eigen_fuzzy_modules() errors when a module exceeds max.size", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  # default modules have 40 nodes each; max.size = 5 is impossible
  expect_error(eigen_fuzzy_modules(fx$x, tm, max.size = 5))
})

# ---------------------------------------------------------------------------
# adj_fuzzy_modules()
# ---------------------------------------------------------------------------

test_that("adj_fuzzy_modules() builds a valid overlapping module object", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)
  adj <- abs(stats::cor(t(fx$x)))

  af <- adj_fuzzy_modules(fx$x, adj, tm, max.size = 60, ratio = 1.5)

  expect_s4_class(af, "module")
  expect_true(af@overlapping)
  expect_true(.module_check(fx$x, af))
})

test_that("adj_fuzzy_modules() grows modules by ratio up to max.size", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)
  adj <- abs(stats::cor(t(fx$x)))

  af <- adj_fuzzy_modules(fx$x, adj, tm, max.size = 60, ratio = 1.5)

  for (i in seq_along(tm@index.list)) {
    expect_true(all(tm@index.list[[i]] %in% af@index.list[[i]]))
  }
  expect_lte(max(lengths(af@index.list)), 60)
  expect_true(all(lengths(af@index.list) >= lengths(tm@index.list)))
})

test_that("adj_fuzzy_modules() errors when a module exceeds max.size", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)
  adj <- abs(stats::cor(t(fx$x)))

  expect_error(adj_fuzzy_modules(fx$x, adj, tm, max.size = 5, ratio = 1.5))
})

# ---------------------------------------------------------------------------
# module_correlation()
# ---------------------------------------------------------------------------

test_that("module_correlation() returns per-module upper-tri matrices and a ggplot", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_correlation(fx$x, tm)

  # list with the documented structure
  expect_type(res, "list")
  expect_named(res, c("correlation", "plot"))
  expect_true(inherits(res$plot, "ggplot"))

  # one correlation matrix per module, carrying the module labels
  expect_equal(length(res$correlation), length(tm@index.list))
  expect_equal(names(res$correlation), names(tm@index.list))

  # each matrix is square, named by feature, with only the upper triangle filled
  for (i in seq_along(tm@index.list)) {
    idx <- tm@index.list[[i]]
    cm  <- res$correlation[[i]]
    expect_equal(dim(cm), c(length(idx), length(idx)))
    expect_equal(rownames(cm), rownames(fx$x)[idx])
    expect_equal(colnames(cm), rownames(fx$x)[idx])
    # diagonal and lower triangle are blanked, upper triangle is populated
    expect_true(all(is.na(cm[lower.tri(cm, diag = TRUE)])))
    expect_false(any(is.na(cm[upper.tri(cm)])))
  }
})

test_that("module_correlation() reproduces stats::cor on a module's features", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_correlation(fx$x, tm)

  # the first module's matrix should equal a direct correlation, with the
  # diagonal and lower triangle masked out
  idx <- tm@index.list[[1]]
  expected <- stats::cor(t(fx$x[idx, , drop = FALSE]))
  expected[lower.tri(expected, diag = TRUE)] <- NA
  expect_equal(res$correlation[[1]], expected)
})

test_that("module_correlation() passes cor.method through to stats::cor", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_correlation(fx$x, tm, cor.method = "spearman")

  idx <- tm@index.list[[1]]
  expected <- stats::cor(t(fx$x[idx, , drop = FALSE]), method = "spearman")
  expected[lower.tri(expected, diag = TRUE)] <- NA
  expect_equal(res$correlation[[1]], expected)
})

test_that("module_correlation() boxplot pools every upper-tri correlation", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_correlation(fx$x, tm)

  # one row per within-module feature pair, summed over all modules
  n.pairs <- sum(vapply(tm@index.list,
                        function(m) choose(length(m), 2),
                        numeric(1)))
  expect_equal(nrow(res$plot$data), n.pairs)
  expect_true(all(c("module", "correlation") %in% names(res$plot$data)))
})

test_that("module_correlation() handles single-feature modules", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  nm <- rownames(fx$x)

  # an overlapping pair where the second module holds a single feature shared
  # with the first (so the overlap contract in .module_check is satisfied)
  m <- methods::new("module",
                    source = "single-feature",
                    overlapping = TRUE,
                    index.list = list(A = c(1, 2, 3), B = c(3)),
                    name.list  = list(A = nm[c(1, 2, 3)], B = nm[3]))

  res <- module_correlation(fx$x, m)

  # the single-feature module yields a 1x1 NA matrix and contributes no pairs
  expect_equal(dim(res$correlation[["B"]]), c(1, 1))
  expect_true(is.na(res$correlation[["B"]][1, 1]))
  expect_true(inherits(res$plot, "ggplot"))
  expect_equal(nrow(res$plot$data), choose(3, 2))
})

test_that("module_correlation() rejects a module that does not match the data", {
  skip_if_not_installed("ggplot2")
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  # corrupt the feature names so .module_check's data consistency check fails
  bad <- tm
  bad@name.list[[1]][1] <- "not_a_real_feature"
  expect_error(module_correlation(fx$x, bad))
})
