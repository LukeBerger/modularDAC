# Tests for the module detection functions in R/01.DetectModules.R
# Scope: true_modules(), true_fuzzy(), find_WGCNA_mods(), find_ICA_mods(),
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
# find_ICA_mods()
# ---------------------------------------------------------------------------

test_that("find_ICA_mods() returns a similarity slot and two valid module objects", {
  skip_if_not_installed("fastICA")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(find_ICA_mods(fx$x, n.comp = 6)))

  expect_type(w, "list")
  expect_named(w, c("similarity", "initial.mods", "final.mods"))

  # no max.size constraint -> no network learned, nothing traded
  expect_null(w$similarity)

  # initial.mods / final.mods: valid non-overlapping module objects
  expect_s4_class(w$initial.mods, "module")
  expect_s4_class(w$final.mods, "module")
  expect_false(w$initial.mods@overlapping)
  expect_true(.module_check(fx$x, w$initial.mods))
  expect_true(.module_check(fx$x, w$final.mods))

  # every feature is assigned to one of the n.comp components (none unassigned)
  expect_equal(length(w$initial.mods@index.vector), nrow(fx$x))
  expect_true(all(w$initial.mods@index.vector %in% 1:6))
  expect_length(w$initial.mods@score.vector, nrow(fx$x))

  # with no constraint the final modules equal the initial modules
  expect_equal(w$final.mods@index.vector, w$initial.mods@index.vector)
})

test_that("find_ICA_mods() respects max.size via eigengene trading without a network", {
  skip_if_not_installed("fastICA")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 6, max.size = 25, trade.by = "eigengene")))

  # eigengene trading needs no similarity matrix
  expect_null(w$similarity)
  # 6 modules over 120 features leave room, so trading meets max.size
  expect_lte(max(lengths(w$final.mods@index.list)), 25)
  expect_true(.module_check(fx$x, w$final.mods))
})

test_that("find_ICA_mods() learns a WGCNA similarity matrix for adjacency trading", {
  skip_if_not_installed("fastICA")
  skip_if_not_installed("WGCNA")
  fx <- make_lfr_fixture()

  # 6 components over 120 features => the largest module is always >= 20 > 18,
  # so adjacency trading always runs and a similarity matrix is learned
  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 6, max.size = 18, trade.by = "adjacency")))

  expect_true(is.matrix(w$similarity))
  expect_equal(dim(w$similarity), c(nrow(fx$x), nrow(fx$x)))
  expect_equal(rownames(w$similarity), rownames(fx$x))
  # trading never grows the largest module
  expect_lte(max(lengths(w$final.mods@index.list)),
             max(lengths(w$initial.mods@index.list)))
  expect_true(.module_check(fx$x, w$final.mods))
})

test_that("find_ICA_mods() splits oversized modules when iterate = TRUE", {
  skip_if_not_installed("fastICA")
  skip_if_not_installed("WGCNA")
  fx <- make_lfr_fixture()

  # 2 components over 120 features => each module starts well above max.size,
  # so iterative ICA must split them into more modules
  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 2, max.size = 40, iterate = TRUE)))

  expect_gt(length(w$final.mods@index.list), 2)
  expect_lte(max(lengths(w$final.mods@index.list)), 40)
  expect_true(.module_check(fx$x, w$final.mods))
})

test_that("find_ICA_mods() merge = TRUE stays valid and never adds modules", {
  skip_if_not_installed("fastICA")
  skip_if_not_installed("WGCNA")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 6, merge = TRUE)))

  # merging can only reduce (or keep) the number of modules
  expect_lte(length(w$initial.mods@index.list), 6)
  expect_true(.module_check(fx$x, w$initial.mods))
  expect_true(.module_check(fx$x, w$final.mods))
})

test_that("find_ICA_mods() can learn an ARACNE similarity matrix", {
  skip_if_not_installed("fastICA")
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 6, max.size = 18,
                  trade.by = "adjacency", network = "ARACNE")))

  expect_true(is.matrix(w$similarity))
  expect_equal(dim(w$similarity), c(nrow(fx$x), nrow(fx$x)))
  expect_true(.module_check(fx$x, w$final.mods))
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

# ---------------------------------------------------------------------------
# module_match()
# ---------------------------------------------------------------------------

test_that("module_match() returns an overlap matrix, best-match table, and overall score", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_match(tm, tm)

  expect_type(res, "list")
  expect_named(res, c("overlap.matrix", "best.matches", "overall.overlap"))

  n <- length(tm@index.list)
  # overlap matrix: n1 x n2, labelled by module
  expect_equal(dim(res$overlap.matrix), c(n, n))
  expect_equal(rownames(res$overlap.matrix), names(tm@index.list))
  expect_equal(colnames(res$overlap.matrix), names(tm@index.list))

  # best-match table: three columns, one row per matched pair
  expect_named(res$best.matches, c("set1", "set2", "percent"))
  expect_equal(nrow(res$best.matches), n)

  # overall score is a single number in 0..100
  expect_length(res$overall.overlap, 1)
  expect_gte(res$overall.overlap, 0)
  expect_lte(res$overall.overlap, 100)
})

test_that("module_match() scores identical module sets as a perfect match", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_match(tm, tm)

  # diagonal of the overlap matrix equals each module's size (modules disjoint)
  expect_equal(unname(diag(res$overlap.matrix)), unname(lengths(tm@index.list)))
  # every node is matched, so each best pair and the overall score are 100%
  expect_true(all(res$best.matches$percent == 100))
  expect_equal(res$overall.overlap, 100)
})

test_that("module_match() overlap matrix counts shared nodes", {
  fx <- make_modular_fixture()
  m1 <- true_modules(fx$g)
  m2 <- true_fuzzy(m1, fx$g)        # overlapping superset of m1, same dataset

  res <- module_match(m1, m2)

  expect_equal(dim(res$overlap.matrix),
               c(length(m1@index.list), length(m2@index.list)))
  # every cell equals the intersection of the two index lists
  for (i in seq_along(m1@index.list)) {
    for (j in seq_along(m2@index.list)) {
      expect_equal(res$overlap.matrix[i, j],
                   length(intersect(m1@index.list[[i]], m2@index.list[[j]])))
    }
  }
})

test_that("module_match() is invariant to module ordering", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  shuffled <- tm
  ord <- rev(seq_along(tm@index.list))
  shuffled@index.list <- tm@index.list[ord]
  shuffled@name.list  <- tm@name.list[ord]

  res <- module_match(tm, shuffled)
  # same partition, just reordered -> still a perfect match
  expect_equal(res$overall.overlap, 100)
  expect_true(all(res$best.matches$percent == 100))
})

test_that("module_match() scores a coarser partition between 0 and 100", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)
  skip_if(length(tm@index.list) < 2)

  # merge the first two modules into one -> partial agreement, one module of tm
  # is then left without a partner
  merged <- tm
  merged@index.list <- c(list(sort(unlist(tm@index.list[1:2]))),
                         tm@index.list[-(1:2)])
  merged@name.list  <- c(list(unlist(tm@name.list[1:2])),
                         tm@name.list[-(1:2)])

  res <- module_match(tm, merged)
  expect_gt(res$overall.overlap, 0)
  expect_lt(res$overall.overlap, 100)
})

test_that("module_match() rejects non-module inputs", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  expect_error(module_match(tm, 42), regexp = "module")
  expect_error(module_match("x", tm), regexp = "module")
})

test_that("module_match() rejects module sets from different datasets", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)
  nm <- rownames(fx$x)

  # a valid module covering only half the nodes -> different index universe
  half <- seq_len(nrow(fx$x) %/% 2)
  m.half <- methods::new("module",
                         source = "half",
                         overlapping = FALSE,
                         index.list = list(half),
                         name.list  = list(nm[half]))
  expect_error(module_match(tm, m.half), regexp = "different node")

  # same indices but a renamed node -> different name universe
  renamed <- tm
  renamed@name.list[[1]][1] <- "not_a_real_feature"
  expect_error(module_match(tm, renamed), regexp = "different node names")
})

# ---------------------------------------------------------------------------
# module_contiguity()
# ---------------------------------------------------------------------------

test_that("module_contiguity() returns an overall score and a per-module edge table", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_contiguity(fx$g, tm)

  expect_type(res, "list")
  expect_named(res, c("overall.contiguity", "module.edges"))

  # overall score: single number in 0..100
  expect_length(res$overall.contiguity, 1)
  expect_gte(res$overall.contiguity, 0)
  expect_lte(res$overall.contiguity, 100)

  # per-module table: one row per module, three named columns, labelled rows
  expect_s3_class(res$module.edges, "data.frame")
  expect_named(res$module.edges, c("within.edges", "between.edges", "percent.within"))
  expect_equal(nrow(res$module.edges), length(tm@index.list))
  expect_equal(rownames(res$module.edges), names(tm@index.list))
})

test_that("module_contiguity() edge counts match a brute-force adjacency calculation", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_contiguity(fx$g, tm)
  adj <- as.matrix(igraph::as_adjacency_matrix(fx$g))

  for (i in seq_along(tm@name.list)) {
    idx <- which(rownames(adj) %in% tm@name.list[[i]])
    within  <- sum(adj[idx, idx, drop = FALSE]) / 2          # both endpoints inside
    between <- sum(adj[idx, -idx, drop = FALSE])             # exactly one inside
    expect_equal(res$module.edges$within.edges[i], within)
    expect_equal(res$module.edges$between.edges[i], between)
  }
})

test_that("module_contiguity() percentages and overall score follow their definitions", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  res <- module_contiguity(fx$g, tm)
  w <- res$module.edges$within.edges
  b <- res$module.edges$between.edges

  expect_equal(res$module.edges$percent.within, 100 * w / (w + b))
  expect_equal(res$overall.contiguity, 100 * sum(w) / sum(w + b))
})

test_that("module_contiguity() reports a single all-encompassing module as fully contiguous", {
  fx <- make_modular_fixture()

  # one module containing every vertex -> no edges leave it
  m.all <- methods::new("module",
                        source = "all",
                        overlapping = FALSE,
                        index.list = list(seq_len(igraph::vcount(fx$g))),
                        name.list  = list(igraph::V(fx$g)$name))

  res <- module_contiguity(fx$g, m.all)
  expect_equal(res$module.edges$within.edges, igraph::gsize(fx$g))
  expect_equal(res$module.edges$between.edges, 0)
  expect_equal(res$module.edges$percent.within, 100)
  expect_equal(res$overall.contiguity, 100)
})

test_that("module_contiguity() handles overlapping modules", {
  fx <- make_modular_fixture()
  tf <- true_fuzzy(true_modules(fx$g), fx$g)   # overlapping modules

  res <- module_contiguity(fx$g, tf)
  expect_equal(nrow(res$module.edges), length(tf@index.list))
  expect_named(res$module.edges, c("within.edges", "between.edges", "percent.within"))
  expect_gte(res$overall.contiguity, 0)
  expect_lte(res$overall.contiguity, 100)
})

test_that("module_contiguity() validates its inputs", {
  fx <- make_modular_fixture()
  tm <- true_modules(fx$g)

  expect_error(module_contiguity("not a graph", tm), regexp = "igraph")
  expect_error(module_contiguity(fx$g, 42), regexp = "module")

  # a module node name absent from the graph is rejected
  bad <- tm
  bad@name.list[[1]][1] <- "not_a_real_node"
  expect_error(module_contiguity(fx$g, bad), regexp = "not present")
})
