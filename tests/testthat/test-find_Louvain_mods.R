# Tests for find_Louvain_mods() in R/find_Louvain_mods.R

test_that("find_Louvain_mods() returns the graph, its adjacency and two module objects", {
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  set.seed(1)
  l <- suppressWarnings(suppressMessages(find_Louvain_mods(fx$x, max.size = 60)))

  expect_type(l, "list")
  expect_named(l, c("graph", "graph.adj", "initial.mods", "final.mods"))
  expect_s3_class(l$graph, "igraph")

  # graph.adj: symmetric, non-negative p x p matrix named by feature
  expect_true(is.matrix(l$graph.adj))
  expect_equal(dim(l$graph.adj), c(nrow(fx$x), nrow(fx$x)))
  expect_equal(rownames(l$graph.adj), rownames(fx$x))
  expect_equal(colnames(l$graph.adj), rownames(fx$x))
  expect_true(isSymmetric(l$graph.adj))
  expect_true(all(l$graph.adj >= 0))
  # one non-zero pair of entries per edge
  expect_equal(sum(l$graph.adj > 0), 2 * igraph::ecount(l$graph))

  # initial.mods / final.mods: valid non-overlapping module objects
  expect_s4_class(l$initial.mods, "module")
  expect_s4_class(l$final.mods, "module")
  expect_false(l$initial.mods@overlapping)
  expect_false(l$final.mods@overlapping)
  expect_true(.module_check(fx$x, l$initial.mods))
  expect_true(.module_check(fx$x, l$final.mods))
})

test_that("find_Louvain_mods() respects min.size and max.size", {
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  set.seed(1)
  l <- suppressMessages(find_Louvain_mods(fx$x, min.size = 10, max.size = 30))
  sizes <- lengths(l$final.mods@index.list)

  # every feature is in exactly one module
  expect_equal(length(l$final.mods@index.vector), nrow(fx$x))
  expect_equal(sum(sizes), nrow(fx$x))
  expect_lte(max(sizes), 30)
  expect_gte(min(sizes), 10)
})

test_that("find_Louvain_mods() labels the final modules 1:N with no gaps", {
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  # a tight max.size forces splitting, and merging hands out and deletes ids
  set.seed(1)
  l <- suppressWarnings(suppressMessages(find_Louvain_mods(fx$x, min.size = 5, max.size = 15)))
  fin <- l$final.mods
  n <- length(fin@index.list)

  expect_setequal(as.integer(names(fin@index.list)), seq_len(n))
  expect_setequal(unique(fin@index.vector), seq_len(n))
  expect_true(.module_check(fx$x, fin))
})

test_that("find_Louvain_mods() is reproducible under a fixed seed", {
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  set.seed(7)
  a <- suppressWarnings(suppressMessages(find_Louvain_mods(fx$x, max.size = 60)))
  set.seed(7)
  b <- suppressWarnings(suppressMessages(find_Louvain_mods(fx$x, max.size = 60)))

  expect_identical(a$final.mods@index.vector, b$final.mods@index.vector)
})

test_that("find_Louvain_mods() accepts a supplied graph and places features absent from it", {
  fx <- make_lfr_fixture()

  # drop 10 vertices: those features must still be placed in a module
  absent <- rownames(fx$x)[1:10]
  g <- igraph::delete_vertices(fx$g, absent)

  set.seed(1)
  l <- suppressMessages(find_Louvain_mods(fx$x, graph = g, min.size = 10, max.size = 60))

  expect_true(.module_check(fx$x, l$final.mods))
  expect_equal(sum(lengths(l$final.mods@index.list)), nrow(fx$x))
  expect_true(all(absent %in% unlist(l$final.mods@name.list)))
  # the absent features have no edges in graph.adj
  expect_true(all(l$graph.adj[absent, ] == 0))
  # every edge of the supplied graph is in graph.adj, and nothing else
  el <- igraph::as_edgelist(g, names = TRUE)
  expect_true(all(l$graph.adj[el] > 0))
  expect_equal(sum(l$graph.adj > 0), 2 * igraph::ecount(g))
})

test_that("find_Louvain_mods() uses absolute weights, so signed weights are allowed", {
  fx <- make_lfr_fixture()
  g <- fx$g
  set.seed(2)
  igraph::E(g)$weight <- stats::runif(igraph::ecount(g), -1, 1)

  set.seed(1)
  l <- suppressMessages(find_Louvain_mods(fx$x, graph = g, min.size = 10, max.size = 60))

  expect_true(all(l$graph.adj >= 0))
  el <- igraph::as_edgelist(g, names = TRUE)
  expect_equal(unname(l$graph.adj[el]), pmax(abs(igraph::E(g)$weight), 1e-12))
})

test_that("adj_fuzzy_modules() on graph.adj with min.adj > 0 recruits only graph neighbours", {
  fx <- make_lfr_fixture()

  set.seed(1)
  l <- suppressMessages(find_Louvain_mods(fx$x, graph = fx$g, min.size = 10, max.size = 40))
  fz <- adj_fuzzy_modules(fx$x, l$graph.adj, l$final.mods,
                          max.size = 80, ratio = 2, min.adj = 1e-12)

  expect_true(.module_check(fx$x, fz))
  for (i in seq_along(fz@index.list)) {
    core <- fz@core.list[[i]]
    aux  <- setdiff(fz@index.list[[i]], core)
    if (length(aux) == 0) next
    # every auxiliary node has an edge to at least one core node
    expect_true(all(colSums(l$graph.adj[core, aux, drop = FALSE] > 0) > 0))
  }
})

test_that("find_Louvain_mods() rejects bad input", {
  fx <- make_lfr_fixture()
  x <- fx$x

  # unnamed or duplicated features
  x.unnamed <- x; rownames(x.unnamed) <- NULL
  expect_error(find_Louvain_mods(x.unnamed, graph = fx$g), "unique row names")
  x.dup <- x; rownames(x.dup)[2] <- rownames(x.dup)[1]
  expect_error(find_Louvain_mods(x.dup, graph = fx$g), "unique row names")

  # size limits that cannot both hold
  expect_error(find_Louvain_mods(x, graph = fx$g, min.size = 50, max.size = 20), "cannot exceed")

  # graphs that cannot be matched to x
  expect_error(find_Louvain_mods(x, graph = "not a graph"), "igraph object")
  expect_error(find_Louvain_mods(x, graph = igraph::as_directed(fx$g)), "undirected")
  g.unnamed <- igraph::delete_vertex_attr(fx$g, "name")
  expect_error(find_Louvain_mods(x, graph = g.unnamed), "named by feature")
  g.other <- fx$g
  igraph::V(g.other)$name <- paste0("Other_", seq_len(igraph::vcount(g.other)))
  expect_error(suppressMessages(find_Louvain_mods(x, graph = g.other)), "match rownames")

  # a percentile outside [0, 1] when the ARACNE graph is to be learned
  skip_if_not_installed("minet")
  expect_error(suppressMessages(find_Louvain_mods(x, percentile.threshold = 1.5)), "between 0 and 1")
})

test_that("find_Louvain_mods() removes zero-variance features with a warning", {
  fx <- make_lfr_fixture()
  x <- fx$x
  x[3, ] <- 1
  const <- rownames(x)[3]

  set.seed(1)
  expect_warning(
    l <- suppressMessages(find_Louvain_mods(x, graph = fx$g, min.size = 10, max.size = 60)),
    "zero variance"
  )
  # the module object describes the reduced matrix, without the constant feature
  expect_equal(l$final.mods@data.dim, c(nrow(x) - 1, ncol(x)))
  expect_false(const %in% unlist(l$final.mods@name.list))
  expect_false(const %in% rownames(l$graph.adj))
  expect_true(.module_check(x[-3, ], l$final.mods))
})

test_that("find_Louvain_mods() keeps every DPI survivor by default and fewer as the percentile rises", {
  skip_if_not_installed("minet")
  fx <- make_lfr_fixture()

  # the DPI survivors, computed straight from minet
  a <- minet::aracne(minet::build.mim(dataset = t(fx$x)), eps = 0.001)
  diag(a) <- 0
  n.survivors <- sum(a[upper.tri(a)] > 0)

  edges <- vapply(c(0, 0.25, 0.5), function(pct) {
    set.seed(1)
    igraph::ecount(suppressMessages(
      find_Louvain_mods(fx$x, max.size = 60, percentile.threshold = pct))$graph)
  }, numeric(1))

  expect_equal(edges[1], n.survivors)
  expect_true(all(diff(edges) < 0))
  # dropping the weakest quarter keeps about three quarters of the survivors
  expect_equal(edges[2] / n.survivors, 0.75, tolerance = 0.02)
})
