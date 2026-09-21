# Tests for module_contiguity() in R/module_contiguity.R

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
