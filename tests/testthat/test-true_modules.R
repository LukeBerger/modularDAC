# Tests for true_modules() in R/true_modules.R

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
