# Tests for adj_fuzzy_modules() in R/adj_fuzzy_modules.R

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
