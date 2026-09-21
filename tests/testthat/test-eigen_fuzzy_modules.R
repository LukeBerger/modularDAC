# Tests for eigen_fuzzy_modules() in R/eigen_fuzzy_modules.R

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
