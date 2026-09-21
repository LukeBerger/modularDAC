# Tests for true_fuzzy() in R/true_fuzzy.R

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
