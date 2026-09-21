# Tests for find_ICA_mods() in R/find_ICA_mods.R

test_that("find_ICA_mods() returns loadings and an overlapping module with owned cores", {
  skip_if_not_installed("fastICA")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(find_ICA_mods(fx$x, n.comp = 6)))

  expect_type(w, "list")
  expect_named(w, c("ica.loadings", "mods"))

  # ica.loadings: p x n.comp metagene loading matrix, feature-named rows
  expect_true(is.matrix(w$ica.loadings))
  expect_equal(dim(w$ica.loadings), c(nrow(fx$x), 6))
  expect_equal(rownames(w$ica.loadings), rownames(fx$x))

  m <- w$mods
  expect_s4_class(m, "module")
  expect_true(m@overlapping)
  expect_true(.module_check(fx$x, m))
  expect_length(m@score.vector, nrow(fx$x))

  # cores partition the features: every feature owned by exactly one component
  expect_false(any(duplicated(unlist(m@core.list))))
  expect_setequal(unlist(m@core.list), seq_len(nrow(fx$x)))

  # multi-membership makes the modules overlap (total memberships exceed p)
  expect_gt(length(unlist(m@index.list)), length(unlist(m@core.list)))
})

test_that("find_ICA_mods() membership.z tunes how much the modules overlap", {
  skip_if_not_installed("fastICA")
  fx <- make_lfr_fixture()

  loose  <- suppressWarnings(suppressMessages(find_ICA_mods(fx$x, n.comp = 6, membership.z = 1)))
  strict <- suppressWarnings(suppressMessages(find_ICA_mods(fx$x, n.comp = 6, membership.z = 3)))

  # a higher loading threshold recruits fewer auxiliary members -> less overlap
  expect_lte(sum(lengths(strict$mods@index.list)),
             sum(lengths(loose$mods@index.list)))
  # the core partition is unaffected by the threshold
  expect_setequal(unlist(strict$mods@core.list), unlist(loose$mods@core.list))
})

test_that("find_ICA_mods() trims auxiliary members to respect max.size", {
  skip_if_not_installed("fastICA")
  fx <- make_lfr_fixture()

  w <- suppressWarnings(suppressMessages(
    find_ICA_mods(fx$x, n.comp = 12, max.size = 20)))

  # cores are never trimmed, so each module is capped at max.size whenever its
  # core set fits; oversized cores (should they occur) are left intact and warned
  sizes <- lengths(w$mods@index.list)
  cores <- lengths(w$mods@core.list)
  expect_true(all(sizes <= pmax(20, cores)))
  expect_true(.module_check(fx$x, w$mods))
})
