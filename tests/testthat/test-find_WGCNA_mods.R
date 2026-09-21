# Tests for find_WGCNA_mods() and its private .relabel_sequential()
# helper in R/find_WGCNA_mods.R

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

test_that("find_WGCNA_mods() labels the final modules 1:N with no gaps", {
  fx <- make_lfr_fixture()

  # a tight max.size forces heavy iterative splitting/merging, the path that
  # used to leave inflated, gappy ids (e.g. 1, 45, 101, 4)
  w <- suppressWarnings(suppressMessages(
    find_WGCNA_mods(fx$x, min.size = 5, max.size = 20)
  ))
  fin <- w$final.mods
  n <- length(fin@index.list)

  # both the labels and the membership ids form a contiguous 1:N sequence
  expect_setequal(as.integer(names(fin@index.list)), seq_len(n))
  expect_setequal(unique(fin@index.vector), seq_len(n))
  # relabelling preserves a valid partition
  expect_true(.module_check(fx$x, fin))
})

test_that(".relabel_sequential() renumbers modules to 1:N and preserves 0", {
  # gappy / inflated ids collapse to a contiguous run, ascending by original id
  expect_equal(.relabel_sequential(c(1, 45, 45, 101, 4)),
               c(1,  3,  3,   4, 2))
  # the unassigned group (0) is left untouched
  expect_equal(.relabel_sequential(c(0, 7, 7, 30, 0, 4)),
               c(0, 2, 2,  3, 0, 1))
  # the grouping itself is unchanged, only the labels
  iv <- c(3, 3, 9, 1, 9, 1)
  expect_equal(unname(split(seq_along(iv), .relabel_sequential(iv))),
               unname(split(seq_along(iv), iv)))
})
