# Tests for module_match() in R/module_match.R

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
  shuffled@core.list  <- tm@core.list[ord]   # keep ownership aligned with the reorder

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
  # drop the (now stale) ownership so the hand-built object stays valid
  merged@core.list <- list()

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
