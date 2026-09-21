# Tests for module_correlation() in R/module_correlation.R

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
