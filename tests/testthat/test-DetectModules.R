test_that("Module Check", {
  # Simulate True Graph and Data
  er <- make_modular_graph()
  x <- sim_graph_data(er, n.samples = 100)

  # true mods
  t <- .true_modules(er)
  expect_true(.module_check(x , t))

  # basic mods
  # w <- suppressMessages(find_WGCNA_mods(t(x), cor.FN = "bicor"))
  # expect_true(.module_check(x , w)) # commenting out to save on run time when not working with this funciton

  i <- find_ICA_mods(x, 3)
  expect_true(.module_check(x , i))

  p <- pragmatic_modules(x,n.mods = 3, max.size = 60)
  expect_true(.module_check(x , p))

  # fuzzy mods
  ef <- eigen_fuzzy_modules(x, p, 80)
  expect_true(.module_check(x , ef))

  nf <- nodewise_fuzzy_modules(x, p, 80)
  expect_true(.module_check(x , nf))

  # overlap mods
  o <- create_overlap_modules(x, p)
  expect_true(.module_check(x , o))

  # complex ICA mods
  ci <- complex_ICA_modules(x, 3)
  expect_true(.module_check(x , ci))

  # assess accuracy and contiguity of modules
  pm <- percent_module_match(t, i)
  expect_type(pm, "double")
  expect_length(pm, 1)
  expect_gte(pm, 0)
  expect_lte(pm, 100)
  expect_equal(percent_module_match(t, t), 100)

  # module contiguity
  mc <- module_contiguity(er, t)
  expect_type(mc, "double")
  expect_length(mc, 1)
  expect_gte(mc, 0)
  expect_lte(mc, 100)

})
