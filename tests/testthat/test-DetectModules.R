test_that("Module Check", {
  # Simulate True Graph and Data
  g <- make_lfr(n = 360 )
  x <- sim_graph_data(g, n.samples = 180)

  # true mods
  t <- true_modules(g)
  expect_true(.module_check(x , t))

  # basic mods
  w <- suppressMessages(find_WGCNA_mods(x))
  expect_true(.module_check(x , w)) # commenting out to save on run time when not working with this function

  i <- find_ICA_mods(x, 10)
  expect_true(.module_check(x , i))

  p <- pragmatic_modules(x,n.mods = 10, max.size = 60)
  expect_true(.module_check(x , p))

  # fuzzy mods
  tf <- true_fuzzy(t, g)
  expect_true(.module_check(x , tf))

  ef <- eigen_fuzzy_modules(x, p, 100)
  expect_true(.module_check(x , ef))

  # nf <- nodewise_fuzzy_modules(x, p, 100)
  # expect_true(.module_check(x , nf))

  # overlap mods
  # o <- create_overlap_modules(x, p)
  # expect_true(.module_check(x , o))
  #
  # # complex ICA mods
  # ci <- complex_ICA_modules(x, 3)
  # expect_true(.module_check(x , ci))

  # assess accuracy and contiguity of modules
  pm <- percent_module_match(t, i)
  expect_type(pm, "double")
  expect_length(pm, 1)
  expect_gte(pm, 0)
  expect_lte(pm, 100)
  expect_equal(percent_module_match(t, t), 100)

  # module contiguity
  mc <- module_contiguity(g, t)
  expect_type(mc, "double")
  expect_length(mc, 1)
  expect_gte(mc, 0)
  expect_lte(mc, 100)

})
