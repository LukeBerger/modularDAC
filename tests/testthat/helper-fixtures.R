# Shared fixtures for the modularDAC test suite.
#
# testthat sources every helper-*.R file before running the tests, so everything
# here is available to all test files. The graph and data matrix are built once
# here rather than per file, since several test files need the same ones.
#
# Module objects are validated in the tests with the internal .module_check(),
# which enforces the structural contract of the "module" S4 class (matching
# index/name lists, correct overlap behaviour, feature names drawn from the
# data). Generators are stochastic, so every fixture seeds the RNG.

# Silence cat()/message()/warning() noise while returning the value. SILGGM and
# WGCNA report progress via cat() rather than message(), so suppressMessages()
# alone does not keep the test log clean.
quiet <- function(expr) {
  utils::capture.output(val <- suppressWarnings(suppressMessages(expr)))
  val
}

# An Erdos-Renyi modular graph + data: cheap, used for the index-based
# (true_*, eigen_*, adj_*) functions that only need consistent node ordering.
make_modular_fixture <- function(seed = 1) {
  set.seed(seed)
  g <- make_modular_graph()
  list(g = g, x = sim_graph_data(g, n.samples = 100))
}

# An LFR graph + data: has the community structure WGCNA needs to actually
# recover modules. find_WGCNA_mods() is the only function that requires it.
make_lfr_fixture <- function(seed = 1) {
  set.seed(seed)
  g <- make_lfr(n = 120)
  list(g = g, x = sim_graph_data(g, n.samples = 100))
}

# The graph, data matrix and overlapping modules shared by the graph-learning
# and divide-and-conquer test files.
set.seed(1)
.g <- make_modular_graph()
.x <- sim_graph_data(.g, n.samples = 100)   # p x n (features x samples)
.fuzzy <- true_fuzzy(true_modules(.g), .g)  # overlapping modules to divide along
