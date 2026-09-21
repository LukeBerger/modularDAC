# Regenerate data/lfr_example.rda, the example dataset used by the README.
#
# Dev-branch only: it uses the benchmarking simulators (make_lfr,
# sim_graph_data) to build a dataset with a known modular structure. The .rda it
# produces ships on both branches; this script does not.
#
#     Rscript data-raw/make_example_data.R
#
# Seed 1 gives a three-community LFR graph with balanced module sizes
# (44 / 40 / 36); other seeds sometimes collapse to two communities.

devtools::load_all(".", quiet = TRUE)
stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))

set.seed(1)

# 1. a benchmark graph with known community structure
g <- make_lfr(n = 120, min.community = 30, max.community = 50)
stopifnot(length(unique(igraph::V(g)$module)) == 3)

# 2. expression-like data whose covariance follows that graph
x <- sim_graph_data(g, n.samples = 60)

# 3. feature and sample names in the shape a user's own data would take
rownames(x) <- sprintf("Gene_%03d", seq_len(nrow(x)))
colnames(x) <- sprintf("Sample_%02d", seq_len(ncol(x)))

# 4. wrap as a SummarizedExperiment, carrying the simulated module labels so the
#    example can be checked against the structure the data was built from
lfr_example <- SummarizedExperiment::SummarizedExperiment(
  assays  = list(expression = x),
  rowData = S4Vectors::DataFrame(
    gene        = rownames(x),
    true.module = as.integer(igraph::V(g)$module),
    row.names   = rownames(x)
  ),
  colData = S4Vectors::DataFrame(
    sample    = colnames(x),
    row.names = colnames(x)
  ),
  metadata = list(
    source = paste("Simulated with modularDAC::make_lfr(n = 120,",
                   "min.community = 30, max.community = 50) and",
                   "modularDAC::sim_graph_data(n.samples = 60), seed 1."),
    n.true.modules = 3L,
    n.true.edges   = igraph::ecount(g)
  )
)

dir.create("data", showWarnings = FALSE)
save(lfr_example, file = "data/lfr_example.rda", compress = "xz")

cat("wrote data/lfr_example.rda:", nrow(x), "features x", ncol(x), "samples,",
    igraph::ecount(g), "true edges,",
    round(file.size("data/lfr_example.rda") / 1024), "KB\n")
