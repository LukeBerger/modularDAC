# modularDAC

Learn large networks by dividing the data into overlapping modules, learning the
modular subnetworks independently, and stitching the results back into a single graph.

---

## The divide-and-conquer algorithm

Network inference methods that estimate partial correlation have to
reason about all `p * (p - 1) / 2` feature pairs at once. Cost grows steeply
with the number of features, and once `p` approaches or exceeds the number of
samples `n` the whole-data estimate is also poorly conditioned: there is not
enough data to pin down that many parameters simultaneously.

modularDAC exploits the fact that biological networks are **modular**. Features
fall into correlated groups, and primarily have edges within their group. If that
holds, the network can be recovered from many small problems instead of one large one.

The algorithm runs in three stages.

**1. Divide.** Partition the features into modules, then *grow* each module into
an overlapping one.

The partition comes from `find_WGCNA_mods()` (co-expression modules via WGCNA)
or `find_ICA_mods()` (independent components, which overlap naturally). The
growth step — `eigen_fuzzy_modules()` or `adj_fuzzy_modules()` — recruits the
features outside a module that correlate most strongly with it, up to a size or
ratio cap.

The overlap is the part that matters, and it is why this is not just clustering
followed by independent fits. Each module **owns** the features of the original
partition; those are its *core* nodes. The recruited features are *auxiliary*:
they are not owned by this module, and are present only so the learner can
condition on them. A partial correlation is trustworthy only when the
conditioning set contains the node's Markov blanket. Cut a module at its
boundary and the features near that boundary get estimated without the
neighbours that explain their correlations, which manufactures edges that are
not there. Growing the module approximates blanket closure for its core nodes.

**2. Conquer.** Learn a graph within each module independently, optionally in
parallel across cores. Any learner can be plugged in; the default is
`learn_SILGGM_graph()` (partial correlation). Also available:
`learn_WGCNA_graph()`, `learn_ARACNE_graph()`, `learn_CLR_graph()`,
`learn_GENIE3_graph()` and `learn_bdgraph_graph()`.

**3. Combine.** Stitch the sub-graphs into one graph over all features.
Ownership decides who gets a say: an edge is credited only from a module in
which at least one endpoint is a core node, so every edge has at most two
authoritative proposers (the owners of its two endpoints), and
auxiliary–auxiliary pairs — whose endpoints' blankets were never guaranteed to
be inside the module — are discarded. Where two modules both propose an edge,
`weight.summary` decides: `"min"` (the default) keeps the edge only if every
expected proposer found it, a consensus rule; `"mean"` averages the signed
weights, so an edge some proposers missed survives with a shrunken weight.

### When to use it

Divide and conquer is designed for high-dimensional datasets, such as
high-throughput RNA-seq, and for analyses that require learning many graphs
or relearning a graph many times (bootstrapping). It is based on the assumption
that the true graph structure underlying the data is modular (as is the case
for many biological systems).

It is the wrong tool when `p` is small enough to fit directly. The worked
example below uses 120 features purely so it runs in seconds. If your real
data is that size you should use the whole-data learner. We recommend
modularDAC for 1500+ feature networks.

---

## Installation

The `main` branch holds the user-facing package and is the default branch:

```r
# install.packages("devtools")
devtools::install_github("LukeBerger/modularDAC")
```

The `dev` branch additionally carries the benchmarking code — graph and data
simulators, ground-truth module constructors, and scoring functions such as
`calc_F1()` and `module_contiguity()` — used to evaluate the method:

```r
devtools::install_github("LukeBerger/modularDAC", ref = "dev")
```

Most dependencies are optional and only checked when you call a function that
needs one. For the example below you need `WGCNA` and `SILGGM`, plus
`SummarizedExperiment` to load the example data:

```r
install.packages(c("WGCNA", "SILGGM", "flashClust", "dynamicTreeCut"))

# WGCNA pulls in a few Bioconductor packages
# install.packages("BiocManager")
BiocManager::install(c("GO.db", "impute", "preprocessCore", "SummarizedExperiment"))
```

---

## A worked example

`lfr_example` ships with the package: 120 features by 60 samples, simulated from
a Lancichinetti–Fortunato–Radicchi benchmark graph that resolves into three
communities of 44, 40 and 36 nodes joined by 378 edges. It is a
`SummarizedExperiment`.

```r
library(modularDAC)
library(SummarizedExperiment)

data(lfr_example)
lfr_example
#> class: SummarizedExperiment
#> dim: 120 60
#> metadata(3): source n.true.modules n.true.edges
#> assays(1): expression
#> rownames(120): Gene_001 Gene_002 ... Gene_119 Gene_120
#> rowData names(2): gene true.module
#> colnames(60): Sample_01 Sample_02 ... Sample_59 Sample_60
#> colData names(1): sample
```

Every modularDAC function takes a plain numeric matrix with **features as rows
and samples as columns**, so pull the assay out first:

```r
x <- assay(lfr_example, "expression")
dim(x)
#> [1] 120  60
```

### 1. Divide — detect modules, then grow them

```r
set.seed(1)
mods <- find_WGCNA_mods(x, min.size = 10, max.size = 60)

lengths(mods$final.mods@index.list)
#>  1  2  3
#> 57 31 32
```

`max.size` caps how large a module may get, which is what bounds the size of
each sub-problem. `find_WGCNA_mods()` returns the WGCNA adjacency alongside two
module objects: `initial.mods` (the natural cut) and `final.mods` (after
splitting anything over `max.size`).

The partition does not overlap yet, so grow it:

```r
fuzzy <- eigen_fuzzy_modules(x, mods$final.mods, max.size = 80, ratio = 1.5)

lengths(fuzzy@index.list)   # core + recruited auxiliary features
#> [1] 80 46 48

lengths(fuzzy@core.list)    # the features each module owns — unchanged
#>  1  2  3
#> 57 31 32
```

Each module grew by up to `ratio` times its original size (capped at
`max.size`), recruiting the outside features most correlated with its
eigengene. The core sets still partition all 120 features: every feature is
owned exactly once, which is what makes the stitch in step 3 well defined.

### 2 & 3. Conquer and combine

```r
dac <- divide_and_conquer(x, fuzzy, n.cores = 1)

dac$graph
#> IGRAPH ... UNW- 120 123 --
#> + attr: name (v/c), weight (e/n)
```

That is the whole pipeline. The result is a weighted, undirected `igraph` over
all 120 features:

```r
names(dac)
#> [1] "module.subgraphs" "graph" "weights" "other.outputs"
```

| element | contents |
|---|---|
| `graph` | the stitched network over every feature |
| `module.subgraphs` | the per-module graphs, before stitching |
| `weights` | a 120 × 120 feature-by-feature weight matrix, reconciled the same way as the graph |
| `other.outputs` | whatever the learner returned per module |

Raise `n.cores` to learn the modules in parallel — the point of the method on
real data. Swap the learner with `graph.learning.func`; for the built-in
`learn_*_graph()` functions you also pass the matching argument wrapper and
output parser (see `?divide_and_conquer`).

### Checking the result

`lfr_example` carries the community each feature was simulated from, so you can
see how the detected modules line up with the truth. Real data has no such
column:

```r
table(detected = mods$final.mods@index.vector,
      true     = rowData(lfr_example)$true.module)
#>         true
#> detected  1  2  3
#>        1  9 10 38
#>        2  7 22  2
#>        3 24  4  4
```

Recovery is partial here, and that is expected rather than a bug: with 60
samples and 120 features there are fewer observations than features, so the
correlations the detection step relies on are noisy. Module detection sharpens
considerably as the sample size grows relative to the feature count. It is worth
checking this on your own data before trusting the modules — `module_match()`
and `module_contiguity()` on the `dev` branch exist for exactly that.

---

## What else is in the package

**Module detection** — `find_WGCNA_mods()`, `find_ICA_mods()`

**Module growth** — `eigen_fuzzy_modules()`, `adj_fuzzy_modules()`

**Graph learning** — `learn_SILGGM_graph()`, `learn_WGCNA_graph()`,
`learn_ARACNE_graph()`, `learn_CLR_graph()`, `learn_GENIE3_graph()`,
`learn_bdgraph_graph()`. Each returns `list(graph = <igraph>, weights = <matrix>)`
and can be used on its own, without divide-and-conquer.

**The driver** — `divide_and_conquer()`

The `module` S4 class carries a module set: which features belong to each module
(`index.list` / `name.list`), which it owns (`core.list`), and whether the
modules overlap. See `?"module-class"`.

---

## License

MIT. See [LICENSE](LICENSE).
