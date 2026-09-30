# modularDAC

Learn large networks by dividing the data into overlapping modules, learning the modular sub-networks independently, and combining the results  into a single graph.

---

## The divide-and-conquer algorithm

Network inference methods that estimate partial correlation have to reason about all `p * (p - 1) / 2` feature pairs at once. Computational cost and runtime grows steeply with the number of features, and once `p` approaches or exceeds the number of samples `n` inference can be unstable.

modularDAC exploits the fact that biological networks are often **modular**. Features
fall into correlated groups, and primarily have edges within their group. If this assumption holds, the network can be recovered from many small problems instead of one large one.

The algorithm runs in three stages.

**1. Divide.** 
Partition the features into modules, then expand them to create overlapping **fuzzy modules**.

The partition comes from `find_WGCNA_mods()` (co-expression modules via WGCNA) or `find_ICA_mods()` (independent components, which overlap naturally). The expansion step uses `eigen_fuzzy_modules()` or `adj_fuzzy_modules()` to recruit  features outside a module that correlate most strongly with it, up to a size or ratio cap.

The expansion is necessary in order  to  1) ensure modules contain overlapping pairs of nodes which can be used to stitch together the final network and 2) provide nodes at the boundaries of the module with the necessary information to condition their partial correlations on.


**2. Conquer.** Learn a graph within each (fuzzy) module independently, optionally in parallel across cores, using a  inference function. The default is `learn_SILGGM_graph()`.

**3. Combine.** Stitch the sub-networks together based on overlapping edges. By default `weight.summary = "min"`: edges are only kept if they are found in every sub-network in which they are possible. Edge **AB** is only included in the final graph if it is found in every fuzzy module that contains both  **A** and **B** 

Using `weight.summary = "mean"` will instead take the average signed weight of **AB** from every possible sub-network, keeping edges at reduced weight if even if weight = 0 in some networks.

### When to use it

Divide and conquer is designed for high-dimensional datasets  and for analyses that require learning many graphs or relearning a graph many times, such as bootstrapping. It is based on the assumption that the true graph structure underlying the data is modular (as is the case for many biological systems).

It is the wrong tool when `p` is small enough to fit directly. The worked example below uses 120 features purely so it runs in seconds. If your real data is that size you should use the whole-data learner. We recommend modularDAC for 1500+ feature networks.

---

## Installation

The `main` branch holds the user-facing package and is the default branch:

```r
# install.packages("devtools")
devtools::install_github("LukeBerger/modularDAC")
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

Raise `n.cores` to learn the modules in parallel. Swap the learner with `graph.learning.func`; for the built-in
`learn_*_graph()` functions you also pass the matching argument wrapper and
output parser (see `?divide_and_conquer`).

---

## What is in the package

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

## Using a custom inference function 

modularDAC allows users to input any inference method of their  choosing to the `divide_and_conquer` function provided its inputs and outputs are formatted properly. Here we use a custom function, which infers totally random edge weights, as an example.

```r
learn_random_graph <- function(x, edge.prob = 0.15){
  features <- rownames(x)   # x has features as rows and samples as columns
  n.features <- length(features)

  # each feature pair gets a random weight with probability edge.prob, and 0 otherwise
  weights <- matrix(0, n.features, n.features, dimnames = list(features, features))
  is.edge <- upper.tri(weights) & (matrix(runif(n.features^2), n.features) < edge.prob)
  weights[is.edge] <- rnorm(sum(is.edge), mean = 0.1, sd = 0.05)
  weights <- weights + t(weights)   # symmetric, with a zero diagonal

  # a weighted, undirected igraph whose vertices are named after the features
  graph <- igraph::graph_from_adjacency_matrix(weights, mode = "undirected",
                                               weighted = TRUE, diag = FALSE)

  return(list(graph = graph, weights = weights, number_edges = sum(is.edge)))
}
```

The learner is called once per module, on that module's slice of `x` (still
features by samples). The one hard requirement is that it produces an `igraph`
whose vertices are named after the rows of `x`. The sub-graphs are stitched back
together by vertex name, so the features the modules share must have the same
name in each. Anything else it returns, like `number_edges` here, is up to
you.

`divide_and_conquer()` does not call the learner directly. It hands each module's
data to an **argument wrapper**, and hands the learner's outputs to an **output
parser**. You write both to fit your learner.

### The argument wrapper

The wrapper receives `sub.x`, a list with one data matrix per module, plus any
extra arguments given to `divide_and_conquer()`. It returns one argument list
per module, and module `i` is learned with
`do.call(graph.learning.func, args[[i]])`, so the names in each list must match
the learner's arguments:

```r
random_arg_wrapper <- function(sub.x, edge.prob = 0.15){
  lapply(sub.x, function(x){
    list(x = x, edge.prob = edge.prob)
  })
}
```

### The output parser

The parser receives a list of the learner's outputs, one per module, and must
return a list with exactly these two elements:

| element | contents |
|---|---|
| `learned.graphs` | a list of the per-module `igraph`s; these are stitched into `graph` |
| `other.outputs` | a list of anything else to keep, one entry per module; returned as `other.outputs` |

```r
random_output_parser <- function(graph.learning.outputs){
  return(
    list(
      learned.graphs = lapply(graph.learning.outputs, function(o) o$graph),
      other.outputs  = lapply(graph.learning.outputs, function(o){
        list(weights = o$weights, number_edges = o$number_edges)
      })
    )
  )
}
```

The combined `weights` matrix is built from `other.outputs` when each entry is
itself a feature-named weight matrix, which is how the default SILGGM parser
works. Otherwise it is built from the edge weights of the learned graphs. Here
each entry is a list, so the graphs are used. That gives the same result, since
`learn_random_graph()` only puts non-zero weights on edges.

### Running it

Pass all three functions to `divide_and_conquer()`, using `x` and `fuzzy` from
the worked example:

```r
set.seed(1)
random.dac <- divide_and_conquer(
  x, fuzzy,
  graph.learning.func = learn_random_graph,
  arg.wrapping.func   = random_arg_wrapper,
  out.parsing.func    = random_output_parser,
  packages.to.each    = "igraph",
  export.to.each      = "learn_random_graph",
  edge.prob = 0.15
)
```

`edge.prob` is not an argument of `divide_and_conquer()`, so it is passed through to
`random_arg_wrapper()` and from there into every call to `learn_random_graph()`.
`packages.to.each` and `export.to.each` only matter when `n.cores > 1`: they
list the packages each parallel worker loads and the functions it is sent. The
defaults are set up for SILGGM, so replace them whenever you swap the learner.

Each module's `number_edges` comes back in `other.outputs`:

```r
sapply(random.dac$other.outputs, function(o) o$number_edges)
#> [1] 473 145 164

random.dac$graph
#> IGRAPH ... UNW- 120 628 --
#> + attr: name (v/c), weight (e/n)
```

Each module draws its edges independently, so the modules disagree about the
features they share. The stitch settles those disagreements with the same
ownership and `weight.summary` rules as any other learner, and 628 of the 782
edges drawn across the three modules survive.

Because this learner is random, `set.seed()` only makes the result reproducible
when `n.cores = 1`. Parallel workers each draw from their own random stream.

---

## License

MIT. See [LICENSE](LICENSE).
