#' Detect co-expression modules by Louvain community detection on a feature graph
#'
#' Learns a sparse feature graph from the data with ARACNE and partitions it with
#' the Louvain algorithm, then enforces the module size limits. ARACNE is run
#' directly through minet: mutual information from [minet::build.mim()], indirect
#' edges pruned by the data processing inequality (DPI) in [minet::aracne()], and
#' the surviving edges kept only where their mutual information is significant
#' (Benjamini-Hochberg adjusted p-value below `max.fdr`), the same edge rule and
#' independence test as [learn_ARACNE_graph()]. `percentile.threshold` replaces
#' that test with a cut of the weakest DPI survivors.
#'
#' The DPI is kept on purpose even though modules do not need indirect edges
#' removed: in data with a strong shared component, the raw significant
#' mutual-information graph links most pairs across modules, and it is the DPI
#' that strips that cross-module association out.
#'
#' A community larger than `max.size` is
#' re-split by Louvain on its own subgraph, falling back to average-linkage
#' clustering on `1 - |cor|` if that subgraph will not split. A community smaller
#' than `min.size` is merged, smallest first, into the community it shares the
#' most edge weight with (the highest mean absolute correlation if it has no
#' edges to any), provided the merged community stays within `max.size`.
#'
#' Unlike [find_WGCNA_mods()], there is no unassigned group: a feature with no
#' edges starts as its own community and is merged like any other undersized
#' one, so every feature lands in a real module and the cut cannot come back
#' empty or as one module spanning the whole graph.
#'
#' The returned modules do not overlap. To build overlapping modules for
#' [divide_and_conquer()], grow them along the same graph with
#' `adj_fuzzy_modules(x, res$graph.adj, res$final.mods, max.size, ratio, min.adj = 1e-12)`:
#' because `graph.adj` is zero between features with no edge, any positive
#' `min.adj` recruits only graph neighbours of a module's core nodes, an
#' approximation to their Markov blankets.
#'
#' Louvain is randomised, so call [set.seed()] first for a reproducible
#' partition.
#'
#' @param x a numeric matrix with p features (rows) and n samples (columns), with unique row names
#' @param graph an igraph object or NULL (the default). If NULL, an ARACNE graph is learned from `x` (see `eps`, `max.fdr` and `percentile.threshold`). Otherwise an undirected graph whose vertex names are feature names; vertices not in `rownames(x)` are dropped, features of `x` absent from the graph start as singleton communities, and edge weights (if any) are used as absolute values, with an unweighted graph treated as all weights 1
#' @param min.size an integer, the minimum number of nodes in a module; smaller communities are merged into a neighbour where the merged size stays within `max.size`, and any that cannot be merged are kept, with a warning
#' @param max.size an integer, or NULL for no cap; the maximum number of nodes allowed in a module
#' @param resolution a numeric, the Louvain resolution passed to [igraph::cluster_louvain()] for the initial partition and for re-splitting communities larger than `max.size`; higher values give more, smaller communities. The default 0.5 sits below igraph's 1 because the DPI-pruned ARACNE graph has almost no triangles, and at resolution 1 Louvain fragments such a near-tree into many small communities that the size rules then have to merge back
#' @param eps a numeric, the DPI tolerance passed to [minet::aracne()]; used only when `graph` is NULL
#' @param max.fdr a numeric between 0 and 1, used only when `graph` is NULL and `percentile.threshold` is NULL: the default edge rule (0.05). A DPI survivor is kept when the Benjamini-Hochberg adjusted p-value of its mutual information, tested against the independence null over all node pairs, is below `max.fdr`; as in [learn_ARACNE_graph()], the test is run on the raw mutual information, before DPI
#' @param percentile.threshold a numeric between 0 and 1, or NULL (the default, leaving this rule off), used only when `graph` is NULL. Supplying it replaces the `max.fdr` test: the ARACNE edges that survive DPI are ranked by mutual information and the weakest `percentile.threshold` fraction of them is dropped, so 0 keeps every survivor and leaves DPI as the only pruning. Unlike `max.fdr`, a positive value is not calibrated across graph sizes. Ranking within the survivors (not all node pairs) keeps the cut responsive at any sparsity

#' @return a list containing: 'graph', the igraph the communities were detected on (restricted to `rownames(x)`, with the absolute edge weights Louvain used); 'graph.adj', the p x p matrix of those weights in `rownames(x)` order (0 where there is no edge), ready to pass as `adj` to [adj_fuzzy_modules()]; 'initial.mods', the raw Louvain partition; and 'final.mods', that partition after splitting and merging to satisfy the size limits

#' @importFrom igraph cluster_louvain membership induced_subgraph V E as_edgelist edge_attr_names ecount vcount is_directed graph_from_adjacency_matrix
#' @importFrom stats cor hclust as.dist cutree setNames
#' @importFrom methods new

#' @export
find_Louvain_mods <- function(x,
                              graph = NULL,
                              min.size = 10,
                              max.size = NULL,
                              resolution = 0.5,
                              eps = 0.001,
                              max.fdr = 0.05,
                              percentile.threshold = NULL) {
  if (is.null(rownames(x)) || anyDuplicated(rownames(x))) {
    stop("x must have unique row names; they identify the features as graph vertices.", call. = FALSE)
  }

  # treat an unconstrained max.size as Inf (so the size checks below are explicit)
  if (is.null(max.size)) max.size <- Inf
  if (min.size > max.size) {
    stop("min.size (", min.size, ") cannot exceed max.size (", max.size, ").", call. = FALSE)
  }

  # drop features with no variance (they cannot be placed in any module)
  rm.r <- which(rowSums((x - rowMeans(x))^2) == 0)
  if (length(rm.r) > 0) {
    warning("Some features have zero variance and were removed; they could not be placed in any module.")
    message("Removed rows: ", paste(rm.r, collapse = ", "))
    x <- x[-rm.r, , drop = FALSE]
  }
  rn <- rownames(x)
  p  <- nrow(x)

  # learn the ARACNE graph, or check the one supplied
  if (is.null(graph)) {
    if (!requireNamespace("minet", quietly = TRUE)) {
      stop("Package minet is required to learn the ARACNE graph. Install with: ",
           "BiocManager::install('minet'), or supply a graph.", call. = FALSE)
    }
    # percentile.threshold, when supplied, replaces the default max.fdr rule
    use.percentile <- !is.null(percentile.threshold)
    if (use.percentile) {
      if (!is.numeric(percentile.threshold) || length(percentile.threshold) != 1 ||
          percentile.threshold < 0 || percentile.threshold > 1) {
        stop("percentile.threshold must be a single number between 0 and 1.", call. = FALSE)
      }
    } else if (!is.numeric(max.fdr) || length(max.fdr) != 1 || max.fdr <= 0 || max.fdr > 1) {
      stop("max.fdr must be a single number in (0, 1], or supply percentile.threshold.",
           call. = FALSE)
    }
    message("Learning the ARACNE graph...")
    # mutual information (minet expects samples x features), then DPI pruning
    mim <- minet::build.mim(dataset = t(x))
    aracne.mat <- minet::aracne(mim, eps = eps)
    diag(aracne.mat) <- 0
    if (use.percentile) {
      # drop the weakest percentile.threshold of the DPI survivors; >= so that 0
      # keeps every survivor
      cut <- .percentile_threshold(aracne.mat, percentile.threshold,
                                   "ARACNE mutual information", nonzero.only = TRUE)
      aracne.mat[aracne.mat < cut] <- 0
    } else {
      # keep the survivors whose raw MI is significant under the independence
      # null (BH over all pairs); tested on mim, not the DPI output, exactly as
      # learn_ARACNE_graph() does
      aracne.mat[.matrix_p_adjust(.mi_pvalue(mim, ncol(x))) >= max.fdr] <- 0
    }
    graph <- igraph::graph_from_adjacency_matrix(aracne.mat, mode = "undirected",
                                                 weighted = TRUE)
  } else {
    if (!inherits(graph, "igraph")) stop("graph must be an igraph object.", call. = FALSE)
    if (igraph::is_directed(graph)) stop("graph must be undirected.", call. = FALSE)
    if (is.null(igraph::V(graph)$name)) {
      stop("graph's vertices must be named by feature (rownames(x)).", call. = FALSE)
    }
  }

  # restrict the graph to the features in x
  graph <- igraph::induced_subgraph(graph, intersect(igraph::V(graph)$name, rn))
  if (igraph::vcount(graph) == 0) {
    stop("None of graph's vertex names match rownames(x).", call. = FALSE)
  }

  # Louvain needs positive weights: use absolute values (1 when unweighted), and
  # give an edge whose weight is 0 or missing a tiny positive one, so it still
  # counts as an edge here and when growing modules with adj_fuzzy_modules()
  w <- if ("weight" %in% igraph::edge_attr_names(graph)) abs(igraph::E(graph)$weight)
       else rep(1, igraph::ecount(graph))
  w[!is.finite(w) | w <= 0] <- 1e-12
  igraph::E(graph)$weight <- w

  # the same weights as a p x p matrix in x's row order, 0 where there is no edge
  graph.adj <- matrix(0, p, p, dimnames = list(rn, rn))
  el <- igraph::as_edgelist(graph, names = TRUE)
  graph.adj[el] <- w
  graph.adj[el[, 2:1, drop = FALSE]] <- w

  # Louvain on the whole graph. Membership is matched to x's rows by name, since
  # the graph's vertex order need not be rownames(x)'s, and features absent from
  # the graph each get a singleton community of their own.
  message("Detecting Louvain communities...")
  lv <- igraph::membership(igraph::cluster_louvain(graph, weights = w, resolution = resolution))
  memb <- stats::setNames(rep(NA_integer_, p), rn)
  memb[igraph::V(graph)$name] <- as.integer(lv)
  miss <- which(is.na(memb))
  memb[miss] <- max(c(0L, memb), na.rm = TRUE) + seq_along(miss)
  initial.index.vector <- as.numeric(match(memb, sort(unique(memb))))

  # Enforce max.size: split one oversized community at a time by Louvain on its
  # own subgraph, falling back to average-linkage on 1 - |cor| (cut into the
  # fewest pieces that fit) when the subgraph does not split. Splitting runs
  # before merging, so every oversized community holds only graph vertices.
  # A split piece can itself exceed max.size, so the sizes are re-checked after
  # every split.
  sz <- table(memb)
  while (any(sz > max.size)) {
    big   <- as.integer(names(sz)[sz > max.size][1])
    nodes <- rn[memb == big]
    sg  <- igraph::induced_subgraph(graph, nodes)
    # membership() follows the subgraph's vertex order; index it by name
    sub <- stats::setNames(
      as.integer(igraph::membership(igraph::cluster_louvain(sg, weights = igraph::E(sg)$weight,
                                                            resolution = resolution))),
      igraph::V(sg)$name)[nodes]
    if (anyNA(sub) || length(unique(sub)) < 2) {
      d  <- 1 - abs(stats::cor(t(x[nodes, , drop = FALSE])))
      d[is.na(d)] <- 1                                # a constant feature is unrelated to all
      hc <- stats::hclust(stats::as.dist(d), method = "average")
      # the fewest pieces whose largest fits (k = length(nodes) always does)
      k   <- ceiling(length(nodes) / max.size)
      sub <- stats::cutree(hc, k = k)
      while (max(table(sub)) > max.size && k < length(nodes)) {
        k   <- k + 1L
        sub <- stats::cutree(hc, k = k)
      }
    }
    message("Splitting a community of ", length(nodes), " nodes to satisfy max.size.")
    # the first piece keeps the id, the others get fresh ids beyond the current max
    memb[nodes] <- ifelse(sub == 1L, big, max(memb) + sub - 1L)
    sz <- table(memb)
  }

  # Enforce min.size: merge undersized communities, smallest first, into the
  # community they share the most edge weight with (the highest mean |cor| if
  # they have no edge to any), as long as the merged size stays within max.size.
  # A community that fits nowhere is kept: small modules are legal, only noisier.
  # Each merge changes the sizes, so the undersized set is recomputed after every
  # pass; a community that fits nowhere joins `stuck` and is not tried again.
  stuck <- integer(0)
  sz    <- table(memb)
  small <- setdiff(as.integer(names(sz)[sz < min.size]), stuck)
  while (length(small) > 0 && length(sz) > 1L) {
    m      <- small[which.min(sz[as.character(small)])]
    nodes  <- rn[memb == m]
    others <- setdiff(as.integer(names(sz)), m)
    # total edge weight from m to every other community, in one pass over its rows
    score <- tapply(colSums(graph.adj[nodes, , drop = FALSE]), memb, sum)[as.character(others)]
    if (all(score <= 0)) {
      # mean |cor| block = mean of the column means, since every column has |nodes| rows
      r <- abs(stats::cor(t(x[nodes, , drop = FALSE]), t(x)))
      r[is.na(r)] <- 0
      score <- tapply(colMeans(r), memb, mean)[as.character(others)]
    }
    score[(sz[as.character(others)] + length(nodes)) > max.size] <- -Inf
    if (all(!is.finite(score))) {
      stuck <- c(stuck, m)
    } else {
      memb[nodes] <- others[which.max(score)]
    }
    sz    <- table(memb)
    small <- setdiff(as.integer(names(sz)[sz < min.size]), stuck)
  }
  if (length(stuck) > 0) {
    warning(length(stuck), " module(s) remain below min.size because merging them would exceed max.size.",
            call. = FALSE)
  }

  # renumber the final modules to a contiguous 1:N sequence
  final.index.vector <- as.numeric(match(memb, sort(unique(memb))))

  # build the module objects (a non-overlapping partition owns all of its own
  # nodes, so core.list == index.list)
  initial.index.list <- split(seq_len(p), initial.index.vector)
  initial.mods <- methods::new("module",
                               source = "find_Louvain_mods initial mods",
                               data.dim = dim(x),
                               overlapping = FALSE,
                               index.vector = initial.index.vector,
                               index.list = initial.index.list,
                               name.list = split(rn, initial.index.vector),
                               core.list = initial.index.list)

  final.index.list <- split(seq_len(p), final.index.vector)
  final.mods <- methods::new("module",
                             source = "find_Louvain_mods sized mods",
                             data.dim = dim(x),
                             overlapping = FALSE,
                             index.vector = final.index.vector,
                             index.list = final.index.list,
                             name.list = split(rn, final.index.vector),
                             core.list = final.index.list)

  return(list(
    graph = graph,
    graph.adj = graph.adj,
    initial.mods = initial.mods,
    final.mods = final.mods
  ))
}
