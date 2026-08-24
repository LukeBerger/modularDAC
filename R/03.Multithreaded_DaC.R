#####################
### DAC Algorithm ###
#####################

#' Learn a network from a data matrix using a divide-and-conquer strategy: data is split into overlapping modules, graphs are learned per module, and the results are stitched back together
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param subgraph.module a module S4 object with overlapping index sets defining which feature subsets to use for each sub-graph
#' @param weight.summary a character, how to reconcile the weight of an edge that is possible in multiple sub-graphs; 'min' (default) keeps the smallest-magnitude weight across the sub-graphs where the edge is possible (an absent edge counts as weight 0, so an edge missing from any possible sub-graph is dropped), while 'mean' keeps the average of the signed weights across those sub-graphs (so an edge present in only some possible sub-graphs survives with a reduced weight)
#' @param n.cores an integer, the number of cores to use for parallel processing
#' @param graph.learning.func a function that accepts a data matrix and returns an igraph object
#' @param arg.wrapping.func a function that packages each module's data subset and additional arguments into a list ready for graph.learning.func
#' @param out.parsing.func a function that extracts a list of igraph objects and any other outputs from the raw output of graph.learning.func
#' @param packages.to.each a character vector of package names to load on each parallel worker
#' @param export.to.each a character vector of function names to export to each parallel worker
#' @param ... additional arguments passed through arg.wrapping.func to graph.learning.func

#' @return a list with elements: modular.subgraphs (list of igraph objects, one per module), final.graph (igraph object combining all sub-graphs), and other.outputs (any additional outputs from graph learning)

#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom igraph is_igraph

#' @export
divide_and_conquer <- function(x,
                               subgraph.module,
                               weight.summary = c("min", "mean"),
                               n.cores = 1,
                               graph.learning.func = learn_SILGGM_graph,
                               arg.wrapping.func = .default_arg_wrapper,
                               out.parsing.func = .default_output_parser,
                               packages.to.each = c("SILGGM", "igraph"),
                               export.to.each = c("learn_SILGGM_graph",
                                                  ".upper_tri_to_matrix",
                                                  ".matrix_p_adjust",
                                                  ".upper_tri_vec",
                                                  ".pcor_zscore",
                                                  ".pvalue"),
                               ...
){

  # check inputs

  # check if data is a numeric matrix
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' should be a p x n numeric matrix")
  }

  # extract index list
  subgraph.index.list <- subgraph.module@index.list

  # check 'subgraph.index.list' is a list of integer vectors with overlapping values
  if (!is.list(subgraph.index.list) || length(subgraph.index.list) == 0) {
    stop("'subgraph.index.list' must be a non-empty list of integer vectors.")
  }
  row.range <- seq_len(nrow(x)) # valid row indices of the data matrix
  # indices shared by more than one subgraph (used for the overlap check below)
  all.idx <- unlist(subgraph.index.list)
  dup.idx <- unique(all.idx[duplicated(all.idx)])
  for (i in seq_along(subgraph.index.list)) { # for each subgraph
    sg <- subgraph.index.list[[i]] # get the subgraph node indices
    # check that all values are numeric
    if (!is.numeric(sg) || any(sg %% 1 != 0)) {
      stop(sprintf("Subgraph %d has non integer values.", i))
    }
    # check that all values correspond to rows in data matrix
    if (any(!(sg %in% row.range))) {
      stop(sprintf("Subgraph %d contains invalid row indices.", i))
    }
    # check that the subgraphs all contain overlaps (necessary for reconnecting after learning)
    if (!any(sg %in% dup.idx)){
      warning(sprintf("Subgraph %d has no overlaps with the other subgraphs.", i))
    }
  }

  # check every row index is in at least one subGraph
  if (!all(row.range %in% unique(all.idx))) {
    stop("Not all rows of data are covered by 'subgraph.index.list'.")
  }

  # resolve the supplied functions (match.fun accepts a function or a name)
  arg.wrapping.func <- match.fun(arg.wrapping.func)
  graph.learning.func <- match.fun(graph.learning.func)
  out.parsing.func <- match.fun(out.parsing.func)

  # set up multithreading

  # set up the parallel backend; always release it on exit -- stop the workers AND
  # reset foreach to the sequential backend -- so later parallel code starts clean
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doParallel::registerDoParallel(cl)
  } else {
    # run sequentially without leaving (or relying on) a stale backend
    foreach::registerDoSEQ()
  }

  # divide data and learn graphs

  # divide data based on subGraph indices
  sub.x  <- lapply(subgraph.index.list, function(sg){
    x[sg, , drop = FALSE]
  })

  # learn graphs in parallel
  # nOTE: this function can be substituted for any graph learning function that returns and Igraph object
  # (foreach with no .combine returns a list of results, one per module)
  args <- arg.wrapping.func(sub.x, ...)
  graph.learning.outputs <- foreach::foreach(i = seq_along(args),
                                             .packages = packages.to.each,
                                             .export = export.to.each
  ) %dopar% {
    do.call(graph.learning.func, args[[i]])
  }

  # separate outputs of graph learning into a list of learned igraph objects and other outputs
  parsed.outputs <- out.parsing.func(graph.learning.outputs)

  if(!all(vapply(parsed.outputs$learned.graphs, igraph::is_igraph, logical(1)))){
    stop("The 'learned graphs' within 'parsed outputs' are not all igraph objects.
         Check your graph learning and output parsing functions
         and ensure parsed.outputs$learned.graphs contains a list of igraphs with overlapping node names")
  }

  # stitch graphs back together

  weight.summary <- match.arg(weight.summary)
  # node ownership (core.list) makes the stitch owner-based: each edge is credited
  # only from the module(s) that own one of its endpoints. Falls back to NULL
  # (ownership-agnostic) when the module carries no core.list.
  core.sets <- if (length(subgraph.module@core.list) > 0) {
    lapply(subgraph.module@core.list, function(idx) rownames(x)[idx])
  } else {
    NULL
  }
  final.graph <- .connect_subgraphs(x, parsed.outputs$learned.graphs, weight.summary, core.sets)

  # return sub graphs and final graphs
  return(
    list(
      module.subgraphs = parsed.outputs$learned.graphs,
      graph = final.graph,
      other.outputs = parsed.outputs$other.outputs
    )
  )
}

#' Combine a list of sub-graphs with overlapping nodes into a single weighted graph by reconciling edge weights
#'
#' Stitching is done on an edge list rather than dense p x p matrices, so cost
#' scales with the number of within-module edges (which is sparse) instead of
#' with p^2. When \code{core.sets} records node ownership, an edge is only
#' credited from a sub-graph in which at least one endpoint is a \emph{core}
#' (owned) node -- discarding auxiliary-auxiliary edges, whose endpoints' Markov
#' blankets are not guaranteed to lie inside the module. Every kept detection
#' therefore comes from an owner of one of its endpoints, so each edge has at
#' most two authoritative proposers (the owners of its two endpoints).
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param sub.graphs a list of weighted igraph objects with overlapping node names
#' @param weight.summary a character, how to reconcile an edge across the owner modules that could propose it; 'min' (default) keeps the edge only if every expected proposer detected it (consensus), storing the smallest-magnitude signed weight, while 'mean' averages the signed weights over the expected proposers so an edge missing from some proposer survives with a shrunk weight
#' @param core.sets a list (aligned with sub.graphs) of character vectors giving the core (owned) node names of each sub-graph, or NULL to treat every node as core (which reproduces the ownership-agnostic behaviour: all sub-graphs containing both endpoints are proposers)

#' @importFrom igraph as_edgelist graph_from_data_frame V E is_weighted ecount

#' @return a weighted igraph object comprised of the combined sub-graphs

#' @keywords internal
.connect_subgraphs <- function(x, sub.graphs, weight.summary = c("min", "mean"), core.sets = NULL){
  weight.summary <- match.arg(weight.summary)

  # the full node set, in the order of the input data
  nodes <- rownames(x)

  # per-sub-graph node sets; default ownership is "every node is core", which
  # makes every sub-graph containing both endpoints a proposer (legacy behaviour)
  node.sets <- lapply(sub.graphs, function(g) igraph::V(g)$name)
  if (is.null(core.sets)) core.sets <- node.sets

  # an empty graph over the full node set, used as the fallback return
  empty.graph <- igraph::graph_from_data_frame(
    data.frame(from = character(0), to = character(0)),
    directed = FALSE, vertices = data.frame(name = nodes))

  # --- collect authoritative edge detections (>= 1 core endpoint per sub-graph) ---
  rows <- vector("list", length(sub.graphs))
  for (s in seq_along(sub.graphs)) {
    g <- sub.graphs[[s]]
    if (igraph::ecount(g) == 0) next
    el <- igraph::as_edgelist(g, names = TRUE)
    w  <- if (igraph::is_weighted(g)) igraph::E(g)$weight else rep(1, nrow(el))
    core <- core.sets[[s]]
    keep <- (el[, 1] %in% core) | (el[, 2] %in% core)
    if (!any(keep)) next
    el <- el[keep, , drop = FALSE]; w <- w[keep]
    a <- pmin(el[, 1], el[, 2]); b <- pmax(el[, 1], el[, 2])   # canonical endpoint order
    rows[[s]] <- data.frame(a = a, b = b, w = w, sub = s, stringsAsFactors = FALSE)
  }
  edges <- do.call(rbind, rows)
  if (is.null(edges) || nrow(edges) == 0) return(empty.graph)

  # index endpoints into the node universe and key each undirected pair
  ai <- match(edges$a, nodes); bi <- match(edges$b, nodes)
  edges$key <- paste(ai, bi, sep = "-")

  # per-sub boolean presence / core membership over the node universe
  in.sub  <- lapply(node.sets, function(ns) { v <- logical(length(nodes)); v[match(ns, nodes)] <- TRUE; v })
  in.core <- lapply(core.sets, function(cs) { v <- logical(length(nodes)); v[match(cs, nodes)] <- TRUE; v })

  # unique candidate edges
  uk <- !duplicated(edges$key)
  ua <- ai[uk]; ub <- bi[uk]; keys <- edges$key[uk]
  U <- length(keys)

  # expected proposer count E: sub-graphs where both endpoints are present and at
  # least one is core there. Computed from membership (not detections) so a missing
  # edge in a proposer can be detected as such below.
  E <- integer(U)
  for (s in seq_along(sub.graphs)) {
    E <- E + ((in.sub[[s]][ua] & in.sub[[s]][ub]) & (in.core[[s]][ua] | in.core[[s]][ub]))
  }

  # detections grouped by edge (in the order of `keys`)
  det.sub <- split(edges$sub, edges$key)[keys]
  det.w   <- split(edges$w,   edges$key)[keys]

  # reconcile each edge over its expected proposers
  D     <- vapply(det.sub, function(s) length(unique(s)), integer(1))   # proposers that detected it
  n.pos <- vapply(det.w, function(w) sum(w > 0), integer(1))
  n.neg <- vapply(det.w, function(w) sum(w < 0), integer(1))
  if (weight.summary == "min") {
    # consensus: keep only when every expected proposer detected the edge
    min.signed <- vapply(det.w, function(w) w[which.min(abs(w))], numeric(1))
    final.w <- ifelse(D >= E, min.signed, 0)
  } else {
    # average signed weight over the expected proposers (absent-in-some -> shrunk)
    sum.signed <- vapply(det.w, sum, numeric(1))
    final.w <- sum.signed / E
  }

  keep <- final.w != 0

  # warn if any retained edge had conflicting signs across its proposers
  n.conflict <- sum((n.pos > 0) & (n.neg > 0) & keep)
  if (n.conflict > 0) {
    warning(sprintf(
      "%d edge(s) in the final graph had different signs across the sub-graphs in which they appeared.",
      n.conflict
    ))
  }

  if (!any(keep)) return(empty.graph)

  # build and return the combined weighted graph on the full node set
  igraph::graph_from_data_frame(
    data.frame(from = nodes[ua[keep]], to = nodes[ub[keep]], weight = final.w[keep]),
    directed = FALSE, vertices = data.frame(name = nodes))
}

#####################################################
### Wrappers for different graph learning methods ###
#####################################################

#' Default argument wrapper for divide_and_conquer: packages each module data subset with SILGGM arguments

#' @param sub.x a list of numeric matrices, each being a p_i x n subset of x for one module (features as rows)
#' @param method a character, the statistical inference method; one of 'B_NW_SL', 'D-S_NW_SL', 'D-S_GL', 'GFC_SL', or 'GFC_L'
#' @param global a logical, if TRUE global inference is performed across all feature pairs
#' @param alpha a numeric or numeric vector, the pre-specified FDR significance level(s) for edge inclusion (default 0.05)
#' @param fdr.filter a logical, if TRUE the adjacency matrix is filtered to retain only edges below max.fdr
#' @param max.fdr a numeric, the FDR q-value threshold above which an edge is set to zero
#' @param pos.cut a numeric, the minimum absolute partial correlation for a positive edge to be retained
#' @param neg.cut a numeric, the minimum absolute partial correlation for a negative edge to be retained

#' @return a list of arguments to feed to graph learning function

#' @keywords internal
.default_arg_wrapper <- function(sub.x,
                                 method  = "B_NW_SL",
                                 global = TRUE,
                                 alpha = 0.05,
                                 fdr.filter = TRUE,
                                 max.fdr = 0.05,
                                 pos.cut = 0,
                                 neg.cut = 0){

  # return a list of packaged args (each module stays feature x sample;
  # learn_SILGGM_graph transposes internally)
  lapply(sub.x, function(x){
    list(x = x,
         method  = method,
         global = global,
         alpha = alpha,
         fdr.filter = fdr.filter,
         max.fdr = max.fdr,
         pos.cut = pos.cut,
         neg.cut = neg.cut)
  })
}

#' Default output parser for divide_and_conquer: splits each learner's output into the learned graph and the weight matrix

#' @param graph.learning.outputs a list of outputs, one per sub-graph, returned by graph.learning.func; each output is a list with a 'graph' (igraph) and a 'weights' (matrix) element

#' @return a list with 'learned.graphs', a list of the learned igraph objects, and 'other.outputs', a list of the corresponding weight matrices

#' @keywords internal
.default_output_parser <- function(graph.learning.outputs){
  return(
    list(
      learned.graphs = lapply(graph.learning.outputs, function(o) o$graph),
      other.outputs  = lapply(graph.learning.outputs, function(o) o$weights)
    )
  )
}



#' Argument wrapper for RSNet: packages each module data subset with RSNet arguments
#' @param sub.x a list of numeric matrices, each being a p_i x n subset of x for one module (features as rows)
#' @param num_iteration an integer, the number of bootstrap/subsampling iterations to perform
#' @param boot a logical, if TRUE bootstrap resampling is used
#' @param sub_ratio a numeric between 0 and 1, the subsampling ratio applied to samples
#' @param sample_class a character vector or NULL, class labels for samples to include; NULL uses all samples
#' @param correlated a logical, if TRUE correlated feature structure is accounted for
#' @param cluster_ratio a numeric, the clustering ratio parameter
#' @param estimate_CI a logical, if TRUE confidence intervals are estimated for each edge
#' @param method a character, the statistical inference method passed to RSNet (e.g. 'B_NW_SL')
#' @param n_cores an integer or NULL, the number of cores for parallel processing; NULL uses a single core
#' @param CI a numeric between 0 and 1, the confidence level for interval estimation
#' @param filter a character, the filtering criterion applied to edges (e.g. 'pval')
#' @param threshold a numeric, the threshold value for the chosen filter

#' @return a list containing the args for rsnet learning

#' @keywords internal
.RSNet_arg_wrapper <- function(sub.x,
                               num_iteration  = 1,
                               boot = FALSE,
                               sub_ratio = 1,
                               sample_class = NULL,
                               correlated = FALSE,
                               cluster_ratio = 1,
                               estimate_CI = TRUE,
                               method = "B_NW_SL",
                               n_cores = NULL,
                               CI = 0.95,
                               filter = "pval",
                               threshold = 0.005){

  # return a list of packaged args (each module stays feature x sample;
  # .RSNet_full_learning transposes internally)
  lapply(sub.x, function(x){
    list(dat = x,
         num_iteration = num_iteration,
         boot = boot,
         sub_ratio = sub_ratio,
         sample_class = sample_class,
         correlated = correlated,
         cluster_ratio = cluster_ratio,
         estimate_CI = estimate_CI,
         method = method,
         n_cores = n_cores,
         CI = CI,
         filter = filter,
         threshold = threshold)
  })
}

#' Graph learning wrapper for RSNet: generates ensemble and consensus GGM networks from a module data subset
#' @param dat a numeric matrix with p features (rows) and n samples (columns), as prepared by .RSNet_arg_wrapper
#' @param num_iteration an integer, the number of bootstrap/subsampling iterations to perform
#' @param boot a logical, if TRUE bootstrap resampling is used
#' @param sub_ratio a numeric between 0 and 1, the subsampling ratio applied to samples
#' @param sample_class a character vector or NULL, class labels for samples to include; NULL uses all samples
#' @param correlated a logical, if TRUE correlated feature structure is accounted for
#' @param cluster_ratio a numeric, the clustering ratio parameter
#' @param estimate_CI a logical, if TRUE confidence intervals are estimated for each edge
#' @param method a character, the statistical inference method passed to RSNet (e.g. 'B_NW_SL')
#' @param n_cores an integer or NULL, the number of cores for parallel processing; NULL uses a single core
#' @param CI a numeric between 0 and 1, the confidence level for interval estimation
#' @param filter a character, the filtering criterion applied to edges (e.g. 'pval')
#' @param threshold a numeric, the threshold value for the chosen filter


#' @return a list of learned graphs (and other consensus_net_ggm outputs)

#' @keywords internal
.RSNet_full_learning <- function(dat,
                                 num_iteration = 1,
                                 boot = FALSE,
                                 sub_ratio = 1,
                                 sample_class = NULL,
                                 correlated = FALSE,
                                 cluster_ratio = 1,
                                 estimate_CI = TRUE,
                                 method = "B_NW_SL",
                                 n_cores = NULL,
                                 CI = 0.95,
                                 filter = "pval",
                                 threshold = 0.05){
  if (!requireNamespace("RSNet", quietly = TRUE)) {
    stop("Package RSNet is required. Install with: install.packages('RSNet')", call. = FALSE)
  }

  # RSNet expects an n samples x p features matrix
  t.dat <- t(dat)

  ens.net <- RSNet::capture_all(RSNet::ensemble_ggm(t.dat,
                                 num_iteration,
                                 boot,
                                 sub_ratio,
                                 sample_class,
                                 correlated,
                                 cluster_ratio,
                                 estimate_CI,
                                 method,
                                 n_cores = n_cores
  ))
  conc.nets <-  RSNet::capture_all(RSNet::consensus_net_ggm(ggm_networks = ens.net,
                                        CI = CI,
                                        filter = filter,
                                        threshold = threshold

  ))
  return(conc.nets)
}

#' Output parser for RSNet: extracts consensus networks and other outputs from RSNet graph learning results
#' @param graph.learning.outputs a list of outputs, one per sub-graph, returned by .RSNet_full_learning

#' @return a list containing the learned subgraphs and any other outputs of graph learning

#' @keywords internal
.RSNet_output_parser <- function(graph.learning.outputs){
  nets <- lapply(graph.learning.outputs, function(o) o$consensus_network)
  other <- lapply(graph.learning.outputs, function(o) o[-1])
  return(
    list(
      learned.graphs = nets,
      other.outputs = other
    )
  )
}

