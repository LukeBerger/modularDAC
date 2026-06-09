#####################
### DAC Algorithm ###
#####################

#' Learn a network from a data matrix using a divide-and-conquer strategy: data is split into overlapping modules, graphs are learned per module, and the results are stitched back together
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param subgraph.module a module S4 object with overlapping index sets defining which feature subsets to use for each sub-graph
#' @param keep.all.edges a logical, if TRUE all edges from all sub-graphs are retained; if FALSE only edges consistent across overlapping sub-graphs are kept
#' @param n.cores an integer, the number of cores to use for parallel processing
#' @param graph.learning.func a function that accepts a data matrix and returns an igraph object
#' @param arg.wrapping.func a function that packages each module's data subset and additional arguments into a list ready for graph.learning.func
#' @param out.parsing.func a function that extracts a list of igraph objects and any other outputs from the raw output of graph.learning.func
#' @param packages.to.each a character vector of package names to load on each parallel worker
#' @param export.to.each a character vector of function names to export to each parallel worker
#' @param ... additional arguments passed through arg.wrapping.func to graph.learning.func

#' @return a list with elements: modular.subgraphs (list of igraph objects, one per module), final.graph (igraph object combining all sub-graphs), and other.outputs (any additional outputs from graph learning)

#' @importFrom foreach foreach %dopar%
#' @importFrom igraph is_igraph

#' @export
divide_and_conquer <- function(x,
                               subgraph.module,
                               keep.all.edges = FALSE,
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
                                                  ".pvalue",
                                                  ".abs_pcor_filter"),
                               ...
){

  # check inputs

  # check if data is a numeric matrix
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'data' should be a p x n numeric matrix ")
  }

  # extract index list
  subgraph.index.list <- subgraph.module@index.list

  # check 'subgraph.index.list' is a list of integer vectors with overlapping values
  if (!is.list(subgraph.index.list) || length(subgraph.index.list) == 0) {
    stop("'subgraph.index.list' must be a non-empty list of integer vectors.")
  }
  row.range <- c(1:nrow(x)) # get rows of data  matrix
  for (i in seq_along(subgraph.index.list)) { # for each subgraph
    sg <- subgraph.index.list[[i]] # get the subgraph node indices
    other.indexs <- unlist(subgraph.index.list[-i]) # and the other graphs indices
    # check that all values are numeric
    if (!is.numeric(sg) || any(sg %% 1 != 0)) {
      stop(sprintf("Subgraph %d has non integer values.", i))
    }
    # check that all values correspond to rows in data matrix
    if (any(!(sg %in% row.range))) {
      stop(sprintf("Subgraph %d contains invalid row indices.", i))
    }
    # check that the subgraphs all contain overlaps (necessary for reconnecting after learning)
    if (!any(sg %in% other.indexs)){
      warning(sprintf("Subgraph %d has no overlaps with the other subgraphs.", i))
    }
  }

  # check every row index is in at least one subGraph
  if (!all(row.range %in% unique(unlist(subgraph.index.list)))) {
    stop("Not all rows of data are covered by 'subgraph.index.list'.")
  }

  # check if both functions exist
  if (!is.function(arg.wrapping.func)){
    stop("Argument wrapping function is not a function.")
  }
  if (!is.function(graph.learning.func)){
    stop("Graph learning function is not a function.")
  }
  if (!is.function(out.parsing.func)){
    stop("Output Parsing function is not a function.")
  }
  arg.wrapping.func <- match.fun(arg.wrapping.func)
  graph.learning.func <- match.fun(graph.learning.func)
  out.parsing.func <- match.fun(out.parsing.func)

  # set up multithreading

  # if given access to multiple cores...
  if (n.cores > 1) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  # divide data and learn graphs

  # divide data based on subGraph indices
  sub.x  <- lapply(subgraph.index.list, function(sg){
    x[sg,]
  })

  # learn graphs in parallel
  # nOTE: this function can be substituted for any graph learning function that returns and Igraph object
  args <- arg.wrapping.func(sub.x, ...)
  graph.learning.outputs <- foreach::foreach(i = seq_along(args),
                                             .combine = 'list',
                                             .multicombine = TRUE,
                                             .maxcombine = length(args),
                                             .packages = packages.to.each,
                                             .export = export.to.each
  ) %dopar% {
    do.call(graph.learning.func, args[[i]])
  }

  # separate outputs of graph learning into a list of learned igraph objects and other outputs
  parsed.outputs <- out.parsing.func(graph.learning.outputs)

  if(!all(unlist(lapply(parsed.outputs$learned.graphs , igraph::is_igraph)))){
    stop("The 'learned graphs' within 'parsed outputs' are not all igraph objects.
         Check your graph learning and output parsing functions
         and ensure parsed.outputs$learned.graphs contains a list of igraphs with overlapping node names")
  }

  # stitch graphs back together

  final.graph <- .connect_subgraphs(x, parsed.outputs$learned.graphs, keep.all.edges)

  # return sub graphs and final graphs
  return(
    list(
      modular.subgraphs = parsed.outputs$learned.graphs,
      final.graph = final.graph,
      other.outputs = parsed.outputs$other.outputs
    )
  )
}

#' Combine a list of sub-graphs with overlapping nodes into a single graph by reconciling edge disagreements
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param sub.graphs a list of igraph objects with overlapping node names
#' @param keep.all.edges a logical, if TRUE all edges are retained; if FALSE edges that appear in only one of two overlapping sub-graphs are removed

#' @importFrom igraph union as_adjacency_matrix graph_from_adjacency_matrix V as_edgelist delete_edges get_edge_ids
#' @importFrom stringr str_remove
#' @importFrom utils combn
#' @importFrom dplyr bind_rows anti_join %>%

#' @return an igraph object comprised of the combine subgraphs

#' @keywords internal
.connect_subgraphs <- function(x, sub.graphs, keep.all.edges = FALSE){
  # new graph containing all edges
  full.graph <- do.call(igraph::union, sub.graphs)

  # if keeping all edges, return that graph
  if(keep.all.edges){return(full.graph)}

  # rebuild from adj matrix to properly reorder nodes and remove unwanted attributes
  adj <- as.matrix(igraph::as_adjacency_matrix(full.graph))
  original.order <- match(rownames(x), rownames(adj))
  adj <- adj[original.order , original.order]
  full.graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

  # remove non overlap edges
  sg.pairs <- utils::combn(seq_along(sub.graphs), m = 2, simplify = F) # get all pairs of subgraphs
  edges.to.remove <- lapply(sg.pairs, function(pr){ # for each pair of subgraphs
    # get index for 'left' and 'right' subgraphs
    left = pr[1]
    right = pr[2]

    # get set of nodes that are in both graphs
    overlapping.nodes <- intersect(
      igraph::V(sub.graphs[[left]])$name, igraph::V(sub.graphs[[right]])$name
    )

    # if there are overlapping nodes
    if(length(overlapping.nodes) > 1){
      # get edge list for these nodes
      left.edges <- igraph::as_edgelist(sub.graphs[[left]])
      right.edges <- igraph::as_edgelist(sub.graphs[[right]])

      # if both graphs contain edges
      if(nrow(left.edges) > 0 && nrow(right.edges) > 0){
        left.edges <-left.edges[
          which(
            apply(left.edges, 1, function(row){
              all(row %in% overlapping.nodes)
            })
          )
          ,] # take only edges between overlap nodes

        right.edges <-right.edges[
          which(
            apply(right.edges, 1, function(row){
              all(row %in% overlapping.nodes)
            })
          )
          ,] # for both sets of edges

        # this will return a vector in the edge case in which there is only one edge with overlap nodes
        if(!is.matrix(left.edges) && length(left.edges) == 2){
          left.edges <- matrix(left.edges, ncol = 2, byrow = T)
        }
        if(!is.matrix(right.edges) && length(right.edges) == 2){
          right.edges <- matrix(right.edges, ncol = 2, byrow = T)
        }

        # if these matrixes have at least one row
        if(nrow(left.edges) > 0 && nrow(right.edges) > 0){
          # get the edges that dont occur in both sets
          sym.diff <- dplyr::bind_rows(
            as.data.frame(left.edges) %>% dplyr::anti_join(as.data.frame(right.edges), by = c("V1","V2")),
            as.data.frame(right.edges) %>% dplyr::anti_join(as.data.frame(left.edges), by = c("V1","V2"))
          )
          # return the matrix of edges to remove
          return(as.matrix(sym.diff))
        }

      }
    }
    return(matrix(nrow = 0 , ncol = 2))
  })

  # convert list of matrixs into single matrix containing edge information
  edges.to.remove <-as.matrix(unique(do.call(rbind, edges.to.remove)))

  # remove these edges from combined graoh
  full.graph <- igraph::delete_edges(full.graph, igraph::get_edge_ids(full.graph, as.vector(t(edges.to.remove)), directed = FALSE))

  return(full.graph)
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
                                 global = T,
                                 alpha = 0.05,
                                 fdr.filter = T,
                                 max.fdr = 0.05,
                                 pos.cut = 0,
                                 neg.cut = 0){

  # return a list of packaged args
  lapply(sub.x, function(x){
    list(x = t(x), # transform to fit n x p matrix required by silggm
         method  = method,
         global = global,
         alpha = alpha,
         fdr.filter = fdr.filter,
         max.fdr = max.fdr,
         pos.cut = pos.cut,
         neg.cut = pos.cut)
  })
}

#' Default output parser for divide_and_conquer: returns the graph learning outputs as learned graphs with no other outputs

#' @param graph.learning.outputs a list of outputs, one per sub-graph, returned by graph.learning.func

#' @return a list containing the learned subgraphs and any other outputs of graph learning

#' @keywords internal
.default_output_parser <- function(graph.learning.outputs){
  return(
    list(
      learned.graphs = graph.learning.outputs,
      other.outputs = NULL
    )
  )
}


#' Argument wrapper for BDgraph: transposes each module data subset for use with BDgraph::bdgraph
#' @param sub.x a list of numeric matrices, each being a p_i x n subset of x for one module (features as rows)

#' @return a list containing the learned subgraphs data transformed for bdgraph

#' @keywords internal
.bd_arg_wrapper <- function(sub.x){
  lapply(sub.x, function(x){
    list(
      x = t(x)
    )
  })
}

#' Graph learning wrapper for BDgraph: fits a Gaussian graphical model and returns an igraph object
#' @param x a numeric matrix with n samples (rows) and p features (columns), as prepared by .bd_arg_wrapper

#' @return a list containing the learned subgraphs data transformed

#' @importFrom igraph graph_from_adjacency_matrix

#' @keywords internal
.bd_wrapper <- function(x){
  if (!requireNamespace("BDgraph", quietly = TRUE)) {
    stop("Package BDgraph is required. Install with: install.packages('BDgraph')", call. = FALSE)
  }
  bdFit <- BDgraph::bdgraph(x)
  postProb <- BDgraph::plinks(bdFit)
  adjMat <- BDgraph::select(postProb)
  igraph::graph_from_adjacency_matrix(adjMat, mode = "undirected") # can change this if I want to test edge directional update for .connect_subgraphs
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
                               n_cores =NULL,
                               CI = 0.95,
                               filter = "pval",
                               threshold = 0.005){

  # return a list of packaged args,
  lapply(sub.x, function(x){
    list(dat = t(x), # transform to fit n x p matrix required by silggm
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
#' @param dat a numeric matrix with n samples (rows) and p features (columns), as prepared by .RSNet_arg_wrapper
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
                                 n_cores =NULL,
                                 CI = 0.95,
                                 filter = "pval",
                                 threshold = 0.05){
  if (!requireNamespace("RSNet", quietly = TRUE)) {
    stop("Package RSNet is required. Install with: install.packages('RSNet')", call. = FALSE)
  }

  ens.net <- RSNet::capture_all(RSNet::ensemble_ggm(dat,
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

