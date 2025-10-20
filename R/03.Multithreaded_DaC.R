# library(igraph)
# library(stringr)
# library(dplyr)
# library(parallel)
# library(foreach)
# library(doParallel)


### Takes a network and list of vectors defining
### subgraphs/modules (overlapping sets of of Nodes)
### Learns the sub graphs as individual graphs
### Reconnects them based on overlapping edges
### And Returns a full graph

## Args:
## data, an p (variable) x n (sample) matrix containing the data to build the network
## subGraphs, the subgraphs/modules defined by vectors of ints within a list

fastDivideAndConquer <- function(x, # input data, p x n matrix
                                 modules, # overlappping subsets of p
                                 keep.all.edges = FALSE, # whether to keep all edges in the overlap graph
                                 graph.learning.func = learn_SILGGM_graph,  # function to learn graphs
                                 arg.wrapping.func = .default_arg_wrapper,  # function to wrap each module of x with args for graph learning function
                                 out.parsing.func = .default_output_parser, # function to parse outputs of graph learning into list of igraphs and other outputs
                                 packages.to.each = c("SILGGM", "igraph"),  # packages for each thread
                                 export.to.each = c("learn_SILGGM_graph",
                                                    "upper_tri_to_matrix",
                                                    "matrix_p_adjust",
                                                    ".upper_tri_vec",
                                                    ".pcor_zscore",
                                                    ".pvalue",
                                                    ".abs_pcor_filter"),  # functions to export to each thread
                                 ... # other args to pass through arg wrapper to graph learning
){
  ##################
  ###Check Inputs###
  ##################

  #Check if data is a numeric matrix
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'data' should be a p x n numeric matrix ")
  }

  #Check 'modules' is a list of integer vectors with overlapping values
  if (!is.list(modules) || length(modules) == 0) {
    stop("'modules' must be a non-empty list of integer vectors.")
  }
  row.range <- c(1:nrow(x)) #get rows of data  matrix
  for (i in seq_along(modules)) { #for each subgraph
    sg <- modules[[i]] #get the subgraph node indices
    other.indexs <- unlist(modules[-i]) #and the other graphs indices
    #check that all values are numeric
    if (!is.numeric(sg) || any(sg %% 1 != 0)) {
      stop(sprintf("Subgraph %d has non integer values.", i))
    }
    #check that all values correspond to rows in data matrix
    if (any(!(sg %in% row.range))) {
      stop(sprintf("Subgraph %d contains invalid row indices.", i))
    }
    #check that the subgraphs all contain overlaps (necessary for reconnecting after learning)
    if (!any(sg %in% other.indexs)){
      warning(sprintf("Subgraph %d has no overlaps with the other subgraphs.", i))
    }
  }

  #Check every row index is in at least one subGraph
  if (!all(row.range %in% unique(unlist(modules)))) {
    stop("Not all rows of data are covered by 'modules'.")
  }

  #Check if both functions exist
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

  ##################################
  ###Divide Data and Learn Graphs###
  ##################################

  #divide data based on subGraph indices
  sub.x  <- lapply(modules, function(sg){
    x[sg,]
  })

  #Learn graphs in parallel
  #NOTE: this function can be substituted for any graph learning function that returns and Igraph object
  args <- arg.wrapping.func(sub.x, ...)
  graph.learning.outputs <- foreach::foreach(i = seq_along(args),
                                             .combine = 'list',
                                             .multicombine = TRUE,
                                             .packages = packages.to.each,
                                             .export = export.to.each
  ) %dopar% {
    do.call(graph.learning.func, args[[i]])
  }

  #Separate outputs of graph learning into a list of learned igraph objects and other outputs
  parsed.outputs <- out.parsing.func(graph.learning.outputs)

  if(!all(unlist(lapply(parsed.outputs$learned.graphs , igraph::is_igraph)))){
    stop("The 'learned graphs' within 'parsed outputs' are not all igraph objects.
         Check your graph learning and output parsing functions
         and ensure parsed.outputs$learned.graphs contains a list of igraphs with overlapping node names")
  }

  #################################
  ###Stitch Graphs Back Together###
  #################################

  final.graph <- .connect_subgraphs(parsed.outputs$learned.graphs, keep.all.edges)

  # return sub graphs and final graphs
  return(
    list(
      modular.subgraphs = parsed.outputs$learned.graphs,
      final.graph = final.graph,
      other.outputs = parsed.outputs$other.outputs
    )
  )
}

.connect_subgraphs <- function(sub.graphs, keep.all.edges = FALSE){
  #new graph containing all edges
  full.graph <- do.call(igraph::union, sub.graphs)

  #if keeping all edges, return that graph
  if(keep.all.edges){return(full.graph)}

  #rebuild from adj matrix to properly reorder nodes and remove unwanted attributes
  adj <- as.matrix(igraph::as_adjacency_matrix(full.graph))
  adj <- adj[
    order(as.numeric(stringr::str_remove(rownames(adj), "Node_"))),
    order(as.numeric(stringr::str_remove(colnames(adj), "Node_")))]
  full.graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")

  #remove non overlap edges
  sg.pairs <- utils::combn(seq_along(sub.graphs), m = 2, simplify = F) #get all pairs of subgraphs
  edges.to.remove <- lapply(sg.pairs, function(pr){ #for each pair of subgraphs
    #get index for 'left' and 'right' subgraphs
    left = pr[1]
    right = pr[2]

    #get set of nodes that are in both graphs
    overlapping.nodes <- BiocGenerics::intersect(
      igraph::V(sub.graphs[[left]])$name, igraph::V(sub.graphs[[right]])$name
    )

    #if there are overlapping nodes
    if(length(overlapping.nodes) > 1){
      #get edge list for these nodes
      left.edges <- igraph::as_edgelist(sub.graphs[[left]])
      right.edges <- igraph::as_edgelist(sub.graphs[[right]])

      #if both graphs contain edges
      if(nrow(left.edges) > 0 && nrow(right.edges) > 0){
        left.edges <-left.edges[
          which(
            apply(left.edges, 1, function(row){
              all(row %in% overlapping.nodes)
            })
          )
          ,] #take only edges between overlap nodes

        right.edges <-right.edges[
          which(
            apply(right.edges, 1, function(row){
              all(row %in% overlapping.nodes)
            })
          )
          ,] #for both sets of edges

        #this will return a vector in the edge case in which there is only one edge with overlap nodes
        if(!is.matrix(left.edges) && length(left.edges) == 2){
          left.edges <- matrix(left.edges, ncol = 2, byrow = T)
        }
        if(!is.matrix(right.edges) && length(right.edges) == 2){
          right.edges <- matrix(right.edges, ncol = 2, byrow = T)
        }

        #if these matrixes have at least one row
        if(nrow(left.edges) > 0 && nrow(right.edges) > 0){
          #get the edges that dont occur in both sets
          sym.diff <- dplyr::bind_rows(
            as.data.frame(left.edges) %>% dplyr::anti_join(as.data.frame(right.edges), by = c("V1","V2")),
            as.data.frame(right.edges) %>% dplyr::anti_join(as.data.frame(left.edges), by = c("V1","V2"))
          )
          #return the matrix of edges to remove
          return(as.matrix(sym.diff))
        }

      }
    }
    return(matrix(nrow = 0 , ncol = 2))
  })

  #convert list of matrixs into single matrix containing edge information
  edges.to.remove <-as.matrix(unique(do.call(rbind, edges.to.remove)))

  #remove these edges from combined graoh
  full.graph <- igraph::delete_edges(full.graph, igraph::get.edge.ids(full.graph, as.vector(t(edges.to.remove)), directed = FALSE))

  return(full.graph)
}

.default_arg_wrapper <- function(sub.x,
                                 method  = "B_NW_SL",
                                 global = T,
                                 alpha = 0.05,
                                 fdr.filter = T,
                                 max.fdr = 0.05,
                                 pos.cut = 0,
                                 neg.cut = 0){

  #return a list of packaged args
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

.default_output_parser <- function(graph.learning.outputs){
  return(
    list(
      learned.graphs = graph.learning.outputs,
      other.outputs = NULL
    )
  )
}

############
### TEST ###
############

# # load graph learning func
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/02.LearnGraphs.R")
#
# # create test data
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/00.SimulateGraphs.R")
#
# g <- make_modular_graph()
# x <- sim_graph_data(g, 150)
#
# #get modules from true graph
# mods <- split(igraph::V(g), igraph::V(g)$module)
#
# #create fuzzy module including all second degree neighbors
# fuzzyMods <- lapply(mods, function(mod){
#   sort(unique(unlist(
#     igraph::neighborhood(
#       g,
#       order = 2,
#       nodes = igraph::V(g)[mod],
#       mode = "all",
#       mindist = 0
#     )
#   )))
# })
#
# #get N threads
# nThreads <- min(detectCores() - 1, length(args))
#
# #create and register the cluster
# cl <- makeCluster(nThreads)
# registerDoParallel(cl)
#
# ### Test Function with Defaults
# test <- fastDivideAndConquer(x, fuzzyMods) # basic run
# calc_F1(g, test$final.graph)
#
#
# ### Test function with bd graph learn
# library(BDgraph)
#
# .other_arg_wrapper <- function(sub.x, n.bootstraps, sub.ratio){
#   lapply(sub.x, function(x){
#     list(
#       x = t(x)
#     )
#   })
# }
#
# bdWrapper <- function(x){
#   bdFit <- BDgraph::bdgraph(x)
#   postProb <- BDgraph::plinks(bdFit)
#   adjMat <- BDgraph::select(postProb)
#   igraph::graph_from_adjacency_matrix(adjMat, mode = "undirected") # can change this if I want to test edge directional update for .connect_subgraphs
# }
#
# test <- fastDivideAndConquer(x, fuzzyMods,
#                              graph.learning.func = bdWrapper,
#                              arg.wrapping.func = .other_arg_wrapper,
#                              out.parsing.func = .default_output_parser,
#                              packages.to.each = c("igraph", "BDgraph"),
#                              export.to.each = c("bdWrapper")
# )
# calc_F1(g, test$final.graph)
#
# ### Test function with Ziweis graph learning
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/RSCGGM_wrapper.R")
# setwd("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning")
# devtools::load_all("RSCGGM-main")
#
# .other_arg_wrapper <- function(sub.x, n.bootstraps, sub.ratio){
#   lapply(sub.x, function(x){
#     list(
#       x = t(x), # transform to fit n x p matrix required by silggm
#       self.load = T,
#       n.bootstraps = n.bootstraps,
#       sub.ratio = sub.ratio
#       )
#   })
# }
#
#
# test <- fastDivideAndConquer(x, fuzzyMods,
#                              graph.learning.func = learn_consensus_graph,
#                              arg.wrapping.func = .other_arg_wrapper,
#                              out.parsing.func = .default_output_parser,
#                              packages.to.each = c("igraph"),
#                              export.to.each = c("learn_consensus_graph"),
#                              n.bootstraps = 100, sub.ratio = 0.9
# )
# calc_F1(g, test$final.graph)
#
#
# # Stop the cluster
# stopCluster(cl)


