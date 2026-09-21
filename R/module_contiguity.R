#' Measure how self-contained each module is within a graph
#' @param g an igraph object whose node names match the node names stored in test.module
#' @param test.module a module S4 object containing module membership assignments

#' @return a list with two elements: 'overall.contiguity', a single score from 0 to 100 giving the edge-weighted percentage of module-incident edges that stay within their module; and 'module.edges', a data frame with one row per module (row names are the module labels) holding 'within.edges' (edges with both endpoints in the module), 'between.edges' (edges with exactly one endpoint in the module) and 'percent.within' (100 * within / (within + between))

#' @importFrom igraph is_igraph as_edgelist V
#' @importFrom methods is

#' @export
module_contiguity <- function(g, test.module){
  # --- validate inputs ---
  if(!igraph::is_igraph(g)){
    stop("'g' must be an igraph object.")
  }
  if(!methods::is(test.module, "module")){
    stop("'test.module' must be a 'module' object.")
  }
  # module membership is matched to the graph by node name
  v.names <- igraph::V(g)$name
  module.names <- test.module@name.list
  if(!all(unlist(module.names) %in% v.names)){
    stop("'test.module' contains node names that are not present in 'g'.")
  }

  # module labels used for the output rows
  labels <- names(module.names)
  if(is.null(labels)) labels <- as.character(seq_along(module.names))

  # graph edges as integer endpoint pairs (each undirected edge appears once)
  edges <- igraph::as_edgelist(g, names = FALSE)

  # --- per-module within / between edge counts ---
  counts <- vapply(module.names, function(nm){
    # membership indicator over the graph's vertices
    member <- v.names %in% nm
    if(nrow(edges) == 0){
      return(c(within = 0, between = 0))
    }
    a.in <- member[edges[, 1]]
    b.in <- member[edges[, 2]]
    c(within  = sum(a.in & b.in),       # both endpoints inside the module
      between = sum(xor(a.in, b.in)))   # exactly one endpoint inside the module
  }, numeric(2))

  within.edges   <- counts["within", ]
  between.edges  <- counts["between", ]
  total.edges    <- within.edges + between.edges
  percent.within <- ifelse(total.edges > 0, 100 * within.edges / total.edges, NA_real_)

  module.edges <- data.frame(
    within.edges   = within.edges,
    between.edges  = between.edges,
    percent.within = percent.within,
    row.names      = labels,
    stringsAsFactors = FALSE
  )

  # --- overall contiguity: edge-weighted fraction of module edges kept internal ---
  overall.contiguity <- if(sum(total.edges) > 0){
    100 * sum(within.edges) / sum(total.edges)
  }else{
    NA_real_
  }

  return(
    list(
      overall.contiguity = overall.contiguity,
      module.edges       = module.edges
    )
  )
}
