#' Expand true modules to fuzzy (overlapping) modules based on second-order graph neighborhood
#' @param m a module S4 object containing module membership assignments
#' @param g an igraph object with a 'module' vertex attribute defining ground-truth module membership for each node

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V neighborhood

#' @export
true_fuzzy <- function(m, g){
  node.names <- igraph::V(g)$name
  # collect all nodes within 2 graph hops of each module
  f.index.list <- lapply(m@index.list, function(idx){
    sort(unique(unlist(
      igraph::neighborhood(g, order = 2, nodes = idx, mode = "all", mindist = 0)
    )))
  })
  # build fuzzy module object with expanded index list; the original (pre-growth)
  # nodes are the core nodes this module owns, the recruited neighbours are auxiliary
  methods::new("module",
               source = "True Module Fuzzy",
               overlapping = TRUE,
               index.list = f.index.list,
               name.list = lapply(f.index.list, function(idx) node.names[idx]),
               core.list = m@index.list
  )
}
