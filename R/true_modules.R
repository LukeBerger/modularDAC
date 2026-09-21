#' Extract true module assignments from a simulated graph
#' @param g an igraph object with a 'module' vertex attribute defining ground-truth module membership for each node

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V

#' @keywords internal
#' @note Benchmarking helper, kept on this branch only to back the test
#'   suite. It is not part of the user-facing API; the exported version
#'   lives on the dev branch.
true_modules <- function(g){
  # build module object from ground-truth module vertex attribute
  module <- igraph::V(g)$module
  index.list <- split(seq_len(igraph::vcount(g)), module)
  methods::new("module",
                source = "True Modules",
                overlapping = FALSE,
                index.vector = module,
                index.list = index.list,
                name.list = split(igraph::V(g)$name, module),
                # a non-overlapping partition owns all of its own nodes
                core.list = index.list
  )
}
