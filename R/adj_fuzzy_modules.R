#' Expand non-overlapping modules to fuzzy (overlapping) modules using a thresholded WGCNA adjacency matrix
#'
#' The nodes of the input partition become the \emph{core} (owned) nodes of each
#' fuzzy module; the recruited neighbours are auxiliary nodes carried along so
#' the module can condition on each core node's Markov blanket. Growth stops at
#' \code{ratio}/\code{max.size}, or earlier if \code{min.adj} is supplied and no
#' external feature reaches that adjacency (an approximate blanket-closure cutoff).
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param adj a p x p numeric adjacency matrix (e.g. from WGCNA::adjacency) defining pairwise feature similarity
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param max.size an integer, the maximum number of nodes allowed in any fuzzy module
#' @param ratio a numeric, the maximum ratio of fuzzy module size to original module size
#' @param min.adj a numeric or NULL; if set, only features whose maximum adjacency to a core node exceeds this value are recruited (blanket-closure cutoff)

#' @return a module object

#' @importFrom methods new

#' @export
adj_fuzzy_modules <- function(x, adj, input.modules, max.size, ratio, min.adj = NULL){
  # growth only defines ownership cleanly when the input is a partition
  if (input.modules@overlapping) {
    stop("input.modules must be a non-overlapping partition (its nodes become the core nodes of each fuzzy module).")
  }
  # check if any modules are too large
  if (any(lengths(input.modules@index.list) > max.size)) {
    stop("Some modules are too large, increase max.size.")
  }

  # create new index list (list of nodes in each fuzzy module by index)
  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod.nodes <- input.modules@index.list[[m]]
    # get number of required fuzzy nodes (capped by max.size and by nodes available)
    f.size <- length(mod.nodes) * ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- min(f.size - length(mod.nodes), nrow(x) - length(mod.nodes))
    if(n.fuzzy.nodes <= 0){ return(sort(mod.nodes)) }

    # get matrix of nodes in module adj with nodes outside module
    in.out.adj <- adj[mod.nodes, -mod.nodes, drop = FALSE]

    # get max adj of node outside module with a node inside module
    out.max <- apply(in.out.adj, 2, max)

    # rank and select fuzzy nodes based on adjacency to the module
    adj.rank <- order(out.max, decreasing = TRUE)[seq_len(n.fuzzy.nodes)]
    if(!is.null(min.adj)){ adj.rank <- adj.rank[out.max[adj.rank] >= min.adj] } # closure cutoff
    fuzzy.nodes <- which(rownames(x) %in% colnames(in.out.adj)[adj.rank])

    # return fuzzy module combining original (core) and recruited (auxiliary) nodes
    sort(c(mod.nodes, fuzzy.nodes))
  })

  # convert to module object and return; original nodes are the owned core nodes
  fuzzy.mods <- methods::new("module",
                             source = paste("adj_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(m){rownames(x)[m]}),
                             core.list = input.modules@index.list
  )
  return(fuzzy.mods)
}
