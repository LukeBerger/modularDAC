#' Plot a modular graph with nodes colored by module membership using visNetwork
#' @param g an igraph object whose nodes have a module vertex attribute
#' @param module.name a character, the name of the vertex attribute storing module membership (default 'module')

#' @return a visNetwork HTML widget
#' @importFrom dplyr %>%
#' @importFrom rlang %||%

#' @export
modular_plot <- function(g,
                         module.name = "module"){
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("Package visNetwork is required. Install with: install.packages('visNetwork')", call. = FALSE)
  }
  # build the edge data frame
  edge.list <- igraph::as_edgelist(g)
  edges <- data.frame(
    from = as.character(edge.list[, 1]),
    to = as.character(edge.list[, 2])
  )

  # module membership for every node (used for both colouring and selection)
  mod <- as.factor(igraph::vertex.attributes(g)[[module.name]])

  # colour and filter nodes by module membership
  nodes <- data.frame(
    id = igraph::V(g)$name %||% as.character(igraph::V(g)),
    label = igraph::V(g)$label %||% as.character(igraph::V(g)),
    group = mod, # color by module
    module = mod # select by module
  )

  # plot with visNetwork
  plot <- visNetwork::visNetwork(nodes, edges) %>%
    visNetwork::visOptions(
      highlightNearest = TRUE,
      nodesIdSelection = list(enabled = TRUE, useLabels = TRUE), # iD/label selector
      selectedBy = list(variable = "module", multiple = TRUE) # module selector
    ) %>%
    visNetwork::visLegend() %>%
    visNetwork::visEdges(smooth = FALSE) %>%
    visNetwork::visInteraction(hover = TRUE) %>%
    visNetwork::visPhysics(stabilization = FALSE) %>%
    visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
    visNetwork::visLayout(randomSeed = 123)

  return(plot)
}
