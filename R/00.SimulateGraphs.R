#####################
### Create Graphs ###
#####################

#' Create a modular graph comprised of multiple modules using either Erdos-Renyi or scale-free structure within modules
#' @param g.type a character, the graph type; 'er' for Erdos-Renyi or 'sf' for scale-free
#' @param n.mods an integer, the number of modules in the output graph
#' @param n.nodes an integer, the total number of nodes in the graph
#' @param n.mod.links an integer, the number of edges between each pair of modules
#' @param no.uncon a logical, if TRUE nodes with zero edges are not permitted (each isolated node is assigned one edge)
#' @param link.all a logical, if TRUE all module pairs receive inter-module edges; if FALSE modules are linked in a chain
#' @param p.edge a numeric between 0 and 1, the probability of an edge between any two nodes within a module (Erdos-Renyi only)
#' @param power an integer, the power parameter for scale-free graph generation via preferential attachment
#' @param z.appeal an integer, the zero-appeal parameter controlling the attractiveness of isolated nodes in scale-free generation
#' @param ... additional arguments passed to the underlying igraph graph-building functions

#' @return an igraph object

#' @importFrom igraph sample_gnp sample_pa V add_edges simplify

#' @export
make_modular_graph <- function(g.type=c("er", "sf"),
                               n.mods=3, n.nodes=120, n.mod.links=3,
                               no.uncon = TRUE, link.all = TRUE,
                               p.edge = 0.05, power=1, z.appeal=1, ...){

  # match args
  g.type <- match.arg(g.type)

  # calculate the number of nodes per module
  n.nodes <- round(n.nodes / n.mods)

  # initialize an empty list of modules
  modules <- vector(mode = "list", length = n.mods)

  # create modules and add them to the list
  for (i in seq_len(n.mods)) {
    if(g.type == "er"){
      temp <- igraph::sample_gnp(n = n.nodes,
                                  p = p.edge,
                                  ...)
    }else if(g.type == "sf"){
      temp <- igraph::sample_pa(n.nodes,
                                 power=power,
                                 zero.appeal=z.appeal,
                                 directed=FALSE,
                                 ...)
    }else{
      stop("invalid graph type: submit 'er' or 'sf' ")
    }

    # add module label
    igraph::V(temp)$module <- i

    if(no.uncon){
      # check for unconnected nodes
      no.edge <- igraph::V(temp)[igraph::degree(temp) == 0]

      # if any exist
      if(length(no.edge) > 0){
        # give them each an edge
        for(node in no.edge){
          # get module
          temp.mod <- igraph::V(temp)$module[node]

          # get targets
          targets <- igraph::V(temp)[igraph::V(temp)$module == temp.mod]
          targets <- targets[!(targets %in% no.edge)]

          # add edge
          temp <- igraph::add_edges(temp, c(node, sample(targets, 1)))
        }
      }
    }

    # add to list
    modules[[i]] <- temp
  }

  # merge modules
  mg <- Reduce("+", modules) # taking advantage of the fact that igraphs are just fancy lists

  if(link.all){
    # add edges between each pair of modules
    mod.pairs <- utils::combn(seq_len(n.mods), m = 2, simplify = FALSE)
  }else{
    mod.chain <- seq_len(n.mods)
    mod.pairs <- mapply(c, mod.chain, c(mod.chain[-1], mod.chain[1]), SIMPLIFY = FALSE)
  }

  # for each pair
  mg.module <- igraph::V(mg)$module # module labels are fixed while we only add edges
  for (pair in mod.pairs) {
    # get node indexes based on modules
    i.nodes <- which(mg.module == pair[1]) # nodes in the first module of the pair
    j.nodes <- which(mg.module == pair[2]) # nodes in the second module of the pair

    # link i and j by adding edges between random nodes
    mg <- igraph::add_edges(mg,
                            as.vector(
                              rbind(
                                sample(i.nodes, n.mod.links, replace = TRUE),
                                sample(j.nodes, n.mod.links, replace = TRUE)
                              )
                            )
    )
  }

  # removes any loops and multiples resulting from randomness (will in some cases result in fewer edges than intended)
  mg <- igraph::simplify(mg, remove.multiple=TRUE, remove.loops=TRUE)

  # label nodes
  igraph::V(mg)$name <- paste0("Node_", igraph::V(mg))

  # return graph
  return(mg)
}



#' Wrapper for netUtils::sample_lfr that generates a Lancichinetti-Fortunato-Radicchi (LFR) benchmark graph and labels nodes
#' @param n an integer, the number of nodes in the graph
#' @param tau1 a numeric, the power-law exponent for the degree distribution
#' @param tau2 a numeric, the power-law exponent for the community size distribution
#' @param mu a numeric between 0 and 1, the mixing parameter controlling the fraction of inter-community edges per node
#' @param average.degree a numeric, the average node degree in the graph
#' @param max.degree an integer, the maximum node degree allowed in the graph
#' @param min.community an integer, the minimum number of nodes in each community
#' @param max.community an integer, the maximum number of nodes in each community

#' @return an igraph object

#' @importFrom igraph V

#' @export
make_lfr <- function(n = 120,
                     tau1 = 3,
                     tau2 = 2,
                     mu = 0.08,
                     average.degree = 6,
                     max.degree = 10,
                     min.community = 30,
                     max.community = 50){
  if (!requireNamespace("netUtils", quietly = TRUE)) {
    stop("Package netUtils is required. Install with: install.packages('netUtils')", call. = FALSE)
  }
  # generate LFR benchmark graph
  g <- netUtils::sample_lfr(
    n = n,
    tau1 = tau1,
    tau2 = tau2,
    mu = mu,
    average_degree = average.degree,
    max_degree = max.degree,
    min_community = min.community,
    max_community = max.community
  )

  # assign sequential node names and copy community membership to module attribute
  igraph::V(g)$name <- paste0("Node_", igraph::V(g))
  igraph::V(g)$module <- igraph::V(g)$membership

  return(g)
}


#####################
### Simulate Data ###
#####################

#' Simulate gene expression data for n.samples samples based on the covariance structure implied by a graph
#' @param g an igraph object defining the dependency structure between features
#' @param n.samples an integer, the number of samples to generate
#' @param mean.vec a numeric vector of length p, the mean expression value for each feature; defaults to a zero vector if NULL
#' @param b a numeric (> 2), the degrees-of-freedom parameter of the G-Wishart distribution used to sample the precision matrix (passed to BDgraph::rgwish)

#' @return a p x n data matrix

#' @importFrom igraph as_adjacency_matrix V

#' @export
sim_graph_data <- function(g, n.samples, mean.vec = NULL, b = 3){
  if (!requireNamespace("BDgraph", quietly = TRUE)) {
    stop("Package BDgraph is required. Install with: install.packages('BDgraph')", call. = FALSE)
  }
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package MASS is required. Install with: install.packages('MASS')", call. = FALSE)
  }
  # if mean.vec is not set, default to a zero vector
  if(is.null(mean.vec)){
    # zero mean for all features
    mean.vec <- rep(0, length(g))
  }
  # derive the covariance matrix implied by the graph structure: sample a
  # precision matrix from the G-Wishart distribution constrained to the graph
  # (forcing zeros at non-edges), then convert it to a correlation matrix
  adj.mat <- as.matrix(igraph::as_adjacency_matrix(g))
  diag(adj.mat) <- 0
  K <- BDgraph::rgwish(adj = adj.mat, b = b, D = diag(nrow(adj.mat)), threshold = 1e-8)
  sigma.mat <- stats::cov2cor(solve(K)) # correlation matrix implied by the graph

  # make n samples with MASS
  sim.data <- t(MASS::mvrnorm(n = n.samples, mu = mean.vec, Sigma =  sigma.mat))

  # # format as eset
  # sim.data <- Biobase::ExpressionSet(t(sim.data))

  # add row and column names
  colnames(sim.data) <- paste("Sample", seq_len(ncol(sim.data)), sep = "_")
  rownames(sim.data) <- igraph::V(g)$name

  # return data
  return(sim.data)
}

#' Plot a modular graph with nodes colored by module membership using visNetwork
#' @param g an igraph object whose nodes have a module vertex attribute
#' @param module.name a character, the name of the vertex attribute storing module membership (default 'module')

#' @return a visNetwork HTML widget
#' @importFrom dplyr %>%

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

  # colour and filter nodes by module membership
  nodes <- data.frame(
    id = igraph::V(g)$name %||% as.character(igraph::V(g)),
    label = igraph::V(g)$label %||% as.character(igraph::V(g)),
    group = as.factor(igraph::vertex.attributes(g)[[module.name]]), # color by module
    module = as.factor(igraph::vertex.attributes(g)[[module.name]]) # select by module
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

