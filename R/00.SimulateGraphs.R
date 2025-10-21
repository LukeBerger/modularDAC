# library(igraph)
# library(netUtils)
# library(BDgraph)
# library(MASS)
# library(dplyr)

# Create graph comprised of scale free modules
make_modular_graph <- function(g.type="er",
                               n.mods=3, n.nodes=120, n.mod.links=3,
                               no.uncon =T, link.all = T,
                               p.edge = 0.05, power=1, z.appeal=1, ...){
  #calculate the number of nodes per module
  n.nodes <- round(n.nodes / n.mods)

  #initialize an empty list of modules
  modules <- vector(mode = "list", length = n.mods)

  #create modules and add them to the list
  for (i in 1:n.mods) {
    if(g.type == "er"){
      temp <-  igraph::sample_gnp(n = n.nodes,
                                  p = p.edge,
                                  ...)
    }else if(g.type == "sf"){
      temp <-  igraph::sample_pa(n.nodes,
                                 power=power,
                                 zero.appeal=z.appeal,
                                 directed=FALSE,
                                 ...)
    }else{
      stop("invalid graph type: submit 'er' or 'sf' ")
    }

    #Add module label
    igraph::V(temp)$module <- i

    if(no.uncon){
      #check for unconnected nodes
      no.edge <-  igraph::V(temp)[igraph::degree(temp) == 0]

      #if any exist
      if(length(no.edge) > 0){
        #give them each an edge
        for(node in no.edge){
          #get module
          temp.mod <- igraph::V(temp)$module[node]

          #get targets
          targets <- igraph::V(temp)[igraph::V(temp)$module == temp.mod]
          targets <- targets[!(targets %in% no.edge)]

          #add edge
          temp <- igraph::add_edges(temp, c(node, sample(targets, 1)))
        }
      }
    }

    #add to list
    modules[[i]] <- temp
  }

  #merge modules
  mg <- BiocGenerics::Reduce("+", modules) #taking advantage the fact that igraphs are just fancy lists

  if(link.all){
    #add edges between each pair of modules
    mod.pairs <- utils::combn(1:n.mods, m = 2, simplify = F)
  }else{
    mod.chain <- c(1:n.mods)
    mod.pairs <- mapply(c, mod.chain, c(mod.chain[-1], mod.chain[1]), SIMPLIFY = FALSE)
  }

  #for each pair
  for (pair in mod.pairs) {
    #get node indexes based on modules
    i.nodes <- BiocGenerics::which(igraph::V(mg)$module == pair[1]) #first half of pair
    j.nodes <- BiocGenerics::which(igraph::V(mg)$module == pair[2]) #second half of nodes

    #link i and j by adding edges between random nodes
    mg <- igraph::add_edges(mg,
                            as.vector(
                              rbind(
                                sample(i.nodes, n.mod.links, replace =  T),
                                sample(j.nodes, n.mod.links, replace =  T)
                              )
                            )
    )
  }

  #removes any loops and multiples resulting from randomness (will in some cases result in less edges then intended)
  mg <- igraph::simplify(mg, remove.multiple=TRUE, remove.loops=TRUE)

  #label nodes
  igraph::V(mg)$name <- paste0("Node_", igraph::V(mg))

  #return graph
  return(mg)
}

## Create a hybrid graph by combining multiple modular graphs
## Args: number of subgraphs,
## number of nodes, number of modules in each subgraph
## average number of edges to other nodes with a module
## number of links between models in subgraphs, number of links between subgraphs
make_submodular_graph <- function(g.type="er",
                                  n.mods=3, n.mod.links=3,
                                  n.sg=3, n.sg.nodes=120, n.sg.links=1,
                                  no.uncon =T, link.all = T,
                                  p.edge = 0.05, power=1, z.appeal=1, ...){
  #initialize and empty list of subgraphs
  subgraphs <- vector(mode = "list", length = n.sg)

  #create subgraphs and add them to the list
  for (i in 1:n.sg) {
    #create modular graph: er or sf
    temp <- make_modular_graph(g.type=g.type,
                               n.mods=n.mods, n.nodes=n.sg.nodes, n.mod.links=n.mod.links,
                               no.uncon =no.uncon, link.all = link.all,
                               p.edge = p.edge, power=power, z.appeal=z.appeal) #, ...
    #add subgraph label
    igraph::V(temp)$subgraph <- i

    #fix module lables
    mod.names <- paste0(i, ".", letters[1:n.mods]) #cant have more than 26 modules or this breaks down
    temp <- .label_mods(temp, n.sg.nodes, n.mods, mod.names)

    #and node names
    igraph::V(temp)$name <- paste0("Subgraph_", as.character(i) , "_",igraph::V(temp)$name)

    if(no.uncon){
      #check for unconnected nodes
      no.edge <-  igraph::V(temp)[igraph::degree(temp) == 0]

      #if any exist
      if(length(no.edge) > 0){
        #give them each an edge
        for(node in no.edge){
          #get module
          temp.mod <- igraph::V(temp)$module[node]

          #get targets
          targets <- igraph::V(temp)[igraph::V(temp)$module == temp.mod]
          targets <- targets[!(targets %in% no.edge)]

          #add edge
          temp <- igraph::add_edges(temp, c(node, sample(targets, 1)))
        }
      }
    }
    #add to list
    subgraphs[[i]] <- temp
  }

  #merge subgraphs
  hybrid <- do.call(igraph::union, subgraphs)
  igraph::V(hybrid)$name <- unlist(lapply(subgraphs, function(x) igraph::V(x)$name))
  igraph::V(hybrid)$subgraph <- unlist(lapply(subgraphs, function(x) igraph::V(x)$subgraph))

  if(link.all){
    #add edges between each pair of modules
    mod.pairs <- utils::combn(1:n.sg, m = 2, simplify = F)
  }else{
    mod.chain <- c(1:n.sg)
    mod.pairs <- mapply(c, mod.chain, c(mod.chain[-1], mod.chain[1]), SIMPLIFY = FALSE)
  }
  for (pair in mod.pairs) {
    #get node indexes based on subgraphs
    i.nodes <- BiocGenerics::which(igraph::V(hybrid)$subgraph == pair[1]) #first half of pair
    j.nodes <- BiocGenerics::which(igraph::V(hybrid)$subgraph == pair[2]) #second half of nodes
    #link i and j by adding edges between random nodes
    hybrid <- igraph::add_edges(hybrid,
                                as.vector(
                                  rbind(
                                    sample(i.nodes, n.sg.links, replace =  T),
                                    sample(j.nodes, n.sg.links, replace =  T)
                                  )
                                )
    )
  }

  hybrid <- igraph::simplify(hybrid, remove.multiple=TRUE, remove.loops=TRUE)
  return(hybrid)
}

## Add labels to modular graphs, given the parameters used to generate the graph
## (the number of modules and the total number of nodes)
.label_mods <- function(g, n.nodes, n.mods, mod.names){
  label.start <- 1
  for(i in 1:n.mods){
    igraph::V(g)[label.start:((n.nodes / n.mods)*i)]$module <- mod.names[i]
    label.start <- (n.nodes / n.mods)*i + 1
  }
  return(g)
}


## Wrapper of netUtils sample_lft
## Takes inputs and returns graph labled to match the above format
##      I.E containing Node_Names and $modules
make_lft <- function(n = 120,        #number of nodes
                     tau1 = 3,       #power-law exponent for degree distribution
                     tau2 = 2,       #power-law exponent for community size distribution
                     mu = 0.08,      #mixing parameter
                     average.degree = 6,   #required: typical node degree
                     max.degree = 10,      #required: highest possible node  degree
                     min.community = 30,   # min communitu size
                     max.community = 50){

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

  #label
  igraph::V(g)$name <- paste0("Node_", igraph::V(g))
  igraph::V(g)$module <- igraph::V(g)$membership

  return(g)
}


## Generates simulated data for a graph
## Used BD Graph to learn covariance matrix and mass to generate samples
## Takes graph, vector of node means, and number of samplesto generate
## returns p x n matrix of data
sim_graph_data <- function(g, n.samples, mean.vec = NULL){
  #if mean.vec isnt set, use vector of all 10s
  if(is.null(mean.vec)){
    #set as a vector of 10s
    mean.vec = rep(10, length(g))
  }

  bdg <- mod.bdgraph.sim(n=1, #makes a single samples worth of data
                         graph=as.matrix(igraph::as_adjacency_matrix(g)),
                         type="Gaussian",
                         mean = mean.vec[1]
  )
  sigma.mat <- bdg$sigma #get the covariance matrix from bdg.graph

  #make n samples with MASS
  sim.data <- t(MASS::mvrnorm(n = n.samples, mu = mean.vec, Sigma =  sigma.mat))

  # #format as eset
  # sim.data <- Biobase::ExpressionSet(t(sim.data))

  #add row and column names
  colnames(sim.data) <- paste("Sample", 1:ncol(sim.data), sep="_")
  rownames(sim.data) <- igraph::V(g)$name

  #return data
  return(sim.data)
}


## Uses visNetwork to plot interactive modular graphs
## Requires a graph with a module column on each vertex
## Can also plot hybrid graphs with a subgraph column
## Names are adjustable to account for graphs with multiple recorded modules
modular_plot <- function(g,
                         module.name = "module",
                         subgraphs = FALSE,
                         subname = "subgraph"){
  #get edges matrix
  edge.list <- igraph::as_edgelist(g)
  edges <- data.frame(
    from = as.character(edge.list[, 1]),
    to = as.character(edge.list[, 2])
  )

  #if there are subgraphs, color by module select by subgraph
  if(subgraphs){
    #get nodes matrix
    nodes <- data.frame(
      id = igraph::V(g)$name %||% as.character(igraph::V(g)),
      label = igraph::V(g)$label %||% as.character(igraph::V(g)),
      group = as.factor(igraph::vertex.attributes(g)[[module.name]]),  #color by module
      subgraph = as.factor(igraph::vertex.attributes(g)[[subname]]) #select by subgraph
    )

    #Plot with visNetwork
    visNetwork::visNetwork(nodes, edges) %>%
      visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = list(enabled = TRUE, useLabels = TRUE),  # ID/label selector
        selectedBy = list(variable = "subgraph", multiple = TRUE) #subgraph selector
      ) %>%
      visNetwork::visLegend() %>%
      visNetwork::visEdges(smooth=FALSE) %>%
      visNetwork::visInteraction(hover = TRUE) %>%
      visNetwork::visPhysics(stabilization = FALSE) %>%
      visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
      visNetwork::visLayout(randomSeed = 123)

    #else color and select by modules
  }else{
    #get nodes matrix
    nodes <- data.frame(
      id = igraph::V(g)$name %||% as.character(igraph::V(g)),
      label = igraph::V(g)$label %||% as.character(igraph::V(g)),
      group = as.factor(igraph::vertex.attributes(g)[[module.name]]), #color by module
      module = as.factor(igraph::vertex.attributes(g)[[module.name]]) #select by module
    )

    #Plot with visNetwork
    visNetwork::visNetwork(nodes, edges) %>%
      visNetwork::visOptions(
        highlightNearest = TRUE,
        nodesIdSelection = list(enabled = TRUE, useLabels = TRUE),  # ID/label selector
        selectedBy = list(variable = "module", multiple = TRUE)      # Group/module selector
      ) %>%
      visNetwork::visLegend() %>%
      visNetwork::visEdges(smooth=FALSE) %>%
      visNetwork::visInteraction(hover = TRUE) %>%
      visNetwork::visPhysics(stabilization = FALSE) %>%
      visNetwork::visIgraphLayout(layout = "layout_with_fr") %>%
      visNetwork::visLayout(randomSeed = 123)
  }
}


## bdgraph function with two lines removed which caused a bug when using a
## adjency matrix as the graph parameter
mod.bdgraph.sim <- function (p = 10, graph = "random", n = 0, type = "Gaussian",
                             prob = 0.2, size = NULL, mean = 0, class = NULL, cut = 4,
                             b = 3, D = diag(p), K = NULL, sigma = NULL, q = exp(-1),
                             beta = 1, vis = FALSE, rewire = 0.05, range.mu = c(3, 5),
                             range.dispersion = c(0.01, 0.1), nu = 1)
{
  if (p < 2)
    stop("'p' must be greater than 1")
  if ((prob < 0) | (prob > 1))
    stop("'prob' must be between ( 0, 1 )")
  if (cut < 2)
    stop("'cut' must be greater than 1")
  if (b <= 2)
    stop("'b' must be greater than 2")
  if ((rewire < 0) | (rewire > 1))
    stop("'rewire' must be between ( 0, 1 )")
  if (length(range.mu) != 2)
    stop("'range.mu' must be a vector with length 2")
  if (length(range.dispersion) != 2)
    stop("'range.dispersion' must be a vector with length 2")
  if (inherits(graph, "graph"))
    graph = unclass(graph)
  if (is.matrix(graph) & is.matrix(K))
    if (nrow(graph) != nrow(K))
      stop("'graph' and 'K' have non-conforming size")
  if (!is.null(size))
    if ((sum(size) < 0) | (sum(size) > (p * (p - 1)/2)))
      stop("'size' must be between ( 0, p*(p-1)/2 )")
  if (is.matrix(K)) {
    if (!isSymmetric(K))
      stop("'K' must be a positive definite matrix")
    graph <- "fixed"
    p <- nrow(K)
  }
  if (type == "normal")
    type = "Gaussian"
  if (type == "non-normal")
    type = "non-Gaussian"
  if ((type == "categorical") & (cut == 2))
    type = "binary"
  if (is.matrix(graph)) {
    if (!isSymmetric(graph))
      stop("'graph' must be symmetric matrix")
    p = nrow(graph)
    if (all(graph %in% c(0,1))) # if ((graph != 0) && (graph != 1))
      #   stop("Elements of matrix 'graph' must be 0 or 1")
      G <- graph
    graph <- "fixed"
  }
  if (!any(graph == c("fixed", "AR1", "AR2", "circle")))
    G <- BDgraph::graph.sim(p = p, graph = graph, prob = prob,
                            size = size, class = class, rewire = rewire)
  if (graph == "AR1") {
    sigma = matrix(0, p, p)
    for (i in 1:(p - 1)) for (j in (i + 1):p) sigma[i, j] = (0.7)^abs(i -
                                                                        j)
    sigma = sigma + t(sigma) + diag(p)
    K = solve(sigma)
    G = 1 * (abs(K) > 0.02)
  }
  if (graph == "AR2") {
    K = stats::toeplitz(c(1, 0.5, 0.25, rep(0, p - 3)))
    G = 1 * (abs(K) > 0.02)
  }
  if (graph == "circle") {
    K <- stats::toeplitz(c(1, 0.5, rep(0, p - 2)))
    K[1, p] <- 0.4
    K[p, 1] <- 0.4
    G = 1 * (abs(K) > 0.02)
  }
  if (n != 0) {
    if (!is.null(sigma))
      K <- solve(sigma)
    if (is.matrix(K)) {
      G = 1 * (abs(K) > 0.02)
      diag(G) = 0
      if (is.null(sigma))
        sigma = stats::cov2cor(solve(K))
    }
    else {
      if (!isSymmetric(D))
        stop("'D' must be a positive definite matrix")
      Ti = chol(solve(D))
      diag(G) = 0
      K = matrix(0, p, p)
      threshold = 1e-08
      K = BDgraph::rgwish(adj = G, b = b, D = Ti, threshold = threshold)
      # K = matrix(result$K, p, p)
      sigma = stats::cov2cor(solve(K))
      K = solve(sigma)
    }
    d <- BDgraph::rmvnorm(n = n, mean = mean, sigma = sigma)
    not.cont = numeric(p)
    if (type == "mixed") {
      ps = floor(p/5)
      col_number <- c(1:ps)
      prob <- stats::pnorm(d[, col_number])
      d[, col_number] <- stats::qpois(p = prob, lambda = 10)
      not.cont[1:ps] = 1
      col_number <- c((ps + 1):(2 * ps))
      prob <- stats::pnorm(d[, col_number])
      d[, col_number] <- stats::qpois(p = prob, lambda = 2)
      not.cont[c((ps + 1):(2 * ps))] = 1
      col_number <- c((2 * ps + 1):(3 * ps))
      prob <- stats::pnorm(d[, col_number])
      d[, col_number] <- stats::qexp(p = prob, rate = 10)
      col_number <- c((3 * ps + 1):(4 * ps))
      prob <- stats::pnorm(d[, col_number])
      d[, col_number] <- stats::qbinom(p = prob, size = 1,
                                       prob = 0.5)
      not.cont[c((3 * ps + 1):(4 * ps))] = 1
    }
    if (type == "non-Gaussian") {
      prob <- stats::pnorm(d)
      d <- stats::qexp(p = prob, rate = 10)
    }
    if (type == "t") {
      tau_gamma = stats::rgamma(n, shape = nu/2, rate = nu/2)
      d = mean + d/sqrt(tau_gamma)
    }
    if (type == "alternative-t") {
      taugamma = stats::rgamma(n * p, shape = nu/2, rate = nu/2)
      taugamma = matrix(taugamma, n, p)
      d = mean + d/sqrt(taugamma)
    }
    if (type == "categorical") {
      not.cont[1:p] = 1
      runif_m = matrix(stats::runif(cut * p), nrow = p,
                       ncol = cut)
      marginals = apply(runif_m, 1, function(x) {
        stats::qnorm(cumsum(x/sum(x))[-length(x)])
      })
      if (cut == 2)
        marginals = matrix(marginals, nrow = 1, ncol = p)
      for (j in 1:p) {
        breaks <- c(min(d[, j]) - 1, marginals[, j],
                    max(d[, j]) + 1)
        d[, j] <- as.integer(cut(d[, j], breaks = breaks,
                                 right = FALSE))
      }
      d = d - 1
    }
    if ((type == "dweibull") | (type == "dw")) {
      if (length(q) == 1)
        q = rep(q, time = p)
      if (length(beta) == 1)
        beta = rep(beta, time = p)
      if (is.vector(q) && (length(q) != p))
        stop("'q', as a vector, has non-conforming size with 'p'")
      if (is.vector(beta) && (length(beta) != p))
        stop("'beta', as a vector, has non-conforming size with 'p'")
      if (is.matrix(q) && any(dim(q) != c(n, p)))
        stop("'q', as a matrix, has non-conforming size with 'n' and 'p'")
      if (is.matrix(beta) && any(dim(beta) != c(n, p)))
        stop("'beta', as a matrix, has non-conforming size with 'n' and 'p'")
      not.cont[1:p] = 1
      Y_data <- matrix(c(0, 1), nrow = n, ncol = p)
      Z = tmvtnorm::rtmvnorm(n = n, mean = rep(mean, p),
                             sigma = sigma, lower = rep(-5, length = p),
                             upper = rep(5, length = p))
      pnorm_Z = stats::pnorm(Z)
      if (is.matrix(q) && is.matrix(beta))
        for (j in 1:p) Y_data[, j] = BDgraph::qdweibull(pnorm_Z[,
                                                                j], q = q[, j], beta = beta[, j], zero = TRUE)
      if (is.vector(q) && is.vector(beta))
        for (j in 1:p) Y_data[, j] = BDgraph::qdweibull(pnorm_Z[,
                                                                j], q = q[j], beta = beta[j], zero = TRUE)
      if (any(apply(Y_data, 2, function(x) {
        all(x %in% 0:1)
      })))
        cat(" Some of the variables are binary \n")
      d = Y_data
    }
    if ((type == "nbinom") | (type == "NB")) {
      not.cont[1:p] = 1
      Y.star <- matrix(c(0, 1), nrow = n, ncol = p)
      while (any(apply(Y.star, 2, function(x) {
        all(x %in% 0:1)
      })) == TRUE) {
        d = tmvtnorm::rtmvnorm(n = n, mean = rep(mean,
                                                 p), sigma = sigma, lower = rep(-5, length = p),
                               upper = rep(5, length = p))
        mu = stats::runif(n = p, min = range.mu[1],
                          max = range.mu[2])
        size = stats::runif(n = p, min = range.dispersion[1],
                            max = range.dispersion[2])
        for (j in 1:p) Y.star[, j] = stats::qnbinom(stats::pnorm(d[,
                                                                   j]), size = size[j], mu = mu[j], lower.tail = TRUE,
                                                    log.p = FALSE)
      }
      d = Y.star
    }
    if ((type == "pois") | (type == "count")) {
      not.cont[1:p] = 1
      Y.star = matrix(c(0, 1), nrow = n, ncol = p)
      while (any(apply(Y.star, 2, function(x) {
        all(x %in% 0:1)
      })) == TRUE) {
        d = tmvtnorm::rtmvnorm(n = n, mean = rep(mean,
                                                 p), sigma = sigma, lower = rep(-5, length = p),
                               upper = rep(5, length = p))
        lambda = stats::runif(n = p, min = range.mu[1],
                              max = range.mu[2])
        for (j in 1:p) Y.star[, j] = stats::qpois(stats::pnorm(d[,
                                                                 j]), lambda = lambda[j], lower.tail = TRUE,
                                                  log.p = FALSE)
      }
      d = Y.star
    }
  }
  if (n != 0) {
    if (type != "dw") {
      simulation <- list(G = G, graph = graph, data = d,
                         sigma = sigma, K = K, type = type, not.cont = not.cont)
    }
    else {
      simulation <- list(G = G, graph = graph, data = d,
                         sigma = sigma, K = K, type = type, not.cont = not.cont,
                         beta = beta, q = q)
    }
  }
  else {
    simulation <- list(G = G, graph = graph)
  }
  if (vis == TRUE)
    BDgraph::plot.graph(G, main = "Graph structure")
  class(simulation) <- "sim"
  return(simulation)
}


############
### TEST ###
############

# # make graphs
# er <- make_modular_graph()
# hybrid <- make_submodular_graph()
#
# #plot
# modular_plot(er)
#
# # sim data
# x <- sim_graph_data(er, n.samples = 10)
