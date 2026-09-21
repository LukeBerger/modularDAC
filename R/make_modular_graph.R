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

#' @importFrom igraph sample_gnp sample_pa V add_edges simplify disjoint_union

#' @export
make_modular_graph <- function(g.type=c("er", "sf"),
                               n.mods=3, n.nodes=120, n.mod.links=3,
                               no.uncon = TRUE, link.all = TRUE,
                               p.edge = 0.05, power=1, z.appeal=1, ...){

  # match args
  g.type <- match.arg(g.type)

  # calculate the number of nodes per module
  n.nodes <- round(n.nodes / n.mods)

  # create modules and add them to the list
  modules <- lapply(seq_len(n.mods), function(i){
    if(g.type == "er"){
      temp <- igraph::sample_gnp(n = n.nodes,
                                 p = p.edge,
                                 ...)
    }else{
      temp <- igraph::sample_pa(n.nodes,
                                power=power,
                                zero.appeal=z.appeal,
                                directed=FALSE,
                                ...)
    }

    # add module label
    igraph::V(temp)$module <- i

    if(no.uncon){
      # give every unconnected node one edge to a connected node in the module
      no.edge <- which(igraph::degree(temp) == 0)
      if(length(no.edge) > 0){
        targets <- setdiff(seq_len(n.nodes), no.edge)
        if(length(targets) > 0){
          # index into targets (avoids sample()'s single-value shortcut)
          partners <- targets[sample.int(length(targets), length(no.edge), replace = TRUE)]
          temp <- igraph::add_edges(temp, as.vector(rbind(no.edge, partners)))
        }
      }
    }

    temp
  })

  # merge modules
  mg <- do.call(igraph::disjoint_union, modules) # disjoint union: combine modules, preserving the module vertex attribute

  if(link.all){
    # add edges between each pair of modules
    mod.pairs <- utils::combn(seq_len(n.mods), m = 2, simplify = FALSE)
  }else{
    mod.chain <- seq_len(n.mods)
    mod.pairs <- mapply(c, mod.chain, c(mod.chain[-1], mod.chain[1]), SIMPLIFY = FALSE)
  }

  # collect inter-module edges for every pair, then add them all in one call
  mg.module <- igraph::V(mg)$module # module labels are fixed while we only add edges
  new.edges <- unlist(lapply(mod.pairs, function(pair){
    i.nodes <- which(mg.module == pair[1]) # nodes in the first module of the pair
    j.nodes <- which(mg.module == pair[2]) # nodes in the second module of the pair
    # link i and j by adding edges between random nodes
    rbind(
      sample(i.nodes, n.mod.links, replace = TRUE),
      sample(j.nodes, n.mod.links, replace = TRUE)
    )
  }))
  mg <- igraph::add_edges(mg, new.edges)

  # removes any loops and multiples resulting from randomness (will in some cases result in fewer edges than intended)
  mg <- igraph::simplify(mg, remove.multiple=TRUE, remove.loops=TRUE)

  # label nodes
  igraph::V(mg)$name <- paste0("Node_", igraph::V(mg))

  # return graph
  return(mg)
}
