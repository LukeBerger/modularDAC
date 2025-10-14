# library(igraph)
# library(dplyr)
# library(BiocGenerics)


## function with a series of helpers to get the percentage of nodes
## matching between two sets of modules
## takes a graph, and the names of two modules (attributes of the vertexs)

percent_module_match <- function(g, m1, m2){
  #get modules node sets
  m1.nodes <- split(as.numeric(igraph::V(g)),igraph::vertex.attributes(g)[[m1]])
  m2.nodes <- split(as.numeric(igraph::V(g)),igraph::vertex.attributes(g)[[m2]])

  #get module sets names
  m1.names <- unique(igraph::vertex.attributes(g)[[m1]])
  m2.names <- unique(igraph::vertex.attributes(g)[[m2]])

  #get matched sets
  matches <- .match_modules(m1.nodes, m2.nodes, m1.names, m2.names)
  #count total matches and total nodes in modules
  matched.nodes <- 0
  n.nodes <- vcount(g)
  for (match in matches) {
    matched.nodes = matched.nodes + as.numeric(match[5])
  }
  return((matched.nodes / n.nodes) * 100)
}

.match_modules <- function(m1.nodes, m2.nodes, m1.names, m2.names){
  matches <- vector(mode = "list",
                    length = min(length(m1.nodes), length(m2.nodes)))
  n.matched <- 0
  matched1 <- rep(10^10, length(m1.nodes)) #wont have effect until replaced with actual matches
  matched2 <- rep(10^10, length(m2.nodes))

  while((length(m1.nodes) > n.matched) && (length(m2.nodes) > n.matched)){
    #find the best match
    max.overlap <- -1
    best <- list(NULL, NULL)
    for (i in (1:length(m1.nodes))[-as.numeric(matched1)]) { #iterate through both sets of nodes
      for (j in (1:length(m2.nodes))[-as.numeric(matched2)]) {
        overlap <- length(
          BiocGenerics::intersect(
            as.numeric(m1.nodes[[i]]),as.numeric(m2.nodes[[j]])
          ) #recording overlap
        )
        if (overlap > max.overlap) {
          max.overlap <- overlap
          best <- c(i, j, m1.names[i], m2.names[j], overlap) #index from set 1, set2, name set 1, set 2, overlap
        }
      }
    } #after finishing iteration, current best match is kept

    matches[[n.matched+1]] <- best
    matched1[n.matched+1] <- best[1] #removes from future consideration
    matched2[n.matched+1] <- best[2]
    n.matched <- n.matched + 1
  }
  return(matches)
}

## Assess module quality based on the number of nodes with more edges outside
## of module then within it
module_contiguity <- function(g, module.name = "module"){
  #get module assignments
  mods <-igraph::vertex.attributes(g)[[module.name]]
  adj <- as.matrix(as_adjacency_matrix(g))

  #get number nodes with more neighbors within module than between modules
  moreWithin <- length(which(unlist(
    lapply(1:nrow(adj), function(row){
      node.mod <- mods[row]
      edge.mod <- mods[which(adj[row,] == 1)]
      length(which(edge.mod == node.mod)) > length(which(edge.mod != node.mod))
    })
  )))
  return(round((moreWithin / length(g)) * 100, 4))
}


eigen_fuzzy_modules <- function(x, modules, max.size){
  lapply(modules, function(mod){
    #get number of required fuzzy nodes
    n.fuzzy.nodes <- max.size - length(mod)

    #get the modules eigen gene
    modPC <- prcomp(t(x[mod,]), scale. = TRUE)
    modEigen <- modPC$x[,1]

    #get genes outside the module the covary with the eigen gene
    eigenCov <- apply(x[!(rownames(x) %in% mod),], 1, function(x) cov(x, modEigen)) #covariance
    covRank <- sort(abs(eigenCov), decreasing = TRUE) #ranked absolute covariance
    fuzzy.nodes <- covRank[1:n.fuzzy.nodes] #n.fuzzy.nodes nodes with the highest ranks

    #return fuzzy module combining original and fuzzy nodes
    sort(c(mod, names(fuzzy.nodes)))
  })
}



############
### TEST ###
############

# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/00.SimulateGraphs.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/01.DetectModules.R")
#
# # make data
# er <- make_modular_graph()
# x <- sim_graph_data(er, n.samples = 10)
#
# # learn mods
# i <- find_ICA_mods(x, 3)
#
# # fuzzy mods
# f <- eigen_fuzzy_modules(x, split(V(er)$name, i), 80)
#
# #assess accuracy
# igraph::V(er)$ica <- i
# percent_module_match(er, "module", "ica")
# module_contiguity(er, "ica")
