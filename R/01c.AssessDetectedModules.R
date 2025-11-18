#' Calculates the percentage of nodes assigned to the same module between to module sets, based on node index, not module name
#' @param m1 a module object
#' @param m2 another module object to compare to

#' @return an integer between 0 and 100

#' @export
percent_module_match <- function(m1, m2){
  #get matched sets
  matches <- .match_modules(m1@index.list, m2@index.list, m1@name.list, m2@name.list)
  #count total matches and total nodes in modules
  matched.nodes <- 0
  n.nodes <- min(c(length(unlist(m1@index.vector)), length(m2@index.vector)))
  for (match in matches) {
    matched.nodes = matched.nodes + as.numeric(match[5])
  }
  return((matched.nodes / n.nodes) * 100)
}

#' Finds the best module matches between two sets of nodes
#' @param m1.nodes Node index vector from m1
#' @param m2.nodes Node index vector from m2
#' @param m1.names Node character vector from m1
#' @param m2.names Node character vector from m2

#' @return a list of all the best modules matched based on node membership

#' @keywords internal
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
          intersect(
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

#' Calculate the percentage of nodes with more edges within their module then outside of it
#' @param g an igraph object
#' @param test.module a module to test for contiguity, whose nodes indexs and names must match the graph

#' @return an integer between 0 and 100

#' @importFrom igraph as_adjacency_matrix

#' @export
module_contiguity <- function(g, test.module){
  #get module assignments
  index.vector <- test.module@index.vector
  adj <- as.matrix(igraph::as_adjacency_matrix(g))

  #get number nodes with more neighbors within module than between modules
  more.within <- length(which(unlist(
    lapply(1:nrow(adj), function(row){
      node.mod <- index.vector[row]
      edge.mod <- index.vector[which(adj[row,] == 1)]
      length(which(edge.mod == node.mod)) > length(which(edge.mod != node.mod))
    })
  )))
  return(round((more.within / length(g)) * 100, 4))
}


############
### TEST ###
############


# # make data
# er <- make_modular_graph()
# x <- sim_graph_data(er, n.samples = 10)
#
# # learn mods
# true_modules <- methods::new("module",
#                                         source = "True Modules",
#                                         overlapping = FALSE,
#                                         index.vector = igraph::V(er)$module,
#                                         index.list = split(1:120 , igraph::V(er)$module),
#                                         name.list = split(igraph::V(er)$name , igraph::V(er)$module)
#                     )
# i <- find_ICA_mods(x, 3)
#
#
# test1 <- methods::new("module",
#                     source = "test",
#                     overlapping = FALSE,
#                     index.vector = c(1,1,2,2,3,3),
#                     index.list = list(c(1,2),c(3,4),c(5,6)),
#                     name.list = list(c("Node_1","Node_2"),c("Node_3","Node_4"),c("Node_5","Node_6"))
# )
# test2 <- methods::new("module",
#              source = "test",
#              overlapping = FALSE,
#              index.vector = c(1,1,2,3,2,3),
#              index.list = list(c(1,2),c(3,5),c(4,6)),
#              name.list = list(c("Node_1","Node_2"),c("Node_3","Node_5"),c("Node_4","Node_6"))
# )
#
#
# #assess accuracy
# percent_module_match(test1, test2)
# module_contiguity(er, true_modules)
