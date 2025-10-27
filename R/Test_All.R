# #libs
# library(dplyr)
# library(WGCNA)
# library(parallel)
# library(foreach)
# library(doParallel)
#
# #load functions
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/00.SimulateGraphs.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01.DetectModules.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01b.complexICA.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01c.AssessDetectedModules.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/02.LearnGraphs.R")
# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/03.Multithreaded_DaC.R")


### Sources and runs tests for all functions in package
test_all <- function(){

  ## 00.SimulateGraphs
  # make graphs
  er <- make_modular_graph()
  hybrid <- make_submodular_graph()
  if(!igraph::is_igraph(er)){stop("make_module_graph failed to produce an igraph")}
  if(!igraph::is_igraph(hybrid)){stop("make_submodular_graph failed to produce an igraph")}

  # plot
  vn <- modular_plot(er)
  if(!all(c("visNetwork", "htmlwidget") %in% class(vn))){stop("modular_plot failed to produce a plot")}

  # sim data
  x <- sim_graph_data(er, n.samples = 100)
  if(!is.matrix(x) || !is.numeric(x)){
    stop("sim_graph_data failed to produce a numeric matrix")
  }else{
    if(nrow(x) != igraph::vcount(er) || ncol(x) != 100){
      stop("sim_graph_data failed to produce a n x p matrix")
    }
  }

  print("00.SimulateGraphs : All functions passed all tests!")


  ## 01.DetectModules
  # learn mods
  w <- suppressMessages(find_WGCNA_mods(t(x), cor.FN = "bicor"))
  i <- find_ICA_mods(x, 3)
  p <- pragmatic_modules(x,n.mods = 3, max.size = 60)

  # fuzzy mods
  ef <- eigen_fuzzy_modules(x, p, 80)
  nf <- nodewise_fuzzy_modules(x, p, 80)

  # overlap mods
  o <- create_overlap_modules(x, p)

  ## 01b.complexICA
  # complex ICA mods
  c <- complex_ICA_modules(x, 3)

  #test they are proper modules
  for(mod in list(w, i ,p, ef, nf, o)){
    if(mod@overlapping){
      #check that overlaps exist between all modules
      if(!all(
        unlist(
          lapply(seq_along(mod@index.list), function(i){
            any(unlist(mod@index.list[i]) %in% unlist(mod@index.list[-i]))
          })
        )
      )){stop(paste(mod@source, "has modules with incomplete overlaps"))}
    }else{
      #check feature number matches input data
      if(length(mod@index.vector) != nrow(x)){
        stop(paste(mod@source, "produced a index vector with the incorrect number of features"))
      }
      #check that there are no overlaps
      if(any(
        unlist(
          lapply(seq_along(mod@index.list), function(i){
            any(unlist(mod@index.list[i]) %in% unlist(mod@index.list[-i]))
          })
        )
      )){stop(paste(mod@source, "has modules with overlaps"))}
    }

    #check that all feature names in modules come from the data
    if(!all(unlist(mod@name.list) %in% rownames(x))){
      stop(paste(mod@source, "feature names do not match input data"))
    }
    #check that each module has the same number of feature indexes and names
    if(
      !all(
        unlist(
          lapply(seq_along(mod@index.list), function(i) length(mod@index.list[[i]]) == length(mod@name.list[[i]]))
          )
        )
    ){stop(paste(mod@source, "produced differnet length index and name lists"))}

  }

  print("01.DetectModules : All functions passed all tests!")


  ## 01c.AssessDetectedModules
  # percent match
  pmm <- percent_module_match(w, i)
  if(!is.numeric(pmm) || !dplyr::between(pmm, 0, 100)){stop("percent_module_match failed to return a numeric between 0 and 100")}

  # module contiguity
  mc <- module_contiguity(er, i)
  if(!is.numeric(mc) || !dplyr::between(mc, 0, 100)){stop("module_contiguity failed to return a numeric between 0 and 100")}

  print("01c.AssessDetectedModule : All functions passed all tests!")


  ## 02.LearnGraphs
  # learn graph with silgmm
  s <- learn_SILGGM_graph(t(x))
  if(!igraph::is_igraph(s)){stop("learn_SILGGM_graph failed to produce an igraph")}
  if(!all(igraph::V(s)$name == row.names(x))){stop("learn_SILGGM_graph produced a graph with inproper node names")}

  # get f1 score
  f1 <-calc_F1(er, s)
  if(!is.numeric(f1$F1)){
    stop("calc_f1 did not produce a numeric F1 score")
  }else{
    if(!dplyr::between(f1$F1, 0, 1)){
      stop("calc_f1 did not produce a F1 score between 0 and 1")
    }
  }

  print("02.LearnGraphs : All functions passed all tests!")


  ## 03.MultiThreaded_DaC
  # Test Function with Defaults
  dac <- divide_and_conquer(x, ef@index.list) # basic run
  if(!igraph::is_igraph(dac$final.graph)){stop("divide_and_conquer failed to produce an igraph with default functions")}
  if(!all(unlist(lapply(dac$modular.subgraphs, igraph::is_igraph)))){stop("divide_and_conquer failed to produce a list of modular igraphs with default functions")}

  #test function with BDgrapg learning
  dac <- divide_and_conquer(x, ef@index.list,
                               graph.learning.func = bdWrapper,
                               arg.wrapping.func = .other_arg_wrapper,
                               out.parsing.func = .default_output_parser,
                               packages.to.each = c("igraph", "BDgraph"),
                               export.to.each = c("bdWrapper")
  )
  if(!igraph::is_igraph(dac$final.graph)){stop("divide_and_conquer failed to produce an igraph with bdgraph learning")}
  if(!all(unlist(lapply(dac$modular.subgraphs, igraph::is_igraph)))){stop("divide_and_conquer failed to produce a list of modular igraphs bdgraph learning")}

  print("03.MultiTrheaded_DaC : All functions passed all tests!")


  print("Package testing complete! modularDAC is fully functional!")
}

### Alternative functions for testing divide_and_conquer
.other_arg_wrapper <- function(sub.x, n.bootstraps, sub.ratio){
  lapply(sub.x, function(x){
    list(
      x = t(x)
    )
  })
}

bdWrapper <- function(x){
  bdFit <- BDgraph::bdgraph(x)
  postProb <- BDgraph::plinks(bdFit)
  adjMat <- BDgraph::select(postProb)
  igraph::graph_from_adjacency_matrix(adjMat, mode = "undirected") # can change this if I want to test edge directional update for .connect_subgraphs
}



# #create and register the cluster
# cl <- parallel::makeCluster(detectCores() -1)
# doParallel::registerDoParallel(cl)
#
# #run tester
# test_all()
#
# # Stop the cluster
# stopCluster(cl)
