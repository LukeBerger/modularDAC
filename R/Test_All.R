#libs
library(dplyr)
library(WGCNA)
library(parallel)
library(foreach)
library(doParallel)

#load functions
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/00.SimulateGraphs.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01.DetectModules.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01b.complexICA.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01c.AssessDetectedModules.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/02.LearnGraphs.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/03.Multithreaded_DaC.R")


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


  ## 01.DetectModules (NEED TO UPDATE WITH MODULES OBJECT)
  # learn mods
  w <- suppressMessages(find_WGCNA_mods(t(x), cor.FN = "bicor"))
  i <- find_ICA_mods(x, 3)
  p <- pragmatic_modules(x,n.mods = 3, max.size = 60)

  # fuzzy mods
  f <- eigen_fuzzy_modules(x, lapply(p, function(x) x$index), 80)
  f <- nodewise_fuzzy_modules(x, lapply(p, function(x) x$index), 80)

  # overlap mods
  o <- create_overlap_modules(x, lapply(p, function(x) x$index))

  print("01.DetectModules : All functions passed all tests!")


  ## 01b.complexICA
  # complex ICA mods
  c <- complexICA(x, 3)

  print("01b.complexICA : All functions passed all tests!")


  ## 01c.AssessDetectedModules
  # add modules to er igraph object
  igraph::V(er)$ica <- i

  # percent match
  pmm <- percent_module_match(er, "module", "ica")
  if(!is.numeric(pmm) || !dplyr::between(pmm, 0, 100)){stop("percent_module_match failed to return a numeric between 0 and 100")}

  # module contiguity
  mc <- module_contiguity(er, "ica")
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


  ## 03.MultiTrheaded_DaC
  # Test Function with Defaults
  dac <- fastDivideAndConquer(x, f) # basic run
  if(!igraph::is_igraph(dac$final.graph)){stop("fastDivideAndConquer failed to produce an igraph with default functions")}
  if(!all(unlist(lapply(dac$modular.subgraphs, igraph::is_igraph)))){stop("fastDivideAndConquer failed to produce a list of modular igraphs with default functions")}

  #test function with BDgrapg learning
  dac <- fastDivideAndConquer(x, f,
                               graph.learning.func = bdWrapper,
                               arg.wrapping.func = .other_arg_wrapper,
                               out.parsing.func = .default_output_parser,
                               packages.to.each = c("igraph", "BDgraph"),
                               export.to.each = c("bdWrapper")
  )
  if(!igraph::is_igraph(dac$final.graph)){stop("fastDivideAndConquer failed to produce an igraph with bdgraph learning")}
  if(!all(unlist(lapply(dac$modular.subgraphs, igraph::is_igraph)))){stop("fastDivideAndConquer failed to produce a list of modular igraphs bdgraph learning")}

  print("03.MultiTrheaded_DaC : All functions passed all tests!")


  print("Package testing complete! modularDAC is fully functional!")
}

### Alternative functions for testing fastDivideAndConquer
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



#create and register the cluster
cl <- parallel::makeCluster(detectCores() -1)
doParallel::registerDoParallel(cl)

#run tester
test_all()

# Stop the cluster
stopCluster(cl)
