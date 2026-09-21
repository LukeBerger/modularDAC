#' Graph learning wrapper for RSNet: generates ensemble and consensus GGM networks from a module data subset
#' @param dat a numeric matrix with p features (rows) and n samples (columns), as prepared by .RSNet_arg_wrapper
#' @param num_iteration an integer, the number of bootstrap/subsampling iterations to perform
#' @param boot a logical, if TRUE bootstrap resampling is used
#' @param sub_ratio a numeric between 0 and 1, the subsampling ratio applied to samples
#' @param sample_class a character vector or NULL, class labels for samples to include; NULL uses all samples
#' @param correlated a logical, if TRUE correlated feature structure is accounted for
#' @param cluster_ratio a numeric, the clustering ratio parameter
#' @param estimate_CI a logical, if TRUE confidence intervals are estimated for each edge
#' @param method a character, the statistical inference method passed to RSNet (e.g. 'B_NW_SL')
#' @param n_cores an integer or NULL, the number of cores for parallel processing; NULL uses a single core
#' @param CI a numeric between 0 and 1, the confidence level for interval estimation
#' @param filter a character, the filtering criterion applied to edges (e.g. 'pval')
#' @param threshold a numeric, the threshold value for the chosen filter


#' @return a list of learned graphs (and other consensus_net_ggm outputs)

#' @keywords internal
.RSNet_full_learning <- function(dat,
                                 num_iteration = 1,
                                 boot = FALSE,
                                 sub_ratio = 1,
                                 sample_class = NULL,
                                 correlated = FALSE,
                                 cluster_ratio = 1,
                                 estimate_CI = TRUE,
                                 method = "B_NW_SL",
                                 n_cores = NULL,
                                 CI = 0.95,
                                 filter = "pval",
                                 threshold = 0.05){
  if (!requireNamespace("RSNet", quietly = TRUE)) {
    stop("Package RSNet is required. Install with: install.packages('RSNet')", call. = FALSE)
  }

  # RSNet expects an n samples x p features matrix
  t.dat <- t(dat)

  ens.net <- RSNet::capture_all(RSNet::ensemble_ggm(t.dat,
                                 num_iteration,
                                 boot,
                                 sub_ratio,
                                 sample_class,
                                 correlated,
                                 cluster_ratio,
                                 estimate_CI,
                                 method,
                                 n_cores = n_cores
  ))
  conc.nets <-  RSNet::capture_all(RSNet::consensus_net_ggm(ggm_networks = ens.net,
                                        CI = CI,
                                        filter = filter,
                                        threshold = threshold

  ))
  return(conc.nets)
}

#' Argument wrapper for RSNet: packages each module data subset with RSNet arguments
#' @param sub.x a list of numeric matrices, each being a p_i x n subset of x for one module (features as rows)
#' @param num_iteration an integer, the number of bootstrap/subsampling iterations to perform
#' @param boot a logical, if TRUE bootstrap resampling is used
#' @param sub_ratio a numeric between 0 and 1, the subsampling ratio applied to samples
#' @param sample_class a character vector or NULL, class labels for samples to include; NULL uses all samples
#' @param correlated a logical, if TRUE correlated feature structure is accounted for
#' @param cluster_ratio a numeric, the clustering ratio parameter
#' @param estimate_CI a logical, if TRUE confidence intervals are estimated for each edge
#' @param method a character, the statistical inference method passed to RSNet (e.g. 'B_NW_SL')
#' @param n_cores an integer or NULL, the number of cores for parallel processing; NULL uses a single core
#' @param CI a numeric between 0 and 1, the confidence level for interval estimation
#' @param filter a character, the filtering criterion applied to edges (e.g. 'pval')
#' @param threshold a numeric, the threshold value for the chosen filter

#' @return a list containing the args for rsnet learning

#' @keywords internal
.RSNet_arg_wrapper <- function(sub.x,
                               num_iteration  = 1,
                               boot = FALSE,
                               sub_ratio = 1,
                               sample_class = NULL,
                               correlated = FALSE,
                               cluster_ratio = 1,
                               estimate_CI = TRUE,
                               method = "B_NW_SL",
                               n_cores = NULL,
                               CI = 0.95,
                               filter = "pval",
                               threshold = 0.005){

  # return a list of packaged args (each module stays feature x sample;
  # .RSNet_full_learning transposes internally)
  lapply(sub.x, function(x){
    list(dat = x,
         num_iteration = num_iteration,
         boot = boot,
         sub_ratio = sub_ratio,
         sample_class = sample_class,
         correlated = correlated,
         cluster_ratio = cluster_ratio,
         estimate_CI = estimate_CI,
         method = method,
         n_cores = n_cores,
         CI = CI,
         filter = filter,
         threshold = threshold)
  })
}

#' Output parser for RSNet: extracts consensus networks and other outputs from RSNet graph learning results
#' @param graph.learning.outputs a list of outputs, one per sub-graph, returned by .RSNet_full_learning

#' @return a list containing the learned subgraphs and any other outputs of graph learning

#' @keywords internal
.RSNet_output_parser <- function(graph.learning.outputs){
  nets <- lapply(graph.learning.outputs, function(o) o$consensus_network)
  other <- lapply(graph.learning.outputs, function(o) o[-1])
  return(
    list(
      learned.graphs = nets,
      other.outputs = other
    )
  )
}
