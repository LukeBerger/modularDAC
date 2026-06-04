
######################
### Define Modules ###
######################

# declare module object
setClass("module",
         slots = list(
           source = "character",
           data.dim = "numeric",
           overlapping = "logical",
           index.vector = "numeric",
           score.vector = "numeric",
           index.list = "list",
           name.list = "list"
         ),
         validity = function(object){
           # check that each module has the same number of feature indexes and names
           if(
             !all(
               unlist(
                 lapply(seq_along(object@index.list), function(i){
                   length(object@index.list[[i]]) == length(object@name.list[[i]])
                 })
               )
             )
           ){return("index list and name list must be the same length")}
           # if not overlapping...
           if(!object@overlapping){
             # check that every module contains a unique index
             if(any(unlist(
               lapply(seq_along(object@index.list), function(i){
                 any(unlist(object@index.list[i]) %in% unlist(object@index.list[-i]))
               })
             ))){return("All modules must contain overlaps with at least on other module")}
           }

         }

)

#' Extract true module assignments from a simulated graph
#' @param g an igraph object with a 'module' vertex attribute defining ground-truth module membership for each node

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V

#' @export
true_modules <- function(g){
  # build module object from ground-truth module vertex attribute
  methods::new("module",
                source = "True Modules",
                overlapping = FALSE,
                index.vector = igraph::V(g)$module,
                index.list = split(1:length(g) , igraph::V(g)$module),
                name.list = split(igraph::V(g)$name , igraph::V(g)$module)
  )
}

#' Expand true modules to fuzzy (overlapping) modules based on second-order graph neighborhood
#' @param m a module S4 object containing module membership assignments
#' @param g an igraph object with a 'module' vertex attribute defining ground-truth module membership for each node

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V neighborhood

#' @export
true_fuzzy <- function(m, g){
  # collect all nodes within 2 graph hops of each module
  f.index.list <- lapply(m@index.list, function(x){
    sort(unique(unlist(
      igraph::neighborhood(
        g,
        order = 2,
        nodes = igraph::V(g)[x],
        mode = "all",
        mindist = 0
      )
    )))
  })
  # build fuzzy module object with expanded index list
  methods::new("module",
               source = "True Module Fuzzy",
               overlapping = TRUE,
               index.list = f.index.list,
               name.list = lapply(f.index.list, function(m) igraph::V(g)$name[m])
  )
}

#' Validate a module object against the data matrix and throw errors if checks fail
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param m a module S4 object containing module membership assignments

#' @return TRUE, if all checks passed

#' @keywords internal
.module_check <- function(x, m){
  if(m@overlapping){
    # check that overlaps exist between all modules
    if(!all(
      unlist(
        lapply(seq_along(m@index.list), function(i){
          any(unlist(m@index.list[i]) %in% unlist(m@index.list[-i]))
        })
      )
    )){stop(paste(m@source, "has mules with incomplete overlaps"))}
  }else{
    # check feature number matches input data
    if(length(m@index.vector) != nrow(x)){
      stop(paste(m@source, "produced a index vector with the incorrect number of features"))
    }
    # check that there are no overlaps
    if(any(
      unlist(
        lapply(seq_along(m@index.list), function(i){
          any(unlist(m@index.list[i]) %in% unlist(m@index.list[-i]))
        })
      )
    )){stop(paste(m@source, "has mules with overlaps"))}
  }

  # check that all feature names in modules come from the data
  if(!all(unlist(m@name.list) %in% rownames(x))){
    stop(paste(m@source, "feature names do not match input data"))
  }
  # check that each module has the same number of feature indexes and names
  if(
    !all(
      unlist(
        lapply(seq_along(m@index.list), function(i) length(m@index.list[[i]]) == length(m@name.list[[i]]))
      )
    )
  ){stop(paste(m@source, "produced differnet length index and name lists"))}
  return(TRUE)
}

##################################
### Module Detection Functions ###
##################################

#' Detect co-expression modules from a data matrix using WGCNA
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param min.size an integer, the minimum number of nodes allowed in a module (passed to cutreeDynamic)
#' @param max.size an integer, the maximum number of nodes allowed in a module
#' @param min.sft a numeric between 0 and 1, the minimum R-squared threshold for soft-thresholding power selection
#' @param beta an integer, the soft-thresholding power for WGCNA::adjacency; if NULL it is selected automatically via pickSoftThreshold
#' @param cor.FN a character, the correlation function to use in WGCNA::adjacency; either 'bicor' or 'cor'
#' @param powers an integer vector, candidate soft-thresholding powers evaluated by WGCNA::pickSoftThreshold
#' @param hclust.method a character, the agglomeration method passed to flashClust
#' @param cut.height a numeric between 0 and 1, the dendrogram cut height used by cutreeDynamic
#' @param merge a logical, if TRUE adjacent modules are merged using WGCNA::mergeCloseModules
#' @param merging.cut a numeric between 0 and 1, the eigengene dissimilarity threshold for WGCNA::mergeCloseModules
#' @param iterate a logical, if TRUE module detection is run iteratively to break up any modules exceeding max.size

#' @return a module object

#' @importFrom stats as.dist
#' @importFrom methods new

#' @export
find_WGCNA_mods <- function(x,
                            min.size=10,
                            max.size=NULL,
                            min.sft=0.85,
                            beta=NULL,
                            cor.FN=c("bicor", "cor"),
                            powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                            hclust.method="average",
                            cut.height = NULL,
                            merge = F,
                            merging.cut = 0.2,
                            iterate = T
) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package WGCNA is required. Install with: install.packages('WGCNA')", call. = FALSE)
  }
  if (!requireNamespace("flashClust", quietly = TRUE)) {
    stop("Package flashClust is required. Install with: install.packages('flashClust')", call. = FALSE)
  }
  if (!requireNamespace("dynamicTreeCut", quietly = TRUE)) {
    stop("Package dynamicTreeCut is required. Install with: install.packages('dynamicTreeCut')", call. = FALSE)
  }
  if (!requireNamespace("MatrixGenerics", quietly = TRUE)) {
    stop("Package MatrixGenerics is required. Install with: install.packages('MatrixGenerics')", call. = FALSE)
  }
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Package matrixStats is required. Install with: install.packages('matrixStats')", call. = FALSE)
  }

  # handle arguments
  cor.FN <- match.arg(cor.FN)

  # merge / iterate relationship
  if(merge & iterate){
    cat("Note: Using both the 'merge' and 'iterate' options, may result in some modules being merged then immediately split, use with caution. \n ")
  }

  # correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # check that all rows have variance
  rv <- MatrixGenerics::rowVars(x)
  if(any(rv == 0)){
    warning("Warning: some matrix elements have zero variance. They will be removed from consideration as they would not be placed in any module")
    rm.r <- which(rv == 0)
    cat("Removed rows: ", rm.r  )
    x <- x[-rm.r, , drop = FALSE]
  }

  # pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=t(x),
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # check selected power
    beta <- .sft_check(sft)
  }

  # if max size is constrained
  if(!is.null(max.size)){
    # determine the minimum number of modules needed to fit all nodes
    min.mods <- ceiling(nrow(x) / max.size)
  }else{ min.mods <- 1}

  # construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=t(x),
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # topological overlap distance transformation
  dis <- WGCNA::TOMdist(adjMat=adj, TOMType="unsigned")

  # fast hierarchical clustering of distance
  dendro <- flashClust::flashClust(d=stats::as.dist(dis), method=hclust.method)

  # module identification using dynamic tree cut algorithm
  cat("Generating initial modules... \n")
  initial.index.vector <- dynamicTreeCut::cutreeDynamic(dendro=dendro,
                                           cutHeight = cut.height,
                                           method="hybrid",
                                           distM=dis,
                                           deepSplit=4,
                                           pamRespectsDendro=FALSE,
                                           minClusterSize=min.size)

  # if this produced less than the minimum number of mods
  n.mods <- length(unique(initial.index.vector[initial.index.vector != 0]))
  if( n.mods < min.mods ){
    cat("Only ", n.mods, " modules produced in first cut. \n",
        min.mods, " will be required with max.size == ", max.size, "\n",
        "Lowering cut height... \n",
        sep = "")

    # if cut height is null
    if(is.null(cut.height)){
      # recalculate the default height used by cuttreeDynamic
      qnt5 <-  quantile(dendro$heigh, probs= 0.05)
      cut.height <- qnt5 + 0.99 * (max(dendro$height) - qnt5)
    }

    # lower height until the desired number of modules is met
    while(n.mods < min.mods){
      cut.height <- cut.height - 0.01
      cat("Lowerin cut height to", cut.height, "and generating new modules... \n")
      initial.index.vector <- dynamicTreeCut::cutreeDynamic(dendro=dendro,
                                                        cutHeight = cut.height,
                                                        method="hybrid",
                                                        distM=dis,
                                                        deepSplit=4,
                                                        pamRespectsDendro=FALSE,
                                                        minClusterSize=min.size)
      # if producing all unassigned
      if(sum(initial.index.vector == 0) == length(initial.index.vector) && cut.height < 0.8){
        cat("Cut height reached " , cut.height, " and is no longer to producing meaningful clusters. \n",
            "Consider increasing max.size or using another method to generate modules. \n",
            sep = "")
        stop("No modules found after lowering cut height")
      }
      n.mods <- length(unique(initial.index.vector[initial.index.vector != 0]))
    }
    cat("Final cut height" , cut.height, "produced", n.mods, "modules. \n")
  }else{
      cat("Initial WGCNA call produced", n.mods, "modules \n")
  }

  # if merging
  if(merge){
    cat("Merging modules based on eigen gene similarity... \n")
    n.merges <- 0
    # while there are at least min.mods modules and they are all bellow maxsize
    while(n.mods > min.mods && all(table(initial.index.vector  ) < max.size ) ){
      # merge similar modules (single iteration)
      merged <- WGCNA::mergeCloseModules(
        exprData=t(x),
        colors=initial.index.vector,
        # mEs= ifelse(exists("merged"), merged$newMEs, NULL),
        unassdColor=0,
        corFnc=cor.FN,
        corOptions=cor.options,
        cutHeight=merging.cut
      )

      # if nothing changes, break loop
      if(all(initial.index.vector == merged$colors)){ # this section of code is outdate, need to revise in later
        break
      }
      # else, use merged modules and restart
      initial.index.vector <- merged$colors
      n.mods <- length(unique(initial.index.vector[initial.index.vector != 0]))
      n.merges <- n.merges + 1
    }
    if(n.merges == 0){
      cat("No modules were small/similiar enough to merge, using initial modules... \n")
    }else{
      cat("Merged" , n.merges, "modules, resulting in", n.mods, "new modules. \n")
    }

  }

  # if iterating
  if(iterate){
    # find modules that are too large
    to.big <- as.numeric(names(which(table(initial.index.vector) > max.size)))
    to.big <- to.big[to.big!=0]

    cat("The first iterations of modules",to.big, "are too large, attempting to split... \n")

    # rerun flast cluster for each module that is to large
    for(m in to.big){
      # subset distance matrix
      m.nodes <- which(initial.index.vector == m)
      m.dis <- dis[m.nodes,m.nodes, drop = F]

      # regenerate modules
      dendro <- flashClust::flashClust(
        d=stats::as.dist(m.dis),
        method=hclust.method)
      m.index.vector <- dynamicTreeCut::cutreeDynamic(
        dendro=dendro,
        cutHeight = cut.height,
        method="hybrid",
        distM=m.dis,
        deepSplit=4,
        pamRespectsDendro=FALSE,
        minClusterSize=min.size,
        verbose = F)

      # check new modules...
      pass = T
      # arent mostly zero
      if(sum(m.index.vector ==0) > length(m.index.vector) / 2){
        pass = F
        cat("Module", m, "failed to split due to high number of unassigned nodes. \n")
      }
      # have at least two non zero values
      if(length(unique(m.index.vector[m.index.vector != 0])) < 2){
        pass = F
        cat("Module", m,"failed to spilt into at least two new modules. \n")
      }

      # if the module passes checks,
      if(pass){
        # split module
        cat("Spliting module" ,m, "into", length(unique(m.index.vector[m.index.vector!=0])), "new modules. \n")
        initial.index.vector[m.nodes] <- ifelse(m.index.vector  ==0,
                                                0,  m.index.vector  + max(initial.index.vector)) # adding max ensures no overlaps)
      }
    }
  }

  # assign remaining unassigned nodes based on adj...
  cat("Assinging unassinged nodes based on WGNCA adjacency... \n")
  diag(adj) <- 0

  # iterative expansion by best neighbor with threshold
  modified.index.vector <- initial.index.vector
  last.mods <- modified.index.vector
  unassigned <- which(modified.index.vector == 0)
  threshold = 0.95
  while(length(unassigned) > 0){ # there might be a more efficient way to check this
    # for each assigned node
    assigned <- which(modified.index.vector != 0)
    unassigned <- which(modified.index.vector == 0)
    for(node in assigned){
      # find its best neighbor based on wgcna adj
      best <- unassigned[which.max(adj[node, unassigned])]
      # if that adj is above the current threshold
      if(length(best) !=0){
        if(adj[node,best] > threshold){
          # assign node
          modified.index.vector[best] <- modified.index.vector[node]
        }
      }
    }
    # if no changes
    if(all(last.mods == modified.index.vector)){
      # and not at min threshold
      if(threshold >= 0){
        # lower threshold
        threshold <- round(threshold - 0.05, 2)
      }else{
        # break if at min threshold and no longer changing
        break
      }
    }
    last.mods <- modified.index.vector
  }

  # fit modules to max size via trading
  to.big  <- sum(table(modified.index.vector) > max.size)
  if(to.big){
    cat(to.big, "modules are larged than", max.size, "nodes. Preforming node trading based on adjacency... \n")
    modified.index.vector <- .fit_to_max_size(
      index.vector = modified.index.vector,
      adj = adj,
      max.size = max.size)
  }

  # create module object
  initial.index.vector <- as.numeric(initial.index.vector)
  modified.index.vector <- as.numeric(modified.index.vector)
  initial.mods  <- methods::new("module",
                     source = "find_WGCNA_mods initial mods",
                     data.dim = dim(x),
                     overlapping = FALSE,
                     index.vector = initial.index.vector,
                     index.list = split(1:nrow(x) , initial.index.vector),
                     name.list = split(rownames(x), initial.index.vector)
  )
  final.mods  <- methods::new("module",
                                source = "find_WGCNA_mods traded mods",
                                data.dim = dim(x),
                                overlapping = FALSE,
                                index.vector = modified.index.vector,
                                index.list = split(1:nrow(x) , modified.index.vector),
                                name.list = split(rownames(x), modified.index.vector)
  )

  return(list(
    wgcna.adj = adj,
    initial.mods = initial.mods,
    final.mods = final.mods
  ))
}

#' Select soft-thresholding power from WGCNA::pickSoftThreshold output
#' @param sft a list, the output of WGCNA::pickSoftThreshold

#' @return an integer, the selected soft-thresholding power

#' @keywords internal
.sft_check <- function(sft) {
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- 6 # default
    cat("Using the following power:", beta, "\n")
  } else {
    cat("Optimal power selected:", beta, "\n")
  }
  return(beta)
}

#' Helper to find_WGCNA_mods that subdivides modules exceeding a maximum size
#' @param index.vector an integer vector of length p, assigning each node to a module (0 = unassigned)
#' @param adj a p x p numeric adjacency matrix produced by WGCNA::adjacency
#' @param max.size an integer, the maximum number of nodes allowed in a module
#' @param verbose a logical, if TRUE progress messages are printed for each split (useful for debugging)

#' @return an integer vector of length p with the updated module assignments

#' @keywords internal

.fit_to_max_size <- function(
    index.vector,
    adj,
    max.size,
    verbose = F

){
  # determine which modules will give and receive nodes
  mod.size <- table(index.vector)

  giving <- as.numeric(names(which(mod.size > max.size)))
  receiving <- as.numeric(names(which(mod.size < max.size)))

  giving.nodes <- which(index.vector %in% giving)
  receiving.nodes <- which(index.vector %in% receiving)

  n.trades <- 0

  while(length(giving) > 0){
    # find the giving node with the highest adj to a  receiving module
    trade.scores <- adj[giving.nodes, receiving.nodes, drop = FALSE]
    trade.giving <- giving.nodes[ which.max(matrixStats::rowMaxs(trade.scores)) ]
    trade.receiving  <-  receiving.nodes[
      which.max( adj[trade.giving,receiving.nodes] )
      ]

    # set giving nodes module to receiving module
    if(verbose){cat("Trading", trade.giving, "to module", index.vector[ trade.receiving ], "\n" )}
    index.vector[ trade.giving ] <- index.vector[ trade.receiving ]

    # update module status
    mod.size <- table(index.vector)

    giving <- as.numeric(names(which(mod.size > max.size)))
    receiving <- as.numeric(names(which(mod.size < max.size)))

    giving.nodes <- which(index.vector %in% giving)
    receiving.nodes <- which(index.vector %in% receiving)

    n.trades <- n.trades + 1
  }
  cat("All modules fit within max size after",  n.trades, "trades. \n")

  return(index.vector)
}

#' Detect co-expression modules from a data matrix using Independent Component Analysis (ICA)
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param n.comp an integer, the number of independent components (modules) to extract
#' @param ... additional arguments passed to fastICA::fastICA

#' @return a module object

#' @importFrom methods new

#' @export
find_ICA_mods <- function(x,
                          n.comp,
                          ...){
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package fastICA is required. Install with: install.packages('fastICA')", call. = FALSE)
  }
  # run ica to decompose features into independent components
  ICA.results <- fastICA::fastICA(X= as.matrix(x),
                                  n.comp = n.comp,
                                  ...
  )
  # assign each feature to the component with the highest absolute loading
  index.vector <- apply(abs(ICA.results$S), 1, which.max)

  # convert to module object and return
  ICA.mods <- methods::new("module",
                     source = "find_ICA_mods",
                     data.dim = dim(x),
                     overlapping = FALSE,
                     index.vector = index.vector,
                     score.vector = apply(abs(ICA.results$S), 1, max),
                     index.list =split(1:nrow(x) , index.vector),
                     name.list = split(rownames(x), index.vector))

  return(ICA.mods)
}

#' Detect co-expression modules from a data matrix using Markov Clustering (MCL)
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param beta an integer, the soft-thresholding power for WGCNA::adjacency
#' @param cut a numeric between 0 and 1, the adjacency threshold used to binarize the weighted adjacency matrix
#' @param expansion an integer, the expansion parameter of the MCL algorithm
#' @param inflation an integer, the inflation parameter of the MCL algorithm
#' @param cor.fn a character, the correlation function to use in WGCNA::adjacency; either 'bicor' or 'cor'
#' @param iter an integer, the maximum number of iterations for the MCL algorithm

#' @return a module object

#' @importFrom methods new

#' @export
find_mcl_mods <- function(x,
                          beta=6,
                          cut=0.05,
                          expansion=2,
                          inflation=1,
                          cor.fn="cor",
                          iter=1000){
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package WGCNA is required. Install with: install.packages('WGCNA')", call. = FALSE)
  }
  if (!requireNamespace("MCL", quietly = TRUE)) {
    stop("Package MCL is required. Install with: install.packages('MCL')", call. = FALSE)
  }
  # compute unsigned and scaled correlation matrix
  mat <- WGCNA::adjacency(datExpr=x,
                          power=beta,
                          corFnc=cor.fn,
                          type="unsigned")

  # convert to adjacency matrix with cut point
  adj <- mat
  adj[mat >= cut] <- 1
  adj[mat < cut] <- 0
  diag(adj) <- 0

  # markov clustering
  clusters <- MCL::mcl(x=as.matrix(adj),
                  addLoops=TRUE,
                  allow1=FALSE,
                  ESM=FALSE,
                  max.iter=iter,
                  expansion=expansion,
                  inflation=inflation)

  # create module object
  mcl.mods <-methods::new("module",
                         source = "mcl",
                         data.dim = dim(t(x)),
                         overlapping = FALSE,
                         index.vector = clusters$Cluster,
                         index.list =split(1:ncol(x) , clusters$Cluster),
                         name.list = split(colnames(x), clusters$Cluster))

  # return Modules
  return(mcl.mods)
}

#' Detect co-expression modules from a data matrix using MEGENA
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param method a character, the correlation method passed to MEGENA::calculate.correlation
#' @param fdr.cutoff a numeric, the FDR cutoff used in MEGENA::calculate.correlation
#' @param mod.pval a numeric, the module significance p-value threshold for MEGENA::do.MEGENA
#' @param hub.pval a numeric, the hub significance p-value threshold for MEGENA::do.MEGENA
#' @param cor.perm an integer, the number of permutations used in MEGENA::calculate.correlation
#' @param hub.perm an integer, the number of permutations used for hub identification in MEGENA::do.MEGENA
#' @param min.size an integer, the minimum number of nodes allowed in a module
#' @param is.signed a logical, if TRUE signed correlation is used in MEGENA::calculate.correlation
#' @param n.cores an integer, the number of cores to use for parallel processing

#' @return a module object

#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom igraph graph.data.frame
#' @importFrom methods new

#' @export
find_megena_mods <- function(x,
                             method="pearson",
                             fdr.cutoff=0.05,
                             mod.pval=0.05,
                             hub.pval=0.05,
                             cor.perm=10,
                             hub.perm=100,
                             min.size=20,
                             is.signed=FALSE,
                             n.cores = 1
                             ){
  if (!requireNamespace("MEGENA", quietly = TRUE)) {
    stop("Package MEGENA is required. Install with: install.packages('MEGENA')", call. = FALSE)
  }

  # if given access to multiple cores...
  do.par <- n.cores > 1
  if (do.par) {
    cl <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  ijw <- MEGENA::calculate.correlation(x,
                               doPerm=cor.perm,
                               doPar=do.par,
                               num.cores=n.cores,
                               FDR.cutoff=fdr.cutoff,
                               method=method,
                               is.signed=is.signed,
                               output.permFDR=FALSE,
                               output.corTable=FALSE)

  # calculate planar filtered network (PFN)
  el <- MEGENA::calculate.PFN(ijw[,1:3], keep.track=FALSE)
  gd <- igraph::graph.data.frame(el, directed=FALSE)

  # create modules
  megena <- MEGENA::do.MEGENA(gd,
                      mod.pval=mod.pval,
                      hub.pval=hub.pval,
                      remove.unsig=TRUE,
                      min.size=min.size,
                      max.size=igraph::vcount(gd)/2,
                      doPar=do.par,
                      num.cores=n.cores,
                      n.perm=hub.perm,
                      save.output=FALSE)

  # compute modules statistics and summary
  meg.mods <- MEGENA::MEGENA.ModuleSummary(megena,
                                  mod.pvalue=mod.pval,
                                  hub.pvalue=hub.pval,
                                  min.size=min.size,
                                  max.size=igraph::vcount(gd)/2,
                                  output.sig=TRUE)

  # create module object
  meg.mods <- methods::new("module",
                         source = "megena",
                         data.dim = dim(x),
                         overlapping = TRUE,
                         index.list =lapply(meg.mods$modules, function(mod) which(rownames(x) %in% mod)),
                         name.list = meg.mods$modules)

  # return modules
  return(meg.mods)
}

#' Detect modules using ICA then subdivide any modules exceeding a maximum size
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param max.size an integer, the maximum number of nodes allowed in any single module
#' @param n.mods an integer, the number of independent components to extract; if NULL it is selected automatically

#' @return a module object

#' @importFrom methods new

#' @export
pragmatic_modules <- function(x, max.size, n.mods = NULL){
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package fastICA is required. Install with: install.packages('fastICA')", call. = FALSE)
  }
  # perform ICA
  # n initial modules of max size will cover the entire network
  n.comp <- ifelse(is.null(n.mods),
                   ceiling(ncol(t(x)) / max.size), n.mods)

  ICA.results <- fastICA::fastICA(X= x,
                                  n.comp = n.comp,
  )
  comp.score <- abs(ICA.results$S)

  # assign each module the best components
  modules <- apply(comp.score, 1, which.max)
  modules <- split(1:ncol(t(x)), modules)

  # trim to max size
  modules <- lapply(seq_along(modules), function(i){
    mod <- modules[[i]] # get module values
    mod <- mod[order(comp.score[mod,i], decreasing = T)] # order by decreasing score
    mod <- mod[1:max.size] # take only the first max.size values
    mod[!is.na(mod)] # return those that arent na (in case it is shorter than max.size)
  })

  # find unassigned nodes
  unassigned <- which(!(1:ncol(t(x)) %in% unlist(modules)))

  if(length(unassigned) > 0){
    # order based on strength of componet scores
    unassigned <- unassigned[order(rowSums(comp.score[unassigned,]), decreasing = T)]

    mod.sizes <- unlist(lapply(modules, length))
    full.mods <- mod.sizes >= max.size
    for(i in unassigned){
      # determine the best, non full mod
      tmp <- comp.score[i,]
      best.mod <- which(tmp == max(comp.score[i, !full.mods]))
      # add to module and increase its size
      modules[[best.mod]][mod.sizes[best.mod] + 1] <- i # i index will be the same as the order of input data
      mod.sizes[best.mod] = mod.sizes[best.mod] + 1

      # determine if module is full
      if(mod.sizes[best.mod] >= max.size){
        full.mods[best.mod] <- TRUE
      }
    }
  }

  # reformat to get index and score vector from modules index list
  index.vector = vector(mode = "numeric", length = nrow(x))
  score.vector = vector(mode = "numeric", length = nrow(x))
  for(i in seq_along(modules)){
    mod <-  modules[[i]]
    index.vector[mod] <- i
    score.vector[mod] <- comp.score[mod,i]
  }

  # convert to module object and return
  prag.mods <- methods::new("module",
                     source = "pragmatic_modules",
                     data.dim = dim(x),
                     overlapping = FALSE,
                     index.vector = index.vector,
                     score.vector = score.vector,
                     index.list = modules,
                     name.list = lapply(modules, function(mod){rownames(x)[mod]})
  )
  return(prag.mods)
}

###############################
### Expand to Fuzzy Modules ###
###############################

#' Expand non-overlapping modules to fuzzy (overlapping) modules by recruiting nodes correlated with each module's eigengene
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param max.size an integer, the maximum number of nodes allowed in any fuzzy module
#' @param n.pc an integer, the number of principal components used to represent each module's eigengene
#' @param ratio a numeric, the maximum ratio of fuzzy module size to original module size

#' @return a module object

#' @importFrom stats prcomp cor
#' @importFrom methods new

#' @export
eigen_fuzzy_modules <- function(x, input.modules, max.size, n.pc = 2, ratio = 1.5){
  # check if any modules are to large
  if (any(lapply(input.modules@index.list, length) > max.size)) {
    stop("Some modules are to large, increase maxsize.")
  }

  # create new index list (list of nodes in each fuzzy module by index)
  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod.nodes <- input.modules@index.list[[m]]
    # get number of required fuzzy nodes
    f.size <- length(mod.nodes)*ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- f.size - length(mod.nodes)

    # get the modules PC
    mod.PC <- stats::prcomp(t(x[mod.nodes,]), scale. = TRUE)
    # var.exp <- lapply(1:n.pc, function(i){
    #   (mod.PC$sdev[i] / sum(mod.PC$sdev))*100
    # })
    # cat(
    #   "First", n.pc, "principle componets used to define fuzzy modules explain",
    #   round(do.call(sum, var.exp)), "% variance in module", m,
    #   "\n"
    # )

    # define eigen gene based on first n.pc principle components
    mod.eigen <- lapply(1:n.pc, function(i){mod.PC$x[,i]})

    # get genes outside the module the correlated with the eigen gene
    eigen.cor <- lapply(mod.eigen, function(pc){
      abs(apply(x[-mod.nodes,], 1, function(x) stats::cor(x, pc)))
    })
    eigen.cor <- Reduce('+', eigen.cor)

    corRank <- sort(eigen.cor, decreasing = TRUE) # ranked absolute correlation
    fuzzy.nodes <- names(corRank[1:n.fuzzy.nodes]) # n.fuzzy.nodes nodes with the highest ranks

    # convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    # return fuzzy module combining original and fuzzy nodes
    sort(c(mod.nodes, fuzzy.nodes))
  })

  # convert to module object and return
  fuzzy.mods <- methods::new("module",
                             source = paste("eigen_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(m){rownames(x)[m]})
  )
  return(fuzzy.mods)

}

#' Expand non-overlapping modules to fuzzy (overlapping) modules using a thresholded WGCNA adjacency matrix
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param adj a p x p numeric adjacency matrix (e.g. from WGCNA::adjacency) defining pairwise feature similarity
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param max.size an integer, the maximum number of nodes allowed in any fuzzy module
#' @param ratio a numeric, the maximum ratio of fuzzy module size to original module size

#' @return a module object

#' @importFrom methods new

#' @export
adj_fuzzy_modules <- function(x, adj, input.modules, max.size, ratio){
  # check if any modules are to large
  if (any(lapply(input.modules@index.list, length) > max.size)) {
    stop("Some modules are to large, increase maxsize.")
  }

  # create new index list (list of nodes in each fuzzy module by index)
  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod.nodes <- input.modules@index.list[[m]]
    # get number of required fuzzy nodes
    f.size <- length(mod.nodes)*ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- f.size - length(mod.nodes)

    # get matrix of nodes in module adj with nodes outside module
    in.out.adj <- adj[mod.nodes, -mod.nodes]

    # get max adj of node outside module with node inside module
    out.max <- apply(in.out.adj,2 ,max)

    # rank and select fuzzy nodes based on adj
    adj.rank <- order(out.max, decreasing = TRUE)
    fuzzy.nodes <- colnames(in.out.adj)[adj.rank[1:n.fuzzy.nodes]]

    # convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    # return fuzzy module combining original and fuzzy nodes
    sort(c(mod.nodes, fuzzy.nodes))
  })

  # convert to module object and return
  fuzzy.mods <- methods::new("module",
                             source = paste("adj_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(m){rownames(x)[m]})
  )
  return(fuzzy.mods)
}

#' Create overlapping modules by merging pairs of non-overlapping modules based on shared nodes
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param use.eigen a logical, if TRUE module pairs are selected based on eigengene correlation; if FALSE nodewise correlation is used
#' @param best.pairs a logical, if TRUE only the best-correlated module pairs are overlapped; if FALSE all pairs are overlapped

#' @return a module object

#' @importFrom stats cor prcomp
#' @importFrom utils combn
#' @importFrom methods new

#' @export
create_overlap_modules <- function(x, input.modules, use.eigen = TRUE, best.pairs = TRUE){
  # extract input modules index list
  input.indexes <- input.modules@index.list

  # get absolute correlation matrix of all node (absolute only in the case of undirected graphs)
  cor.matrix <- abs(stats::cor(t(x)))

  # create potential pairs list
  p.pairs <- as.data.frame(t(utils::combn(length(input.indexes), 2)))
  if(use.eigen){
    # get eigen genes (first PC)
    mod.eigens <- lapply(input.indexes, function(mod){
      mod.PC <- stats::prcomp(t(x[mod,]), scale. = TRUE)
      mod.PC$x[,1]
    })

    # get eigen covar
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      # get overall mean correlation between each pair of modules
      stats::cor(mod.eigens[[r[1]]], mod.eigens[[r[2]]])

    })
  }else{
    # get score representing the abs correlation between nodes in the modules
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      # get overall mean correlation between each pair of modules
      # mean(cor.matrix[modules[[r[1]]], modules[[r[2]]]])
      max(
        mean(apply(cor.matrix[input.indexes[[r[1]]], input.indexes[[r[2]]]]), 1, max), # row/ col matrix are non symmetric
        mean(apply(cor.matrix[input.indexes[[r[1]]], input.indexes[[r[2]]]]), 2, max) # gives rough measure of overall cor between mods
      )

    })
  }
  p.pairs <- p.pairs[order(p.pairs$cor, decreasing = T),]

  # if only using best pairs
  if(best.pairs){
    # track which items are aleardy in pairs, which rows to keep
    in.pairs = vector(mode = "numeric", length = length(input.indexes))
    keep.rows <- vector(mode = "logical", length = nrow(p.pairs))

    # for each row in the the pair table
    in.pairs[p.pairs[1,1]] <- 1 # set first item to one, so it will remain open to close loop
    for(i in 1:nrow(p.pairs)){
      # check if either item is already in two pairs
      used.twice <- in.pairs[c(p.pairs[i,1],p.pairs[i,2])] >= 2

      # if they are both open
      if(all(!used.twice)){
        # keep the row
        keep.rows[i] <- TRUE

        # iterate the pair
        in.pairs[c(p.pairs[i,1],p.pairs[i,2])] <- in.pairs[c(p.pairs[i,1],p.pairs[i,2])] + 1
      }
    }
    # add a final row: the pair between the two modules that only occure once, the first module and the last (by ordered eigen gene cor)
    p.pairs <- rbind(p.pairs[keep.rows,], c(p.pairs[1,1], which(in.pairs == 1), NA))
  }

  # get overlap mods
  overlap.modules <- lapply(1:nrow(p.pairs), function(i){
    left <- input.indexes[[as.numeric(p.pairs[i,1])]] # get two modules in pair
    right <- input.indexes[[as.numeric(p.pairs[i,2])]]

    # get the correlation of nodes between the two modules
    cor.sub.matrix <- cor.matrix[left, right]

    # get the best nodes from each module: the upper half of row
    left.nodes <- left[order(apply(cor.sub.matrix, 1, max))[1:floor(length(left)/2)]]
    right.nodes <- right[order(apply(cor.sub.matrix, 2, max))[1:floor(length(right)/2)]]

    # return overlap mods
    list(
      index = sort(c(left.nodes,right.nodes)),
      overlap = c(p.pairs[i,1],p.pairs[i,2])
    )
  })
  index.list <- lapply(overlap.modules, function(mod){
    mod$index
  })
  names(index.list) <- unlist(
    lapply(overlap.modules, function(mod){
      as.numeric(paste(mod$overlap[1], mod$overlap[2], sep = "."))
    })
  )

  # combine with initial mods
  index.list <- append(input.modules@index.list, index.list)

  # convert to module object and return
  overlap.modules <- methods::new("module",
                       source = paste("create_overlap_modules", "generated from", input.modules@source),
                       data.dim = dim(x),
                       overlapping = TRUE,
                       index.list = index.list,
                       name.list = lapply(index.list, function(mod){rownames(x)[mod]})
  )
  return(overlap.modules)
}

############################
### Check Module Quality ###
############################

#' Calculate the percentage of nodes assigned to the same module across two module sets, compared by node index
#' @param m1 a module S4 object containing module membership assignments
#' @param m2 a module S4 object containing module membership assignments to compare against m1

#' @return an integer between 0 and 100

#' @export
percent_module_match <- function(m1, m2){
  # get matched sets
  matches <- .match_modules(m1@index.list, m2@index.list, m1@name.list, m2@name.list)
  # count total matches and total nodes in modules
  matched.nodes <- 0
  n.nodes <- min(c(length(unlist(m1@index.vector)), length(m2@index.vector)))
  for (match in matches) {
    matched.nodes = matched.nodes + as.numeric(match[5])
  }
  return((matched.nodes / n.nodes) * 100)
}

#' Find the best-matching module pairs between two sets of modules based on node membership overlap
#' @param m1.nodes an integer vector, the node index assignments from m1
#' @param m2.nodes an integer vector, the node index assignments from m2
#' @param m1.names a character vector, the node names from m1
#' @param m2.names a character vector, the node names from m2

#' @return a list of all the best modules matched based on node membership

#' @keywords internal
.match_modules <- function(m1.nodes, m2.nodes, m1.names, m2.names){
  matches <- vector(mode = "list",
                    length = min(length(m1.nodes), length(m2.nodes)))
  n.matched <- 0
  matched1 <- rep(10^10, length(m1.nodes)) # wont have effect until replaced with actual matches
  matched2 <- rep(10^10, length(m2.nodes))

  while((length(m1.nodes) > n.matched) && (length(m2.nodes) > n.matched)){
    # find the best match
    max.overlap <- -1
    best <- list(NULL, NULL)
    for (i in (1:length(m1.nodes))[-as.numeric(matched1)]) { # iterate through both sets of nodes
      for (j in (1:length(m2.nodes))[-as.numeric(matched2)]) {
        overlap <- length(
          intersect(
            as.numeric(m1.nodes[[i]]),as.numeric(m2.nodes[[j]])
          ) # recording overlap
        )
        if (overlap > max.overlap) {
          max.overlap <- overlap
          best <- c(i, j, m1.names[i], m2.names[j], overlap) # index from set 1, set2, name set 1, set 2, overlap
        }
      }
    } # after finishing iteration, current best match is kept

    matches[[n.matched+1]] <- best
    matched1[n.matched+1] <- best[1] # removes from future consideration
    matched2[n.matched+1] <- best[2]
    n.matched <- n.matched + 1
  }
  return(matches)
}

#' Calculate the fraction of nodes in a module that have more edges within the module than outside of it
#' @param g an igraph object whose node names match the node names stored in test.module
#' @param test.module a module S4 object containing module membership assignments

#' @return an integer between 0 and 100

#' @importFrom igraph as_adjacency_matrix

#' @export
module_contiguity <- function(g, test.module){
  # get module assignments
  index.vector <- test.module@index.vector
  adj <- as.matrix(igraph::as_adjacency_matrix(g))

  # get number nodes with more neighbors within module than between modules
  more.within <- length(which(unlist(
    lapply(1:nrow(adj), function(row){
      node.mod <- index.vector[row]
      edge.mod <- index.vector[which(adj[row,] == 1)]
      length(which(edge.mod == node.mod)) > length(which(edge.mod != node.mod))
    })
  )))
  return(round((more.within / length(g)) * 100, 4))
}

