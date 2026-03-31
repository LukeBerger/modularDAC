# Declare module object
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
           #check that each module has the same number of feature indexes and names
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

#' Returns true modules from a simulated graph with V(g)$modules
#' @param g a igraph object with values for each node in the V(g)$module field

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V

#' @export
true_modules <- function(g){
  methods::new("module",
                source = "True Modules",
                overlapping = FALSE,
                index.vector = igraph::V(g)$module,
                index.list = split(1:length(g) , igraph::V(g)$module),
                name.list = split(igraph::V(g)$name , igraph::V(g)$module)
  )
}

#' Returns fuzzy moduls based on neighborhood in true graph
#' @param m a module object
#' @param g a igraph object with values for each node in the V(g)$module field

#' @return a module object

#' @importFrom methods new
#' @importFrom igraph V neighborhood

#' @export
true_fuzzy <- function(m, g){
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
  methods::new("module",
               source = "True Module Fuzzy",
               overlapping = TRUE,
               index.list = f.index.list,
               name.list = lapply(f.index.list, function(m) igraph::V(g)$name[m])
  )
}


#'Checks if a module object is properly formated and throws errors if not
#' @param x data matrix used to learn the modules
#' @param m a module object

#' @return TRUE, if all checks passed

#' @keywords internal
.module_check <- function(x, m){
  if(m@overlapping){
    #check that overlaps exist between all modules
    if(!all(
      unlist(
        lapply(seq_along(m@index.list), function(i){
          any(unlist(m@index.list[i]) %in% unlist(m@index.list[-i]))
        })
      )
    )){stop(paste(m@source, "has mules with incomplete overlaps"))}
  }else{
    #check feature number matches input data
    if(length(m@index.vector) != nrow(x)){
      stop(paste(m@source, "produced a index vector with the incorrect number of features"))
    }
    #check that there are no overlaps
    if(any(
      unlist(
        lapply(seq_along(m@index.list), function(i){
          any(unlist(m@index.list[i]) %in% unlist(m@index.list[-i]))
        })
      )
    )){stop(paste(m@source, "has mules with overlaps"))}
  }

  #check that all feature names in modules come from the data
  if(!all(unlist(m@name.list) %in% rownames(x))){
    stop(paste(m@source, "feature names do not match input data"))
  }
  #check that each module has the same number of feature indexes and names
  if(
    !all(
      unlist(
        lapply(seq_along(m@index.list), function(i) length(m@index.list[[i]]) == length(m@name.list[[i]]))
      )
    )
  ){stop(paste(m@source, "produced differnet length index and name lists"))}
  return(TRUE)
}

#' Uses WGCNA to detect modules from a data matrix
#' @param x a p x n matrix of features
#' @param min.size an integer, the min size of modules produced by dynamicTreeCut
#' @param max.size an integer, the maximum size of modules
#' @param min.sft an integer, the min sft used by the pickSoftThreshold function
#' @param beta an integer, the power value used by WGCNA::adjacency
#' @param cor.FN a character, the cor.options WGCNA::adjacency
#' @param powers an integer vector, the powers vector used by WGCNA::pickSoftThreshold
#' @param hclust.method a character, the method used by flashClust
#' @param cut.height an integer between 0 and 1, the cut height used by cutreeDynamic
#' @param merge a boolean, whether or not to attempt to merge modules using WGCNA::mergeCloseModules
#' @param merging.cut a numeric between 0 and 1, the cutheight for WGCNA::mergeCloseModules

#' @return a module object

#' @importFrom WGCNA pickSoftThreshold adjacency TOMdist bicor mergeCloseModules
#' @importFrom flashClust flashClust
#' @importFrom stats as.dist
#' @importFrom dynamicTreeCut cutreeDynamic
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
                            merging.cut = 0.2
) {

  # Handle arguments
  cor.FN <- match.arg(cor.FN)

  # Correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=t(x),
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # Check selected power
    beta <- .sft_check(sft)
  }

  # If max size is constrained
  if(!is.null(max.size)){
    # Determine the minimum number of modules needed to fit all nodes
    min.mods <- ceiling(nrow(x) / max.size)
  }else{ min.mods <- 1}

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=t(x),
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # Topological overlap dissimilarity transformation
  dis <- WGCNA::TOMdist(adjMat=adj, TOMType="unsigned")

  # Fast hierarchical clustering of dissimilarity
  dendro <- flashClust::flashClust(d=stats::as.dist(dis), method=hclust.method)

  # Module identification using dynamic tree cut algorithm
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
      cat("Initial WGCNA call produced", n.mods, "\n")
  }


  # If merging
  if(merge){
    cat("Merging modules based on eigen gene similarity... \n")
    n.merges <- 0
    # While there are at least min.mods modules and they are all bellow maxsize
    while(n.mods > min.mods && all(table(initial.index.vector < max.size )) ){
      # Merge similar modules (single iteration)
      merged <- WGCNA::mergeCloseModules(
        exprData=t(x),
        colors=initial.index.vector,
        # MEs= ifelse(exists("merged"), merged$newMEs, NULL),
        unassdColor=0,
        corFnc=cor.FN,
        corOptions=cor.options,
        cutHeight=merging.cut
      )

      # if nothing changes, break loop
      if(all(initial.index.vector == merged$colors)){
        break
      }
      # else, use merged modules and restart
      initial.index.vector <- merged$colors
      n.mods <- length(unique(initial.index.vector[initial.index.vector != 0]))
      n.merges <- n.merges + 1
    }
    if(n.merges == 0){
      cat("No modules were similiar enough to merge, using initial modules... \n")
    }else{
      cat("Merged" , n.merges, "modules, resulting in", n.mods, "new modules. \n")
    }

  }

  # Assign remaining unassigned nodes based on adj...
  cat("Assinging unassinged nodes based on WGNCA adjacency... \n")
  diag(adj) <- 0

  # iterative expansion by best neighbor with threshold
  modified.index.vector <- initial.index.vector
  last.mods <- modified.index.vector
  unassigned <- which(modified.index.vector == 0)
  threshold = 0.95
  while(length(unassigned) > 0){ #there might be a more efficient way to check this
    #for each assigned node
    assigned <- which(modified.index.vector != 0)
    unassigned <- which(modified.index.vector == 0)
    for(node in assigned){
      #find its best neighbor based on wgcna adj
      best <- unassigned[which.max(adj[node, unassigned])]
      #if that adj is above the current threshold
      if(length(best) !=0){
        if(adj[node,best] > threshold){
          #assign node
          modified.index.vector[best] <- modified.index.vector[node]
        }
      }
    }
    #if no changes
    if(all(last.mods == modified.index.vector)){
      #and not at min threshold
      if(threshold >= 0){
        #lower threshold
        threshold <- round(threshold - 0.05, 2)
      }else{
        #break if at min threshold and no longer changing
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

  # Create module object
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
  traded.mods  <- methods::new("module",
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
    traded.mods = traded.mods
  ))
}

#' Helper to WGCNA calling functions to select soft thresholding power in WGCNA::adjacency
#' @param sft output of WGCNA::pickSoftThreshold

#' @return an int power level

#' @keywords internal
.sft_check <- function(sft) {
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    beta <- 6 # Default
    cat("Using the following power:", beta, "\n")
  } else {
    cat("Optimal power selected:", beta, "\n")
  }
  return(beta)
}

#' Helper to find_WGCNA_mods that forces all modules to be within some max size
#' @param index.vector index vector defining a module for each node or 0 for unassigned
#' @param adj adjacency matrix produced by WGCNA::adjacency
#' @param max.size integer for max number of nodes in each module

#' @return index vector with each module fitted to max size

#' @keywords internal

.fit_to_max_size <- function(
    index.vector,
    adj,
    max.size

){

  # find modules which are two large
  mod.size <- table(index.vector)
  if(any(mod.size > max.size)){
    # determine which will give and receive nodes
    giving = as.numeric(names(which(mod.size > max.size)))
    receiving = as.numeric(names(which(mod.size < max.size)))

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

      # set giving nodes module to receiving nods
      index.vector[ trade.giving ] <- index.vector[ trade.receiving ]

      # update module status
      mod.size = table(index.vector)
      giving = which(mod.size > max.size)
      receiving = which(mod.size < max.size)
      giving.nodes <- which(index.vector %in% giving)
      receiving.nodes <- which(index.vector %in% receiving)
      n.trades <- n.trades + 1
    }
    cat("All modules fit within max size after",  n.trades, "trades. \n")
  }

  return(index.vector)
}


#' Uses ICA to detect modules from a data matrix
#' @param x a p x n  matrix of features
#' @param n.comp the number of components to search for in the data
#' @param ... other arguments to fastICA

#' @return a module object

#' @importFrom fastICA fastICA
#' @importFrom methods new

#' @export
find_ICA_mods <- function(x,
                          n.comp,
                          ...){
  ICA.results <- fastICA::fastICA(X= as.matrix(x),
                                  n.comp = n.comp,
                                  ...
  )
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

#' Uses MCL to detect modules from a data matrix
#' @param x a n x p matrix of features
#' @param beta an integer, the power value for WGCNA::adjacency
#' @param cut a numeric < 1, the cut threshold for converting from correlation to binary adjacency
#' @param expansion an integer, the expansion value of the mcl function
#' @param inflation an integer, the inflation value of the mcl function
#' @param cor.fn the correlation function to use in WGCNA::adjacency
#' @param iter an integer, the maximum number of iterations in the mcl function

#' @return a module object

#' @importFrom WGCNA adjacency
#' @importFrom MCL mcl
#' @importFrom methods new

#' @export
find_mcl_mods <- function(x,
                          beta=6,
                          cut=0.05,
                          expansion=2,
                          inflation=1,
                          cor.fn="cor",
                          iter=1000){
  # Compute unsigned and scaled correlation matrix
  mat <- WGCNA::adjacency(datExpr=x,
                          power=beta,
                          corFnc=cor.fn,
                          type="unsigned")

  # Convert to adjacency matrix with cut point
  adj <- mat
  adj[mat >= cut] <- 1
  adj[mat < cut] <- 0
  diag(adj) <- 0

  # Markov clustering
  clusters <- MCL::mcl(x=as.matrix(adj),
                  addLoops=TRUE,
                  allow1=FALSE,
                  ESM=FALSE,
                  max.iter=iter,
                  expansion=expansion,
                  inflation=inflation)

  # Create module object
  mcl.mods <-methods::new("module",
                         source = "mcl",
                         data.dim = dim(t(x)),
                         overlapping = FALSE,
                         index.vector = clusters$Cluster,
                         index.list =split(1:ncol(x) , clusters$Cluster),
                         name.list = split(colnames(x), clusters$Cluster))

  # Return Modules
  return(mcl.mods)
}


#' Uses megena to detect modules from a data matrix
#' @param x a n x p matrix of features
#' @param method method for the calculate.correlation function
#' @param fdr.cutoff fdr cut off value in calculate.correlation
#' @param mod.pval module p value for the do.MEGENA function
#' @param hub.pval hub p value for the do.MEGENA function
#' @param cor.perm number of  permutations in calculate.correlation
#' @param hub.perm the number of hub permutations in do.MEGENA
#' @param min.size the minimum number of nodes in each module
#' @param is.signed indicate using signed/unsigned correlation in calculate.correlation
#' @param n.cores number of cores to use for multithreading

#' @return a module object

#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom MEGENA calculate.correlation calculate.PFN do.MEGENA
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

  # If given access to multiple cores...
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

  # Calculate planar filtered network (PFN)
  el <- MEGENA::calculate.PFN(ijw[,1:3], keep.track=FALSE)
  gd <- igraph::graph.data.frame(el, directed=FALSE)

  # Create modules
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

  # Compute modules statistics and summary
  meg.mods <- MEGENA::MEGENA.ModuleSummary(megena,
                                  mod.pvalue=mod.pval,
                                  hub.pvalue=hub.pval,
                                  min.size=min.size,
                                  max.size=igraph::vcount(gd)/2,
                                  output.sig=TRUE)

  # Create module object
  meg.mods <- methods::new("module",
                         source = "megena",
                         data.dim = dim(x),
                         overlapping = TRUE,
                         index.list =lapply(meg.mods$modules, function(mod) which(rownames(x) %in% mod)),
                         name.list = meg.mods$modules)

  # Return modules
  return(meg.mods)
}

#' Uses ICA to detect modules from a data matrix then reassigned nodes to fit all modules to a max size
#' @param x  a p x n matrix of features
#' @param max.size the maximum number of members in a data matrix
#' @param n.mods the number of modules to search for in the data

#' @return a module object

#' @importFrom fastICA fastICA
#' @importFrom methods new

#' @export
pragmatic_modules <- function(x, max.size, n.mods = NULL){
  # Perform ICA
  # n initial modules of max size will cover the entire network
  n.comp <- ifelse(is.null(n.mods),
                   ceiling(ncol(t(x)) / max.size), n.mods)


  ICA.results <- fastICA::fastICA(X= x,
                                  n.comp = n.comp,
  )
  comp.score <- abs(ICA.results$S)

  #assign each module the best components
  modules <- apply(comp.score, 1, which.max)
  modules <- split(1:ncol(t(x)), modules)

  #trim to max size
  modules <- lapply(seq_along(modules), function(i){
    mod <- modules[[i]] #get module values
    mod <- mod[order(comp.score[mod,i], decreasing = T)] # order by decreasing score
    mod <- mod[1:max.size] # take only the first max.size values
    mod[!is.na(mod)] # return those that arent na (in case it is shorter than max.size)
  })

  #find unassigned nodes
  unassigned <- which(!(1:ncol(t(x)) %in% unlist(modules)))

  if(length(unassigned) > 0){
    #order based on strength of componet scores
    unassigned <- unassigned[order(rowSums(comp.score[unassigned,]), decreasing = T)]


    mod.sizes <- unlist(lapply(modules, length))
    full.mods <- mod.sizes >= max.size
    for(i in unassigned){
      #determine the best, non full mod
      tmp <- comp.score[i,]
      best.mod <- which(tmp == max(comp.score[i, !full.mods]))
      #add to module and increase its size
      modules[[best.mod]][mod.sizes[best.mod] + 1] <- i # i index will be the same as the order of input data
      mod.sizes[best.mod] = mod.sizes[best.mod] + 1

      #determine if module is full
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

#' Creates overlapping sets of fuzzy modules based on correlation of nodes outside the module with its eigen gene
#' @param x a p x n matrix of features
#' @param input.modules a module object
#' @param max.size The maximum number of members in a data matrix
#' @param n.pc number of principle componets use to define correlation of fuzzy nodes to the module
#' @param ratio the ratio of original nodes in the module to nodes in the fuzzy module

#' @return a module object

#' @importFrom stats prcomp cor
#' @importFrom methods new

#' @export
eigen_fuzzy_modules <- function(x, input.modules, max.size, n.pc = 2, ratio = 1.5){
  # check if any modules are to large
  if (any(lapply(input.modules@index.list, length) > max.size)) {
    stop("Some modules are to large, increase maxsize.")
  }


  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod <- input.modules@index.list[[m]]
    #get number of required fuzzy nodes
    f.size <- length(mod)*ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- f.size - length(mod)

    #get the modules PC
    mod.PC <- stats::prcomp(t(x[mod,]), scale. = TRUE)
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

    #get genes outside the module the correlated with the eigen gene
    eigen.cor <- lapply(mod.eigen, function(pc){
      abs(apply(x[-mod,], 1, function(x) stats::cor(x, pc)))
    })
    eigen.cor <- Reduce('+', eigen.cor)


    corRank <- sort(eigen.cor, decreasing = TRUE) #ranked absolute correlation
    fuzzy.nodes <- names(corRank[1:n.fuzzy.nodes]) #n.fuzzy.nodes nodes with the highest ranks

    #convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    #return fuzzy module combining original and fuzzy nodes
    sort(c(mod, fuzzy.nodes))
  })

  # convert to module object and return
  fuzzy.mods <- methods::new("module",
                             source = paste("eigen_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(mod){rownames(x)[mod]})
  )
  return(fuzzy.mods)

}

#' Creates overlapping sets of fuzzy modules based on nodewise correlation of nodes outside the module with any of the nodes within it
#' @param x a p x n matrix of features
#' @param input.modules a module object
#' @param max.size The maximum number of members in a data matrix

#' @return a module object

#' @importFrom stats cor
#' @importFrom methods new

#' @export
nodewise_fuzzy_modules <- function(x, input.modules, max.size){
  #get cor matrix
  cor.matrix <- abs(stats::cor(t(x)))

  #create fuzzy modules
  modules <- lapply(input.modules@index.list, function(mod){
    #get number of required fuzzy nodes
    n.fuzzy.nodes <- max.size - length(mod)

    #get only covariance of genes in the module with genes not in the module
    cor.sub.matrix <- cor.matrix[mod,-mod] #mod-genes by not-mod-genes matrix

    #get column maxes for absolute covariance and use them to select fuzzy nodes
    column.maxes <- apply(cor.sub.matrix,2 ,max)
    corRank <- order(column.maxes, decreasing = TRUE)
    fuzzy.nodes <- colnames(cor.sub.matrix)[corRank[1:n.fuzzy.nodes]]

    #convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    #return fuzzy module combining original and fuzzy nodes
    sort(c(mod, fuzzy.nodes))

  })

  # convert to module object and return
  fuzzy_modules <- methods::new("module",
                       source = paste("nodewise_fuzzy_modules", "generated from", input.modules@source),
                       data.dim = dim(x),
                       overlapping = TRUE,
                       index.list = modules,
                       name.list = lapply(modules, function(mod){rownames(x)[mod]})
  )
  return(fuzzy_modules)
}

#' Create a set of overlapping modules from data matrix and strict, none overlapping modules
#' @param x a p x n matrix of features
#' @param input.modules a module object
#' @param use.eigen whether or not to use the eigen gene to determine overlaps, otherwise nodewise correlation is used
#' @param best.pairs whether or not to overlap the best pairs of modules based on eigen genes only, otherwise all pairs of models are overlapped

#' @return a module object

#' @importFrom stats cor prcomp
#' @importFrom utils combn
#' @importFrom methods new

#' @export
create_overlap_modules <- function(x, input.modules, use.eigen = TRUE, best.pairs = TRUE){
  #extract input modules index list
  input.indexes <- input.modules@index.list

  # get absolute correlation matrix of all node (absolute only in the case of undirected graphs)
  cor.matrix <- abs(stats::cor(t(x)))

  #create potential pairs list
  p.pairs <- as.data.frame(t(utils::combn(length(input.indexes), 2)))
  if(use.eigen){
    # get eigen genes (first PC)
    mod.eigens <- lapply(input.indexes, function(mod){
      mod.PC <- stats::prcomp(t(x[mod,]), scale. = TRUE)
      mod.PC$x[,1]
    })

    # get eigen covar
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      #get overall mean correlation between each pair of modules
      stats::cor(mod.eigens[[r[1]]], mod.eigens[[r[2]]])

    })
  }else{
    # get score representing the abs correlation between nodes in the modules
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      #get overall mean correlation between each pair of modules
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
    #track which items are aleardy in pairs, which rows to keep
    in.pairs = vector(mode = "numeric", length = length(input.indexes))
    keep.rows <- vector(mode = "logical", length = nrow(p.pairs))

    #for each row in the the pair table
    in.pairs[p.pairs[1,1]] <- 1 #set first item to one, so it will remain open to close loop
    for(i in 1:nrow(p.pairs)){
      #check if either item is already in two pairs
      used.twice <- in.pairs[c(p.pairs[i,1],p.pairs[i,2])] >= 2

      #if they are both open
      if(all(!used.twice)){
        #keep the row
        keep.rows[i] <- TRUE

        #iterate the pair
        in.pairs[c(p.pairs[i,1],p.pairs[i,2])] <- in.pairs[c(p.pairs[i,1],p.pairs[i,2])] + 1
      }
    }
    # add a final row: the pair between the two modules that only occure once, the first module and the last (by ordered eigen gene cor)
    p.pairs <- rbind(p.pairs[keep.rows,], c(p.pairs[1,1], which(in.pairs == 1), NA))
  }

  #get overlap mods
  overlap.modules <- lapply(1:nrow(p.pairs), function(i){
    left <- input.indexes[[as.numeric(p.pairs[i,1])]] # get two modules in pair
    right <- input.indexes[[as.numeric(p.pairs[i,2])]]

    #get the correlation of nodes between the two modules
    cor.sub.matrix <- cor.matrix[left, right]

    #get the best nodes from each module: the upper half of row
    left.nodes <- left[order(apply(cor.sub.matrix, 1, max))[1:floor(length(left)/2)]]
    right.nodes <- right[order(apply(cor.sub.matrix, 2, max))[1:floor(length(right)/2)]]

    #return overlap mods
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


