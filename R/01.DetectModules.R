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
#' @param x a n x p matrix of features
#' @param min.size an integer, the min size parameter for dynamicTreeCut
#' @param min.sft an integer, the min sft used by the pickSoftThreshold function
#' @param beta an integer, the power value used by WGCNA::adjacency
#' @param cor.FN a character, the cor.options WGCNA::adjacency
#' @param powers an integer vector, the powers vector used by WGCNA::pickSoftThreshold
#' @param hclust.method a character, the method used by flashClust
#' @param cut.height an integer between 0 and 1, the cut height used by cutreeDynamic

#' @return a module object

#' @importFrom WGCNA pickSoftThreshold adjacency TOMdist bicor
#' @importFrom flashClust flashClust
#' @importFrom stats as.dist
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom methods new

#' @export
find_WGCNA_mods <- function(x,
                            min.size=10,
                            min.sft=0.85,
                            beta=NULL,
                            cor.FN=c("bicor", "cor"),
                            powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                            hclust.method="average",
                            cut.height = NULL
) {

  # Handle arguments
  cor.FN <- match.arg(cor.FN)

  # Correlation options
  if (cor.FN == "cor") cor.options = list(use="p")
  if (cor.FN == "bicor") cor.options = list(pearsonFallback="individual")

  # Pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data=x,
                                    corFnc=cor.FN,
                                    RsquaredCut=min.sft,
                                    powerVector=powers)

    # Check selected power
    beta <- .sft_check(sft)
  }

  # Construct co-expression similarity
  adj <- WGCNA::adjacency(datExpr=x,
                          power=beta,
                          corFnc=cor.FN,
                          type="unsigned",
                          corOptions=cor.options)

  # Topological overlap dissimilarity transformation
  dis <- WGCNA::TOMdist(adjMat=adj, TOMType="unsigned")

  # Fast hierarchical clustering of dissimilarity
  dendro <- flashClust::flashClust(d=stats::as.dist(dis), method=hclust.method)

  # Module identification using dynamic tree cut algorithm
  index.vector <- dynamicTreeCut::cutreeDynamic(dendro=dendro,
                                           cutHeight = cut.height,
                                           method="hybrid",
                                           distM=dis,
                                           deepSplit=4,
                                           pamRespectsDendro=FALSE,
                                           minClusterSize=min.size)

  #create module object holding the object
  index.vector = as.numeric(index.vector)
  WGCNA.mods  <- methods::new("module",
                     source = "find_WGCNA_mods",
                     data.dim = dim(t(x)),
                     overlapping = FALSE,
                     index.vector = index.vector,
                     index.list = split(1:ncol(x) , index.vector),
                     name.list = split(colnames(x), index.vector)
  )

  return(WGCNA.mods)
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

#' Uses ICA to detect modules from a data matrix
#' @param x a p x n  matrix of features
#' @param n.comp the number of components to search for in the data
#' @param ... other arguments to fastICA

#' @return a module object

#' @importFrom fastICA fastICA
#' @importFrom methods new

#' @export
find_ICA_mods <- function(x, #exprs(eset)
                          n.comp,
                          ...){
  set.seed(123) #removes randomness of outputs
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

#' @return a module object

#' @importFrom stats prcomp cor
#' @importFrom methods new

#' @export
eigen_fuzzy_modules <- function(x, input.modules, max.size){
  index.list <- lapply(input.modules@index.list, function(mod){
    #get number of required fuzzy nodes
    n.fuzzy.nodes <- max.size - length(mod)
    if(n.fuzzy.nodes < 1){
      stop("Warning, modules exist larger than max size")
    }

    #get the modules eigen gene
    mod.PC <- stats::prcomp(t(x[mod,]), scale. = TRUE)
    mod.eigen <- mod.PC$x[,1]

    #get genes outside the module the covary with the eigen gene
    eigen.cor <- abs(apply(x[-mod,], 1, function(x) stats::cor(x, mod.eigen))) #correlation
    corRank <- sort(eigen.cor, decreasing = TRUE) #ranked absolute covariance
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

  #create poteintial pairs list
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

