# library(WGCNA)
# library(fastICA)

## Takes a n x p  matrix of  data and list of args to cuttreeDynamic
## Returns dendo plot and list of modules
## Based on mods.detect from shin
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
  args <- as.list(environment())
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
  dendro <- flashClust::flashClust(d=as.dist(dis), method=hclust.method)

  # Module identification using dynamic tree cut algorithm
  modules <- dynamicTreeCut::cutreeDynamic(dendro=dendro,
                                           cutHeight = cut.height,
                                           method="hybrid",
                                           distM=dis,
                                           deepSplit=4,
                                           pamRespectsDendro=FALSE,
                                           minClusterSize=min.size)
  return(as.numeric(modules))
}

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


## Simple module detection based on ICA components
## takes a n x p matrix of data, a number of components to search for, and args to fastICA
## returns a list of modules assigned based on maximum abs competent score
find_ICA_mods <- function(x, #exprs(eset)
                          n.comp,
                          ...){
  set.seed(123) #removes randomness of outputs
  ICA.results <- fastICA::fastICA(X= as.matrix(x),
                                  n.comp = n.comp,
                                  ...
  )
  modules <- apply(abs(ICA.results$S), 1, which.max)
  return(modules)
}

## Detect initial modules based on ICA components
## Convert to n.mods modules of max.size nodes
## returns a list contian module indexes, node names, and ica scores
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

  # add additional information
  modules <- lapply(seq_along(modules), function(i){
    mod <-  modules[[i]]
    names(mod) <- colnames(t(x))[mod]
    list(
      index = mod,
      name = colnames(t(x))[mod],
      ica.score = comp.score[mod,i]
    )
  })

  return(modules)
}

## Create overlapping fuzzy modules by adding nodes to prexisiting modules up to a max size,
## new nodes are selected based on correlation with the modules eigen gene (first pc)
eigen_fuzzy_modules <- function(x, modules, max.size){
  lapply(modules, function(mod){
    #get number of required fuzzy nodes
    n.fuzzy.nodes <- max.size - length(mod)

    #get the modules eigen gene
    modPC <- prcomp(t(x[mod,]), scale. = TRUE)
    modEigen <- modPC$x[,1]

    #get genes outside the module the covary with the eigen gene
    eigenCor <- abs(apply(x[-mod,], 1, function(x) cor(x, modEigen))) #correlation
    corRank <- sort(eigenCor, decreasing = TRUE) #ranked absolute covariance
    fuzzy.nodes <- names(corRank[1:n.fuzzy.nodes]) #n.fuzzy.nodes nodes with the highest ranks

    #convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    #return fuzzy module combining original and fuzzy nodes
    sort(c(mod, fuzzy.nodes))
  })
}

## Create overlapping fuzzy modules by adding nodes to prexisiting modules up to a max size,
## new nodes are selected based on correlation with individual nodes in the module
nodewise_fuzzy_modules <- function(x, modules, max.size){
  #get cor matrix
  corMat <- abs(cor(t(x)))

  #create fuzzy modules
  lapply(modules, function(mod){
    #get number of required fuzzy nodes
    n.fuzzy.nodes <- max.size - length(mod)

    #get only covariance of genes in the module with genes not in the module
    corSub <- corMat[mod,-mod] #mod-genes by not-mod-genes matrix

    #get column maxes for absolute covariance and use them to select fuzzy nodes
    cMx <- apply(corSub,2 ,max)
    corRank <- order(cMx, decreasing = TRUE)
    fuzzy.nodes <- colnames(corSub)[corRank[1:n.fuzzy.nodes]]

    #convert fuzzy nodes to numerics (stored naturally as names)
    fuzzy.nodes <- which(rownames(x) %in% fuzzy.nodes)
    names(fuzzy.nodes) <- rownames(x[fuzzy.nodes,])

    #return fuzzy module combining original and fuzzy nodes
    sort(c(mod, fuzzy.nodes))

  })
}



###########
## TEST ###
###########

source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/00.SimulateGraphs.R")

# make data
er <- make_modular_graph()
x <- sim_graph_data(er, n.samples =100)

# learn mods
w <- find_WGCNA_mods(t(x), cor.FN = "bicor")
i <- find_ICA_mods(x, 3)
p <- pragmatic_modules(x,n.mods = 3, max.size = 60)

# fuzzy mods
f <- eigen_fuzzy_modules(x, lapply(p, function(x) x$index), 80)
f <- nodewise_fuzzy_modules(x, lapply(p, function(x) x$index), 80)
