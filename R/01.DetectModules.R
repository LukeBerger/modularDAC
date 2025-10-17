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
    mod.PC <- prcomp(t(x[mod,]), scale. = TRUE)
    mod.eigen <- mod.PC$x[,1]

    #get genes outside the module the covary with the eigen gene
    eigen.cor <- abs(apply(x[-mod,], 1, function(x) cor(x, mod.eigen))) #correlation
    corRank <- sort(eigen.cor, decreasing = TRUE) #ranked absolute covariance
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
  cor.matrix <- abs(cor(t(x)))

  #create fuzzy modules
  lapply(modules, function(mod){
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
}

## Get a new set of modules, completely overlapping the old set
## To ensure the learned graph can stitch together vary nicely
create_overlap_modules <- function(x, modules, use.eigen = TRUE, best.pairs = TRUE){
  # get absolute correlation matrix of all node (absolute only in the case of undirected graphs)
  cor.matrix <- abs(cor(t(x)))


  #create poteintial pairs list
  p.pairs <- as.data.frame(t(combn(length(modules), 2)))
  if(use.eigen){
    # get eigen genes (first PC)
    mod.eigens <- lapply(modules, function(mod){
      mod.PC <- prcomp(t(x[mod,]), scale. = TRUE)
      mod.PC$x[,1]
    })

    # get eigen covar
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      #get overall mean correlation between each pair of modules
      cor(mod.eigens[[r[1]]], mod.eigens[[r[2]]])

    })
  }else{
    # get score representing the abs correlation between nodes in the modules
    p.pairs$cor <- apply(p.pairs, 1, function(r){
      #get overall mean correlation between each pair of modules
      # mean(cor.matrix[modules[[r[1]]], modules[[r[2]]]])
      max(
        mean(apply(cor.matrix[modules[[r[1]]], modules[[r[2]]]]), 1, max), # row/ col matrix are non symmetric
        mean(apply(cor.matrix[modules[[r[1]]], modules[[r[2]]]]), 2, max) # gives rough measure of overall cor between mods
      )

    })
  }
  p.pairs <- p.pairs[order(p.pairs$cor, decreasing = T),]

  # if only using best pairs
  if(best.pairs){
    #track which items are aleardy in pairs, which rows to keep
    in.pairs = vector(mode = "numeric", length = length(modules))
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
    left <- modules[[as.numeric(p.pairs[i,1])]] # get two modules in pair
    right <- modules[[as.numeric(p.pairs[i,2])]]

    #get the correlation of nodes between the two modules
    cor.sub.matrix <- cor.matrix[left, right]

    #get the best nodes from each module: the upper half of row
    left.nodes <- left[order(apply(cor.sub.matrix, 1, max))[1:floor(length(left)/2)]]
    right.nodes <- right[order(apply(cor.sub.matrix, 2, max))[1:floor(length(right)/2)]]

    #return overlap mods
    list(
      index = sort(c(left.nodes,right.nodes)),
      overalp = c(p.pairs[i,1],p.pairs[i,2])
    )

  })

  return(overlap.modules)
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

# overlap mods
o <- create_overlap_modules(x, lapply(p, function(x) x$index))


