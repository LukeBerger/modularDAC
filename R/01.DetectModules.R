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
#' @param assign.by a character, how unassigned nodes are recruited to modules: 'adjacency' uses WGCNA adjacency, 'eigengene' uses correlation with module eigengenes

#' @return a list containing: the WGCNA thresholded adj matrix, the initial modules after merge/iteration, the final modules after max.size fitting

#' @importFrom WGCNA pickSoftThreshold adjacency TOMdist bicor mergeCloseModules
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
                            merge = FALSE,
                            merging.cut = 0.2,
                            iterate = TRUE,
                            assign.by = c("adjacency", "eigengene")
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
  assign.by <- match.arg(assign.by)

  # treat an unconstrained max.size as Inf so the size checks below are explicit
  if (is.null(max.size)) max.size <- Inf

  # merge / iterate relationship
  if (merge && iterate) {
    message("Note: using both 'merge' and 'iterate' may merge then immediately split some modules; use with caution.")
  }

  # correlation options
  cor.options <- if (cor.FN == "cor") list(use = "p") else list(pearsonFallback = "individual")

  # drop features with no variance (they cannot be placed in any module)
  rv <- MatrixGenerics::rowVars(x)
  if (any(rv == 0)) {
    rm.r <- which(rv == 0)
    warning("Some features have zero variance and were removed; they could not be placed in any module.")
    message("Removed rows: ", paste(rm.r, collapse = ", "))
    x <- x[-rm.r, , drop = FALSE]
  }

  # transpose once (WGCNA expects samples x features) and reuse
  datExpr <- t(x)

  # pick soft threshold via scale-free fit
  if (is.null(beta)) {
    sft <- WGCNA::pickSoftThreshold(data = datExpr,
                                    corFnc = cor.FN,
                                    RsquaredCut = min.sft,
                                    powerVector = powers)
    beta <- .sft_check(sft)
  }

  # minimum number of modules needed to fit all nodes within max.size
  min.mods <- max(1, ceiling(nrow(x) / max.size))

  # construct co-expression similarity, TOM distance, and hierarchical clustering
  adj <- WGCNA::adjacency(datExpr = datExpr,
                          power = beta,
                          corFnc = cor.FN,
                          type = "unsigned",
                          corOptions = cor.options)
  dis <- WGCNA::TOMdist(adjMat = adj, TOMType = "unsigned")
  dendro <- flashClust::flashClust(d = stats::as.dist(dis), method = hclust.method)

  # initial module identification via dynamic tree cut
  message("Generating initial modules...")
  initial.index.vector <- dynamicTreeCut::cutreeDynamic(dendro = dendro,
                                                        cutHeight = cut.height,
                                                        method = "hybrid",
                                                        distM = dis,
                                                        deepSplit = 4,
                                                        pamRespectsDendro = FALSE,
                                                        minClusterSize = min.size)

  # lower the cut height until enough modules are produced
  lowered <- .lower_cut_height(index.vector = initial.index.vector,
                               dendro = dendro,
                               dis = dis,
                               min.mods = min.mods,
                               min.size = min.size,
                               cut.height = cut.height,
                               max.size = max.size)
  initial.index.vector <- lowered$index.vector
  cut.height <- lowered$cut.height

  # merge similar modules based on eigengene similarity
  if (merge) {
    initial.index.vector <- .merge_modules(index.vector = initial.index.vector,
                                           datExpr = datExpr,
                                           cor.FN = cor.FN,
                                           cor.options = cor.options,
                                           merging.cut = merging.cut,
                                           min.mods = min.mods,
                                           max.size = max.size)
  }

  # split any modules that exceed max.size
  if (iterate) {
    initial.index.vector <- .split_large_modules(index.vector = initial.index.vector,
                                                 dis = dis,
                                                 max.size = max.size,
                                                 min.size = min.size,
                                                 hclust.method = hclust.method,
                                                 cut.height = cut.height)
  }

  # assign remaining unassigned nodes to a module
  message("Assigning unassigned nodes by ", assign.by, " method...")
  diag(adj) <- 0
  modified.index.vector <- .assign_unassigned(index.vector = initial.index.vector,
                                              adj = adj,
                                              x = x,
                                              method = assign.by)

  # fit modules to max size via node trading
  n.too.big <- sum(table(modified.index.vector) > max.size)
  if (n.too.big > 0) {
    message(n.too.big, " module(s) larger than ", max.size, " nodes; performing node trading based on adjacency...")
    modified.index.vector <- .fit_to_max_size(index.vector = modified.index.vector,
                                              adj = adj,
                                              max.size = max.size)
  }

  # build module objects
  initial.index.vector <- as.numeric(initial.index.vector)
  modified.index.vector <- as.numeric(modified.index.vector)
  initial.mods <- methods::new("module",
                               source = "find_WGCNA_mods initial mods",
                               data.dim = dim(x),
                               overlapping = FALSE,
                               index.vector = initial.index.vector,
                               index.list = split(1:nrow(x), initial.index.vector),
                               name.list = split(rownames(x), initial.index.vector)
  )
  final.mods <- methods::new("module",
                             source = "find_WGCNA_mods traded mods",
                             data.dim = dim(x),
                             overlapping = FALSE,
                             index.vector = modified.index.vector,
                             index.list = split(1:nrow(x), modified.index.vector),
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
    message("Using the following power: ", beta)
  } else {
    message("Optimal power selected: ", beta)
  }
  return(beta)
}

#' Helper that trims modules exceeding a maximum size by trading nodes to smaller modules
#' @param index.vector an integer vector of length p, assigning each node to a module (0 = unassigned)
#' @param max.size an integer, the maximum number of nodes allowed in a module
#' @param adj a p x p numeric similarity matrix (required when method = 'adjacency')
#' @param x a numeric matrix with p features (rows) and n samples (columns) (required when method = 'eigengene')
#' @param method a character, how a giving node is scored against receiving modules: 'adjacency' uses the similarity matrix, 'eigengene' uses correlation with module eigengenes
#' @param verbose a logical, if TRUE progress messages are printed for each trade (useful for debugging)

#' @return an integer vector of length p with the updated module assignments

#' @importFrom stats prcomp cor
#' @keywords internal
.fit_to_max_size <- function(index.vector,
                             max.size,
                             adj = NULL,
                             x = NULL,
                             method = c("adjacency", "eigengene"),
                             verbose = FALSE){
  method <- match.arg(method)
  if (method == "adjacency" && is.null(adj)) {
    stop("method = 'adjacency' requires a similarity matrix 'adj'.")
  }
  if (method == "eigengene" && is.null(x)) {
    stop("method = 'eigengene' requires the data matrix 'x'.")
  }

  mod.size <- table(index.vector)
  giving <- as.numeric(names(which(mod.size > max.size)))
  n.trades <- 0

  while (length(giving) > 0) {
    giving.nodes <- which(index.vector %in% giving)
    receiving <- as.numeric(names(which(mod.size < max.size)))
    if (length(receiving) == 0) {
      warning("No module has room to receive nodes; some modules remain above max.size.")
      break
    }

    if (method == "adjacency") {
      # move the giving node most similar to any receiving module
      receiving.nodes <- which(index.vector %in% receiving)
      trade.scores <- adj[giving.nodes, receiving.nodes, drop = FALSE]
      gi <- which.max(matrixStats::rowMaxs(trade.scores))
      trade.giving <- giving.nodes[gi]
      trade.to.module <- index.vector[ receiving.nodes[ which.max(adj[trade.giving, receiving.nodes]) ] ]
    } else {
      # score each giving node by its absolute correlation with each receiving module eigengene
      eig <- vapply(receiving, function(m) {
        stats::prcomp(t(x[which(index.vector == m), , drop = FALSE]), scale. = TRUE)$x[, 1]
      }, numeric(ncol(x)))
      cors <- abs(stats::cor(t(x[giving.nodes, , drop = FALSE]), eig))
      cors[is.na(cors)] <- 0
      best.col <- max.col(cors, ties.method = "first")
      gi <- which.max(cors[cbind(seq_along(giving.nodes), best.col)])
      trade.giving <- giving.nodes[gi]
      trade.to.module <- receiving[best.col[gi]]
    }

    if (verbose) { message("Trading ", trade.giving, " to module ", trade.to.module) }
    index.vector[trade.giving] <- trade.to.module

    mod.size <- table(index.vector)
    giving <- as.numeric(names(which(mod.size > max.size)))
    n.trades <- n.trades + 1
  }
  message("All modules fit within max size after ", n.trades, " trades.")

  return(index.vector)
}

#' Helper to find_WGCNA_mods that lowers the dynamic-tree-cut height until a
#' minimum number of modules is produced
#' @param index.vector an integer vector of length p, the modules from the initial cut (0 = unassigned)
#' @param dendro a dendrogram (flashClust/hclust object) used for the cut
#' @param dis a p x p TOM distance matrix passed to cutreeDynamic
#' @param min.mods an integer, the minimum number of modules required
#' @param min.size an integer, the minimum number of nodes allowed in a module
#' @param cut.height a numeric or NULL, the starting cut height; recomputed from the dendrogram when NULL
#' @param max.size a numeric, the maximum module size (used only for messaging)
#' @param min.cut.height a numeric, the floor below which the cut height will not be lowered

#' @return a list with the updated index.vector and the cut.height that produced it

#' @importFrom stats quantile
#' @keywords internal
.lower_cut_height <- function(index.vector,
                              dendro,
                              dis,
                              min.mods,
                              min.size,
                              cut.height = NULL,
                              max.size = Inf,
                              min.cut.height = 0) {
  n.mods <- length(unique(index.vector[index.vector != 0]))
  if (n.mods >= min.mods) {
    message("Initial WGCNA call produced ", n.mods, " modules.")
    return(list(index.vector = index.vector, cut.height = cut.height))
  }

  message("Only ", n.mods, " modules produced in first cut; ", min.mods,
          " required with max.size == ", max.size, ". Lowering cut height...")

  # recompute the default height used by cutreeDynamic when none was supplied
  if (is.null(cut.height)) {
    qnt5 <- stats::quantile(dendro$height, probs = 0.05)
    cut.height <- qnt5 + 0.99 * (max(dendro$height) - qnt5)
  }

  while (n.mods < min.mods) {
    cut.height <- cut.height - 0.01
    # stop once the cut height hits its floor, even if min.mods was not reached
    if (cut.height <= min.cut.height) {
      message("Cut height reached its floor (", min.cut.height,
              ") without producing ", min.mods, " modules; proceeding with the current cut.")
      break
    }
    message("Lowering cut height to ", round(cut.height, 4), " and generating new modules...")
    index.vector <- dynamicTreeCut::cutreeDynamic(dendro = dendro,
                                                  cutHeight = cut.height,
                                                  method = "hybrid",
                                                  distM = dis,
                                                  deepSplit = 4,
                                                  pamRespectsDendro = FALSE,
                                                  minClusterSize = min.size)
    # bail out if cutting is no longer producing meaningful clusters
    if (all(index.vector == 0) && cut.height < 0.8) {
      message("Cut height reached ", round(cut.height, 4),
              " and is no longer producing meaningful clusters. ",
              "Consider increasing max.size or using another method to generate modules.")
      stop("No modules found after lowering cut height")
    }
    n.mods <- length(unique(index.vector[index.vector != 0]))
  }
  message("Final cut height ", round(cut.height, 4), " produced ", n.mods, " modules.")
  list(index.vector = index.vector, cut.height = cut.height)
}

#' Helper to find_WGCNA_mods that splits modules exceeding max.size by
#' re-clustering each oversized module
#' @param index.vector an integer vector of length p assigning each node to a module (0 = unassigned)
#' @param dis a p x p TOM distance matrix
#' @param max.size a numeric, the maximum number of nodes allowed in a module
#' @param min.size an integer, the minimum number of nodes allowed in a module
#' @param hclust.method a character, the agglomeration method passed to flashClust
#' @param cut.height a numeric or NULL, the cut height passed to cutreeDynamic

#' @return an integer vector of length p with oversized modules split where possible

#' @keywords internal
.split_large_modules <- function(index.vector,
                                 dis,
                                 max.size,
                                 min.size,
                                 hclust.method = "average",
                                 cut.height = NULL) {
  # find modules that are too large (ignoring the unassigned group, 0)
  to.big <- as.numeric(names(which(table(index.vector) > max.size)))
  to.big <- to.big[to.big != 0]
  if (length(to.big) == 0) return(index.vector)

  message("Modules ", paste(to.big, collapse = ", "),
          " are too large; attempting to split...")

  # re-cluster each oversized module on its own distance submatrix
  for (m in to.big) {
    m.nodes <- which(index.vector == m)
    m.dis <- dis[m.nodes, m.nodes, drop = FALSE]

    dendro <- flashClust::flashClust(d = stats::as.dist(m.dis), method = hclust.method)
    m.index.vector <- dynamicTreeCut::cutreeDynamic(dendro = dendro,
                                                    cutHeight = cut.height,
                                                    method = "hybrid",
                                                    distM = m.dis,
                                                    deepSplit = 4,
                                                    pamRespectsDendro = FALSE,
                                                    minClusterSize = min.size,
                                                    verbose = FALSE)

    # only accept the split if it is mostly assigned and yields at least two modules
    pass <- TRUE
    if (sum(m.index.vector == 0) > length(m.index.vector) / 2) {
      pass <- FALSE
      message("Module ", m, " failed to split due to a high number of unassigned nodes.")
    }
    if (length(unique(m.index.vector[m.index.vector != 0])) < 2) {
      pass <- FALSE
      message("Module ", m, " failed to split into at least two new modules.")
    }

    if (pass) {
      message("Splitting module ", m, " into ",
              length(unique(m.index.vector[m.index.vector != 0])), " new modules.")
      # offset new labels by the current max so they do not collide with existing modules
      index.vector[m.nodes] <- ifelse(m.index.vector == 0,
                                      0, m.index.vector + max(index.vector))
    }
  }
  index.vector
}

#' Helper to find_WGCNA_mods that assigns unassigned nodes (module 0) to an
#' existing module, either by WGCNA adjacency or by correlation with module eigengenes
#' @param index.vector an integer vector of length p, assigning each node to a module (0 = unassigned)
#' @param adj a p x p numeric adjacency matrix produced by WGCNA::adjacency (diagonal expected to be 0); required when method = 'adjacency'
#' @param x a numeric matrix with p features (rows) and n samples (columns); required when method = 'eigengene'
#' @param method a character, how to score candidate modules for an unassigned node: 'adjacency' uses the maximum adjacency to an assigned node, 'eigengene' uses the maximum absolute correlation with a module eigengene
#' @param start.threshold a numeric, the initial score a node must exceed to be recruited to a module
#' @param step a numeric, the amount the threshold is lowered after a pass that recruits no nodes
#' @param verbose a logical, if TRUE the number of nodes recruited at each threshold is printed

#' @return an integer vector of length p with unassigned nodes assigned where possible

#' @importFrom stats prcomp cor
#' @keywords internal
.assign_unassigned <- function(index.vector,
                               adj = NULL,
                               x = NULL,
                               method = c("adjacency", "eigengene"),
                               start.threshold = 0.95,
                               step = 0.05,
                               verbose = FALSE){
  method <- match.arg(method)
  if(method == "adjacency" && is.null(adj)){
    stop("method = 'adjacency' requires an adjacency matrix 'adj'.")
  }
  if(method == "eigengene" && is.null(x)){
    stop("method = 'eigengene' requires the data matrix 'x'.")
  }

  # iteratively recruit each unassigned node to its best-scoring module,
  # lowering the threshold whenever a pass recruits nothing
  threshold <- start.threshold
  repeat{
    unassigned <- which(index.vector == 0)
    if(length(unassigned) == 0) break

    assigned <- which(index.vector != 0)
    if(length(assigned) == 0) break # no module to grow from

    if(method == "adjacency"){
      # for each unassigned node, find its single best-connected assigned node
      score <- adj[unassigned, assigned, drop = FALSE]
      best.col <- max.col(score, ties.method = "first")
      best.score <- score[cbind(seq_along(unassigned), best.col)]
      best.module <- index.vector[assigned[best.col]]
    }else{
      # represent each module by its eigengene (first PC of member expression),
      # then score each unassigned node by its absolute correlation to each eigengene
      mods <- sort(unique(index.vector[assigned]))
      eigengenes <- vapply(mods, function(m){
        stats::prcomp(t(x[which(index.vector == m), , drop = FALSE]), scale. = TRUE)$x[, 1]
      }, numeric(ncol(x)))
      score <- abs(stats::cor(t(x[unassigned, , drop = FALSE]), eigengenes))
      score[is.na(score)] <- 0 # genes with no signal cannot be scored
      best.col <- max.col(score, ties.method = "first")
      best.score <- score[cbind(seq_along(unassigned), best.col)]
      best.module <- mods[best.col]
    }

    # recruit only the nodes whose best score clears the current threshold
    recruit <- best.score > threshold
    if(any(recruit)){
      index.vector[unassigned[recruit]] <- best.module[recruit]
      if(verbose){
        message("Assigned ", sum(recruit), " node(s) at ", method, " threshold ", threshold)
      }
    }else{
      # nothing clears the threshold; lower it, stopping once it reaches 0
      if(threshold <= 0) break
      threshold <- round(threshold - step, 2)
    }
  }

  # any node still in module 0 could not be scored against any module
  leftover <- sum(index.vector == 0)
  if(leftover > 0){
    warning(leftover, " node(s) could not be connected to any module and remain unassigned (module 0).")
  }

  return(index.vector)
}

#' Helper that merges similar modules by eigengene similarity (WGCNA::mergeCloseModules)
#' @param index.vector an integer vector of length p assigning each node to a module (0 = unassigned)
#' @param datExpr a samples x features numeric matrix (i.e. t(x))
#' @param cor.FN a character, the correlation function ('bicor' or 'cor')
#' @param cor.options a list of correlation options passed to mergeCloseModules
#' @param merging.cut a numeric, the eigengene dissimilarity threshold for merging
#' @param min.mods an integer, the minimum number of modules to retain
#' @param max.size a numeric, modules are only merged while all are below this size

#' @return an integer vector of length p with merged module assignments

#' @importFrom WGCNA mergeCloseModules
#' @keywords internal
.merge_modules <- function(index.vector, datExpr, cor.FN, cor.options, merging.cut, min.mods, max.size) {
  message("Merging modules based on eigengene similarity...")
  n.mods <- length(unique(index.vector[index.vector != 0]))
  n.merges <- 0
  module.eigengenes <- NULL # reused across rounds to skip recomputing module eigengenes
  while (n.mods > min.mods && all(table(index.vector) < max.size)) {
    merged <- WGCNA::mergeCloseModules(
      exprData = datExpr,
      colors = index.vector,
      MEs = module.eigengenes,
      unassdColor = 0,
      corFnc = cor.FN,
      corOptions = cor.options,
      cutHeight = merging.cut
    )
    if (all(index.vector == merged$colors)) {
      break
    }
    index.vector <- merged$colors
    module.eigengenes <- merged$newMEs
    n.mods <- length(unique(index.vector[index.vector != 0]))
    n.merges <- n.merges + 1
  }
  if (n.merges == 0) {
    message("No modules were small/similar enough to merge; using the initial modules.")
  } else {
    message("Merged modules in ", n.merges, " round(s), resulting in ", n.mods, " modules.")
  }
  index.vector
}

#' Simple helper that extracts a p x p feature-similarity matrix for node trading
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param method a character, 'WGCNA' (soft-thresholded co-expression adjacency) or 'ARACNE' (mutual information matrix from minet, the basis of ARACNE)
#' @param cor.FN a character, the correlation function for WGCNA ('bicor' or 'cor')
#' @param min.sft a numeric, the minimum R-squared for WGCNA soft-threshold selection
#' @param beta an integer or NULL, the WGCNA soft-thresholding power (auto-selected when NULL)
#' @param powers an integer vector, candidate powers for WGCNA::pickSoftThreshold

#' @return a p x p numeric similarity matrix, rows/cols ordered as the rows of x

#' @importFrom WGCNA pickSoftThreshold adjacency
#' @keywords internal
.feature_similarity <- function(x,
                                method = c("WGCNA", "ARACNE"),
                                cor.FN = c("bicor", "cor"),
                                min.sft = 0.85,
                                beta = NULL,
                                powers = c(seq(1, 10, by = 1), seq(12, 20, by = 2))) {
  method <- match.arg(method)
  cor.FN <- match.arg(cor.FN)
  if (method == "WGCNA") {
    if (!requireNamespace("WGCNA", quietly = TRUE)) {
      stop("Package WGCNA is required. Install with: install.packages('WGCNA')", call. = FALSE)
    }
    cor.options <- if (cor.FN == "cor") list(use = "p") else list(pearsonFallback = "individual")
    datExpr <- t(x)
    if (is.null(beta)) {
      sft <- WGCNA::pickSoftThreshold(data = datExpr, corFnc = cor.FN,
                                      RsquaredCut = min.sft, powerVector = powers)
      beta <- .sft_check(sft)
    }
    sim <- WGCNA::adjacency(datExpr = datExpr, power = beta, corFnc = cor.FN,
                            type = "unsigned", corOptions = cor.options)
  } else {
    if (!requireNamespace("minet", quietly = TRUE)) {
      stop("Package minet is required for method = 'ARACNE'. Install with: install.packages('minet')", call. = FALSE)
    }
    sim <- minet::build.mim(dataset = t(x))
  }
  rownames(sim) <- colnames(sim) <- rownames(x)
  sim
}

#' Helper to find_ICA_mods that splits modules exceeding max.size by re-running ICA
#' @param index.vector an integer vector of length p assigning each node to a module
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param max.size a numeric, the maximum number of nodes allowed in a module
#' @param ... additional arguments passed to fastICA::fastICA

#' @return an integer vector of length p with oversized modules split where possible

#' @keywords internal
.split_large_modules_ICA <- function(index.vector, x, max.size, ...) {
  to.big <- as.numeric(names(which(table(index.vector) > max.size)))
  to.big <- to.big[to.big != 0]
  if (length(to.big) == 0) return(index.vector)

  message("Modules ", paste(to.big, collapse = ", "), " are too large; attempting to split with ICA...")
  for (m in to.big) {
    m.nodes <- which(index.vector == m)
    # extract just enough components to bring sub-modules within max.size
    n.split <- min(ceiling(length(m.nodes) / max.size), length(m.nodes), ncol(x))
    if (n.split < 2) next

    ica <- fastICA::fastICA(X = as.matrix(x[m.nodes, , drop = FALSE]), n.comp = n.split, ...)
    sub.index <- apply(abs(ica$S), 1, which.max)
    if (length(unique(sub.index)) < 2) {
      message("Module ", m, " failed to split into at least two new modules.")
      next
    }
    message("Splitting module ", m, " into ", length(unique(sub.index)), " new modules.")
    # offset new labels by the current max so they do not collide with existing modules
    index.vector[m.nodes] <- sub.index + max(index.vector)
  }
  index.vector
}

#' Detect co-expression modules from a data matrix using Independent Component Analysis (ICA)
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param n.comp an integer, the number of independent components (modules) to extract
#' @param max.size an integer or NULL, the maximum number of nodes allowed in a module (NULL = unconstrained)
#' @param merge a logical, if TRUE similar modules are merged by eigengene similarity (WGCNA::mergeCloseModules)
#' @param merging.cut a numeric, the eigengene dissimilarity threshold for merging
#' @param iterate a logical, if TRUE modules exceeding max.size are split by re-running ICA on their features
#' @param cor.FN a character, the correlation function used for merging and WGCNA-based trading ('bicor' or 'cor')
#' @param trade.by a character, how oversized modules are trimmed to max.size: 'adjacency' (a learned similarity matrix) or 'eigengene' (correlation with module eigengenes)
#' @param network a character, the method used to learn the similarity matrix when trade.by = 'adjacency': 'WGCNA' or 'ARACNE'
#' @param ... additional arguments passed to fastICA::fastICA

#' @return a list with: the learned similarity matrix (or NULL), the initial modules, and the final modules after splitting/trading

#' @importFrom methods new
#' @importFrom stats prcomp cor

#' @export
find_ICA_mods <- function(x,
                          n.comp,
                          max.size = NULL,
                          merge = FALSE,
                          merging.cut = 0.2,
                          iterate = FALSE,
                          cor.FN = c("bicor", "cor"),
                          trade.by = c("adjacency", "eigengene"),
                          network = c("WGCNA", "ARACNE"),
                          ...) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package fastICA is required. Install with: install.packages('fastICA')", call. = FALSE)
  }
  cor.FN <- match.arg(cor.FN)
  trade.by <- match.arg(trade.by)
  network <- match.arg(network)
  if (is.null(max.size)) max.size <- Inf
  cor.options <- if (cor.FN == "cor") list(use = "p") else list(pearsonFallback = "individual")

  # decompose features into independent components and assign each to its top component
  ICA.results <- fastICA::fastICA(X = as.matrix(x), n.comp = n.comp, ...)
  initial.index.vector <- apply(abs(ICA.results$S), 1, which.max)
  score.vector <- apply(abs(ICA.results$S), 1, max)
  message("ICA produced ", length(unique(initial.index.vector)), " modules.")

  # minimum number of modules needed to satisfy max.size
  min.mods <- max(1, ceiling(nrow(x) / max.size))

  # merge similar modules based on eigengene similarity
  if (merge) {
    initial.index.vector <- .merge_modules(index.vector = initial.index.vector,
                                           datExpr = t(x),
                                           cor.FN = cor.FN,
                                           cor.options = cor.options,
                                           merging.cut = merging.cut,
                                           min.mods = min.mods,
                                           max.size = max.size)
  }

  # split modules that exceed max.size by re-running ICA on their features
  if (iterate) {
    initial.index.vector <- .split_large_modules_ICA(index.vector = initial.index.vector,
                                                     x = x, max.size = max.size, ...)
  }

  # trade nodes to bring every module within max.size
  modified.index.vector <- initial.index.vector
  similarity <- NULL
  if (sum(table(modified.index.vector) > max.size) > 0) {
    if (trade.by == "adjacency") {
      message("Learning a ", network, " similarity matrix for adjacency-based trading...")
      similarity <- .feature_similarity(x, method = network, cor.FN = cor.FN)
    }
    message("Trimming oversized modules to max.size by ", trade.by, " trading...")
    modified.index.vector <- .fit_to_max_size(index.vector = modified.index.vector,
                                              max.size = max.size,
                                              adj = similarity,
                                              x = x,
                                              method = trade.by)
  }

  # build module objects
  initial.index.vector <- as.numeric(initial.index.vector)
  modified.index.vector <- as.numeric(modified.index.vector)
  initial.mods <- methods::new("module",
                               source = "find_ICA_mods initial mods",
                               data.dim = dim(x),
                               overlapping = FALSE,
                               index.vector = initial.index.vector,
                               score.vector = score.vector,
                               index.list = split(1:nrow(x), initial.index.vector),
                               name.list = split(rownames(x), initial.index.vector))
  final.mods <- methods::new("module",
                             source = "find_ICA_mods traded mods",
                             data.dim = dim(x),
                             overlapping = FALSE,
                             index.vector = modified.index.vector,
                             score.vector = score.vector,
                             index.list = split(1:nrow(x), modified.index.vector),
                             name.list = split(rownames(x), modified.index.vector))

  return(list(
    similarity = similarity,
    initial.mods = initial.mods,
    final.mods = final.mods
  ))
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



############################
### Check Module Quality ###
############################

#' Compare two module sets from the same dataset by counting shared nodes
#' @param m1 a module S4 object containing module membership assignments
#' @param m2 a module S4 object containing module membership assignments to compare against m1, generated from the same dataset

#' @return a list with three elements: 'overlap.matrix', an n1 x n2 matrix whose (i, j) cell counts the nodes shared between module i of m1 (row) and module j of m2 (column); 'best.matches', a three-column data frame giving each greedily matched module pair ('set1', 'set2') and their percent (Jaccard) overlap ('percent'); and 'overall.overlap', a single score from 0 (no agreement) to 100 (identical sets) summarising agreement between the two sets (a size-weighted mean of the best-match Jaccards, with unmatched modules counted as zero)

#' @importFrom methods is validObject

#' @export
module_match <- function(m1, m2){
  # --- validate inputs ---
  # both inputs must be module objects
  if(!methods::is(m1, "module") || !methods::is(m2, "module")){
    stop("Both 'm1' and 'm2' must be 'module' objects.")
  }
  # enforce the module class invariants (matching index / name list lengths, etc.)
  methods::validObject(m1)
  methods::validObject(m2)

  # both sets must describe the same dataset: identical node-index universe ...
  nodes1 <- sort(unique(unlist(m1@index.list)))
  nodes2 <- sort(unique(unlist(m2@index.list)))
  if(!identical(nodes1, nodes2)){
    stop("'m1' and 'm2' cover different node indices; they must be generated from the same dataset.")
  }
  # ... identical node names ...
  names1 <- sort(unique(unlist(m1@name.list)))
  names2 <- sort(unique(unlist(m2@name.list)))
  if(!identical(names1, names2)){
    stop("'m1' and 'm2' cover different node names; they must be generated from the same dataset.")
  }
  # ... and a consistent index-to-name mapping across the two sets
  map1 <- unique(paste(unlist(m1@index.list), unlist(m1@name.list), sep = "\t"))
  map2 <- unique(paste(unlist(m2@index.list), unlist(m2@name.list), sep = "\t"))
  if(!setequal(map1, map2)){
    stop("'m1' and 'm2' disagree on the index-to-name mapping; they must be generated from the same dataset.")
  }

  # --- set up ---
  idx1 <- m1@index.list
  idx2 <- m2@index.list
  n1 <- length(idx1)
  n2 <- length(idx2)
  # module labels come from the index list names, falling back to position
  lab1 <- names(idx1); if(is.null(lab1)) lab1 <- as.character(seq_len(n1))
  lab2 <- names(idx2); if(is.null(lab2)) lab2 <- as.character(seq_len(n2))
  size1 <- unname(lengths(idx1))
  size2 <- unname(lengths(idx2))

  # --- OUTPUT 1: overlap matrix of shared-node counts ---
  overlap.matrix <- matrix(0L, nrow = n1, ncol = n2,
                           dimnames = list(lab1, lab2))
  for(i in seq_len(n1)){
    for(j in seq_len(n2)){
      overlap.matrix[i, j] <- length(intersect(idx1[[i]], idx2[[j]]))
    }
  }

  # --- OUTPUT 2: best-match table ---
  # greedily pair modules using the precomputed overlap matrix
  matched <- .match_modules(overlap.matrix)
  inter <- matched[, "overlap"]
  union.size <- size1[matched[, "row"]] + size2[matched[, "col"]] - inter
  pair.jaccard <- ifelse(union.size > 0, inter / union.size, 0)
  best.matches <- data.frame(
    set1    = lab1[matched[, "row"]],
    set2    = lab2[matched[, "col"]],
    percent = 100 * pair.jaccard,
    stringsAsFactors = FALSE
  )

  # --- OUTPUT 3: overall overlap metric ---
  # size-weighted mean Jaccard over the best matches; modules in the larger set
  # left without a partner contribute a Jaccard of 0, weighted by their size, so
  # that differing module counts are penalised
  unmatched.rows <- setdiff(seq_len(n1), matched[, "row"])
  unmatched.cols <- setdiff(seq_len(n2), matched[, "col"])
  weight.total <- sum(union.size) + sum(size1[unmatched.rows]) + sum(size2[unmatched.cols])
  overall.overlap <- if(weight.total > 0){
    100 * sum(pair.jaccard * union.size) / weight.total
  }else{
    0
  }

  # return the overlap matrix, best-match table and overall metric
  return(
    list(
      overlap.matrix  = overlap.matrix,
      best.matches    = best.matches,
      overall.overlap = overall.overlap
    )
  )
}

#' Greedily pair modules between two sets using a precomputed overlap matrix
#' @param overlap.matrix an n1 x n2 integer matrix of shared-node counts; rows are set-1 modules and columns are set-2 modules

#' @return an integer matrix with one row per matched pair and columns 'row' (set-1 module index), 'col' (set-2 module index) and 'overlap' (shared-node count)

#' @keywords internal
.match_modules <- function(overlap.matrix){
  n1 <- nrow(overlap.matrix)
  n2 <- ncol(overlap.matrix)
  k <- min(n1, n2)

  matched <- matrix(NA_integer_, nrow = k, ncol = 3,
                    dimnames = list(NULL, c("row", "col", "overlap")))
  avail.row <- rep(TRUE, n1)
  avail.col <- rep(TRUE, n2)

  # repeatedly take the largest remaining overlap, fixing one module from each
  # set per step, until the smaller set is exhausted
  for(m in seq_len(k)){
    masked <- overlap.matrix
    masked[!avail.row, ] <- -1L
    masked[, !avail.col] <- -1L

    best <- which.max(masked)            # column-major linear index of the max cell
    i <- ((best - 1L) %% n1) + 1L
    j <- ((best - 1L) %/% n1) + 1L

    matched[m, ] <- c(i, j, overlap.matrix[i, j])
    avail.row[i] <- FALSE
    avail.col[j] <- FALSE
  }

  return(matched)
}

#' Measure how self-contained each module is within a graph
#' @param g an igraph object whose node names match the node names stored in test.module
#' @param test.module a module S4 object containing module membership assignments

#' @return a list with two elements: 'overall.contiguity', a single score from 0 to 100 giving the edge-weighted percentage of module-incident edges that stay within their module; and 'module.edges', a data frame with one row per module (row names are the module labels) holding 'within.edges' (edges with both endpoints in the module), 'between.edges' (edges with exactly one endpoint in the module) and 'percent.within' (100 * within / (within + between))

#' @importFrom igraph is_igraph as_edgelist V
#' @importFrom methods is

#' @export
module_contiguity <- function(g, test.module){
  # --- validate inputs ---
  if(!igraph::is_igraph(g)){
    stop("'g' must be an igraph object.")
  }
  if(!methods::is(test.module, "module")){
    stop("'test.module' must be a 'module' object.")
  }
  # module membership is matched to the graph by node name
  v.names <- igraph::V(g)$name
  module.names <- test.module@name.list
  if(!all(unlist(module.names) %in% v.names)){
    stop("'test.module' contains node names that are not present in 'g'.")
  }

  # module labels used for the output rows
  labels <- names(module.names)
  if(is.null(labels)) labels <- as.character(seq_along(module.names))

  # graph edges as integer endpoint pairs (each undirected edge appears once)
  edges <- igraph::as_edgelist(g, names = FALSE)

  # --- per-module within / between edge counts ---
  counts <- lapply(module.names, function(nm){
    # membership indicator over the graph's vertices
    member <- v.names %in% nm
    if(nrow(edges) == 0){
      return(c(within = 0, between = 0))
    }
    a.in <- member[edges[, 1]]
    b.in <- member[edges[, 2]]
    c(within  = sum(a.in & b.in),       # both endpoints inside the module
      between = sum(xor(a.in, b.in)))   # exactly one endpoint inside the module
  })

  within.edges   <- vapply(counts, function(z) as.numeric(z["within"]),  numeric(1))
  between.edges  <- vapply(counts, function(z) as.numeric(z["between"]), numeric(1))
  total.edges    <- within.edges + between.edges
  percent.within <- ifelse(total.edges > 0, 100 * within.edges / total.edges, NA_real_)

  module.edges <- data.frame(
    within.edges   = within.edges,
    between.edges  = between.edges,
    percent.within = percent.within,
    row.names      = labels,
    stringsAsFactors = FALSE
  )

  # --- overall contiguity: edge-weighted fraction of module edges kept internal ---
  overall.contiguity <- if(sum(total.edges) > 0){
    100 * sum(within.edges) / sum(total.edges)
  }else{
    NA_real_
  }

  return(
    list(
      overall.contiguity = overall.contiguity,
      module.edges       = module.edges
    )
  )
}


#' Calculate within-module feature correlations and plot their distribution
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param input.modules a module S4 object containing the module assignments to evaluate
#' @param cor.method a character, the correlation method passed to stats::cor; one of 'pearson', 'spearman' or 'kendall'

#' @return a list containing: 'correlation', a named list holding the upper triangular feature correlation matrix for each module (diagonal and lower triangle set to NA); and 'plot', a ggplot boxplot of all correlation values grouped by module

#' @importFrom stats cor

#' @export
module_correlation <- function(x, input.modules, cor.method = c("pearson", "spearman", "kendall")){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package ggplot2 is required. Install with: install.packages('ggplot2')", call. = FALSE)
  }

  # handle arguments
  cor.method <- match.arg(cor.method)

  # validate the module object against the input data
  .module_check(x, input.modules)

  # extract the module index list and assign each module a label
  index.list <- input.modules@index.list
  mod.names <- names(index.list)
  if(is.null(mod.names)){
    mod.names <- as.character(seq_along(index.list))
  }

  # for each module get the upper triangular feature correlation matrix
  cor.list <- lapply(index.list, function(mod){
    # correlate features (rows) across samples (columns)
    cm <- stats::cor(t(x[mod, , drop = FALSE]), method = cor.method)
    # keep only the upper triangle, blanking the diagonal and lower triangle
    cm[lower.tri(cm, diag = TRUE)] <- NA
    cm
  })
  names(cor.list) <- mod.names

  # gather all upper triangular correlations into one long data frame for plotting
  plot.df <- do.call(rbind, lapply(seq_along(cor.list), function(i){
    vals <- cor.list[[i]][upper.tri(cor.list[[i]])]
    # modules with a single feature contribute no pairwise correlations
    if(length(vals) == 0){return(NULL)}
    data.frame(
      module = mod.names[i],
      correlation = vals,
      stringsAsFactors = FALSE
    )
  }))
  # preserve module ordering along the x axis
  plot.df$module <- factor(plot.df$module, levels = mod.names)

  # boxplot of the within-module feature correlations
  cor.boxplot <- ggplot2::ggplot(
    plot.df,
    ggplot2::aes(x = module, y = correlation)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(
      x = "Module",
      y = "Feature correlation",
      title = paste("Within-module feature correlation:", input.modules@source)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

  # return the upper triangular correlation matrices and the boxplot
  return(list(
    correlation = cor.list,
    plot = cor.boxplot
  ))
}
