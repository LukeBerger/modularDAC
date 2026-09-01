######################
### Define Modules ###
######################

# silence R CMD check NOTE for ggplot2 aes() non-standard evaluation in module_correlation
utils::globalVariables(c("module", "correlation"))

# declare module object
#'
#' The optional `core.list` slot records, for each module, the subset of its
#' nodes that the module *owns* ("core" nodes, whose Markov blanket is expected
#' to lie inside the module so their edges are trustworthy) as opposed to
#' auxiliary nodes that are only present to condition on for someone else. When
#' populated, the core sets must partition the node universe: every node is a
#' core (owned) node of exactly one module. This ownership is what makes the
#' divide-and-conquer stitching well defined (see \code{.connect_subgraphs}).
#' @importFrom methods setClass
setClass("module",
         slots = list(
           source = "character",
           data.dim = "numeric",
           overlapping = "logical",
           index.vector = "numeric",
           score.vector = "numeric",
           index.list = "list",
           name.list = "list",
           core.list = "list"
         ),
         prototype = list(core.list = list()),
         validity = function(object){
           # check that each module has the same number of feature indexes and names
           if(!all(lengths(object@index.list) == lengths(object@name.list))){
             return("index list and name list must be the same length")
           }
           # if not overlapping, no index may appear in more than one module
           if(!object@overlapping){
             if(any(duplicated(unlist(object@index.list)))){
               return("Non-overlapping modules must not share indices")
             }
           }
           # core.list, when supplied, records node ownership for stitching
           if(length(object@core.list) > 0){
             if(length(object@core.list) != length(object@index.list)){
               return("core.list must have one entry per module")
             }
             # each module's core must be a subset of its own nodes
             in.mod <- mapply(function(cr, idx) all(cr %in% idx),
                              object@core.list, object@index.list)
             if(!all(in.mod)){
               return("each module's core.list must be a subset of its index.list")
             }
             # cores must partition the node universe: owned by exactly one module
             all.core <- unlist(object@core.list)
             if(any(duplicated(all.core))){
               return("core nodes must be owned by exactly one module (no duplicate cores)")
             }
             if(!setequal(all.core, unique(unlist(object@index.list)))){
               return("every node must be a core (owned) node of exactly one module")
             }
           }

           TRUE
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
  module <- igraph::V(g)$module
  index.list <- split(seq_len(igraph::vcount(g)), module)
  methods::new("module",
                source = "True Modules",
                overlapping = FALSE,
                index.vector = module,
                index.list = index.list,
                name.list = split(igraph::V(g)$name, module),
                # a non-overlapping partition owns all of its own nodes
                core.list = index.list
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
  node.names <- igraph::V(g)$name
  # collect all nodes within 2 graph hops of each module
  f.index.list <- lapply(m@index.list, function(idx){
    sort(unique(unlist(
      igraph::neighborhood(g, order = 2, nodes = idx, mode = "all", mindist = 0)
    )))
  })
  # build fuzzy module object with expanded index list; the original (pre-growth)
  # nodes are the core nodes this module owns, the recruited neighbours are auxiliary
  methods::new("module",
               source = "True Module Fuzzy",
               overlapping = TRUE,
               index.list = f.index.list,
               name.list = lapply(f.index.list, function(idx) node.names[idx]),
               core.list = m@index.list
  )
}

#' Validate a module object against the data matrix and throw errors if checks fail
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param m a module S4 object containing module membership assignments

#' @return TRUE, if all checks passed

#' @keywords internal
.module_check <- function(x, m){
  if(!m@overlapping){
    # check feature number matches input data
    if(length(m@index.vector) != nrow(x)){
      stop(paste(m@source, "produced an index vector with the incorrect number of features"))
    }
    # no index may appear in more than one module
    if(any(duplicated(unlist(m@index.list)))){
      stop(paste(m@source, "has modules with overlaps"))
    }
  }

  # check that all feature names in modules come from the data
  if(!all(unlist(m@name.list) %in% rownames(x))){
    stop(paste(m@source, "feature names do not match input data"))
  }
  # check that each module has the same number of feature indexes and names
  if(!all(lengths(m@index.list) == lengths(m@name.list))){
    stop(paste(m@source, "produced different length index and name lists"))
  }
  # when node ownership is recorded, cores must partition the node universe:
  # every node is a core (owned) node of exactly one module
  if(length(m@core.list) > 0){
    all.core <- unlist(m@core.list)
    if(any(duplicated(all.core))){
      stop(paste(m@source, "has core nodes owned by more than one module"))
    }
    if(!setequal(all.core, unique(unlist(m@index.list)))){
      stop(paste(m@source, "cores do not cover every node exactly once"))
    }
  }
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
#' @param iterate a logical, if TRUE oversized modules are recursively split so every module fits within max.size
#' @param assign.by a character, how unassigned nodes are recruited to modules: 'adjacency' uses WGCNA adjacency, 'eigengene' uses correlation with module eigengenes
#' @param deep.split an integer 0-4 passed to \code{dynamicTreeCut::cutreeDynamic}, controlling how aggressively the dendrogram is split; 4 (the default, and the previously hard-coded value) splits most aggressively, 0 least
#' @param min.natural.mods an integer, the minimum number of modules the dynamic tree cut must produce before the partition is accepted; if fewer are found and \code{retry.betas} is TRUE the soft-threshold power is re-tried (see \code{retry.betas}). Fewer than 1 is always an error, since a cut that assigns every feature to the unassigned group cannot be repaired downstream
#' @param retry.betas a logical; if TRUE (default) and the soft-threshold power was chosen automatically, a natural cut yielding fewer than \code{min.natural.mods} modules is re-tried at the remaining \code{powers}, best scale-free fit first, keeping the first power that resolves the network into modules. Ignored when \code{beta} is supplied explicitly

#' @return a list containing: 'wgcna.adj', the WGCNA co-expression adjacency matrix; 'initial.mods', the natural modules after merging and unassigned-node recruitment; and 'final.mods', those modules after recursive splitting to satisfy max.size

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
                            assign.by = c("adjacency", "eigengene"),
                            deep.split = 4,
                            min.natural.mods = 2L,
                            retry.betas = TRUE
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

  # treat an unconstrained max.size as Inf (so the size checks below are explicit)
  if (is.null(max.size)) max.size <- Inf

  # merge / iterate relationship
  if (merge && iterate) {
    message("Note: using both 'merge' and 'iterate' may merge then immediately split some modules; use with caution.")
  }

  # correlation options
  cor.options <- if (cor.FN == "cor") list(use = "p") else list(pearsonFallback = "individual")

  # drop features with no variance (they cannot be placed in any module)
  rm.r <- which(MatrixGenerics::rowVars(x) == 0)
  if (length(rm.r) > 0) {
    warning("Some features have zero variance and were removed; they could not be placed in any module.")
    message("Removed rows: ", paste(rm.r, collapse = ", "))
    x <- x[-rm.r, , drop = FALSE]
  }

  # transpose once (WGCNA expects samples x features) and reuse
  t.x <- t(x)

  # pick soft threshold via scale-free fit
  sft <- NULL
  auto.beta <- is.null(beta)
  if (auto.beta) {
    # corOptions and networkType are passed explicitly so the power is chosen
    # under the SAME correlation and network settings the adjacency below is
    # built with. They previously diverged: pickSoftThreshold fell back to its
    # own default corOptions = list(use = "p") while adjacency used
    # list(pearsonFallback = "individual") for bicor, so the selected power was
    # fitted to a slightly different matrix than the one it was then applied to.
    sft <- WGCNA::pickSoftThreshold(data = t.x,
                                    corFnc = cor.FN,
                                    corOptions = cor.options,
                                    networkType = "unsigned",
                                    RsquaredCut = min.sft,
                                    powerVector = powers)
    beta <- .sft_check(sft)
  }

  # minimum number of modules to retain when merging (a floor for .merge_modules)
  min.mods <- max(1, ceiling(nrow(x) / max.size))

  # build the co-expression adjacency, TOM distance, dendrogram and natural
  # (dynamic tree cut) partition at one soft-threshold power
  cut.at.power <- function(b) {
    adj <- WGCNA::adjacency(datExpr = t.x,
                            power = b,
                            corFnc = cor.FN,
                            type = "unsigned",
                            corOptions = cor.options)
    dis <- WGCNA::TOMdist(adjMat = adj, TOMType = "unsigned")
    dendro <- flashClust::flashClust(d = stats::as.dist(dis), method = hclust.method)
    iv <- dynamicTreeCut::cutreeDynamic(dendro = dendro,
                                        cutHeight = cut.height,
                                        method = "hybrid",
                                        distM = dis,
                                        deepSplit = deep.split,
                                        pamRespectsDendro = FALSE,
                                        minClusterSize = min.size)
    list(adj = adj, dis = dis, dendro = dendro, index.vector = iv,
         n.mods = length(unique(iv[iv != 0])), beta = b)
  }

  message("Generating initial modules...")
  fit <- cut.at.power(beta)
  n.powers.tried <- 1L

  # RETRY ACROSS POWERS. pickSoftThreshold selects the SMALLEST power whose
  # signed scale-free R^2 clears min.sft, which says nothing about whether the
  # resulting TOM actually resolves into modules -- and on a weakly scale-free
  # network that crossing point is noise-dominated. When the natural cut comes
  # back with fewer than min.natural.mods modules, walk the remaining powers
  # (best scale-free fit first) and keep the first that does resolve.
  #
  # This matters because the degenerate case is silent: a cut that assigns EVERY
  # node to the unassigned group (0) cannot be repaired downstream --
  # .assign_unassigned has no module to recruit into, and .split_to_max_size
  # deliberately never splits group 0 -- so the whole graph would be returned as
  # a single "module 0". Retrying, then failing loudly, replaces that.
  if (fit$n.mods < min.natural.mods && auto.beta && retry.betas) {
    r2 <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
    candidates <- sft$fitIndices$Power[order(r2, decreasing = TRUE)]
    candidates <- setdiff(candidates, beta)
    for (b in candidates) {
      message("Natural cut gave ", fit$n.mods, " module(s) at power ", fit$beta,
              "; retrying at power ", b, ".")
      trial <- cut.at.power(b)
      n.powers.tried <- n.powers.tried + 1L
      if (trial$n.mods >= min.natural.mods) { fit <- trial; break }
      fit <- trial   # keep the most recent attempt so the error below is honest
    }
  }

  if (fit$n.mods < 1) {
    stop("WGCNA module detection failed: the dynamic tree cut assigned every ",
         "feature to the unassigned group ",
         if (n.powers.tried > 1)
           paste0("at all ", n.powers.tried, " soft-threshold powers tried")
         else if (auto.beta && !retry.betas)
           paste0("at power ", fit$beta, " (retry.betas = FALSE; set it TRUE to try other powers)")
         else
           paste0("at the supplied power ", fit$beta, " (pass beta = NULL to search other powers)"),
         ". Use find_ICA_mods(), supply a partition directly, or check the input ",
         "for degenerate correlation structure.", call. = FALSE)
  }
  if (fit$n.mods < min.natural.mods) {
    warning("Natural module detection produced only ", fit$n.mods,
            " module(s) (min.natural.mods = ", min.natural.mods,
            "); the partition may be degenerate.", call. = FALSE)
  }

  beta   <- fit$beta
  adj    <- fit$adj
  dis    <- fit$dis
  dendro <- fit$dendro
  initial.index.vector <- fit$index.vector

  # merge similar modules based on eigengene similarity
  if (merge) {
    initial.index.vector <- .merge_modules(index.vector = initial.index.vector,
                                           t.x = t.x,
                                           cor.FN = cor.FN,
                                           cor.options = cor.options,
                                           merging.cut = merging.cut,
                                           min.mods = min.mods,
                                           max.size = max.size)
  }

  # assign remaining unassigned nodes to a module
  message("Assigning unassigned nodes by ", assign.by, " method...")
  diag(adj) <- 0
  initial.index.vector <- .assign_unassigned(index.vector = initial.index.vector,
                                             adj = adj,
                                             x = x,
                                             method = assign.by)

  # enforce max.size by recursively splitting oversized modules. A cutree split
  # of an oversized module's own sub-dendrogram can always reach max.size without
  # evicting nodes into a foreign module (unlike node trading), so we never break
  # up genuinely connected features to satisfy the cap. Any pieces left below
  # min.size are then merged back into their nearest neighbouring module.
  final.index.vector <- initial.index.vector
  if (iterate) {
    final.index.vector <- .split_to_max_size(index.vector = final.index.vector,
                                             dis = dis,
                                             max.size = max.size,
                                             min.size = min.size,
                                             hclust.method = hclust.method)
    final.index.vector <- .merge_small_modules(index.vector = final.index.vector,
                                               dis = dis,
                                               min.size = min.size,
                                               max.size = max.size)
  } else {
    n.too.big <- sum(table(final.index.vector) > max.size)
    if (n.too.big > 0) {
      warning(n.too.big, " module(s) exceed max.size; set iterate = TRUE to split them.")
    }
  }

  # renumber the final modules to a contiguous 1:N sequence. Iterative splitting
  # (.split_to_max_size) hands out inflated ids and small-module merging
  # (.merge_small_modules) deletes ids, so the labels arrive full of gaps
  # (e.g. 1, 45, 101, 4); relabel so the returned modules read 1:N. The
  # unassigned group (0) is preserved and the labels are deliberately not tied
  # to the initial module ids.
  final.index.vector <- .relabel_sequential(final.index.vector)

  # The unassigned group (0) is NOT a module. split() below would happily emit it
  # as one, which is how a total detection collapse used to be returned as a
  # single valid-looking module covering the whole graph. The retry above makes
  # that path unreachable via cutreeDynamic; this asserts the contract for every
  # other route (an explicit `beta`, retry.betas = FALSE, or a node that scores
  # zero against every module in .assign_unassigned).
  leftover <- sum(final.index.vector == 0)
  if (leftover > 0) {
    stop(leftover, " feature(s) could not be assigned to any module and would be ",
         "returned as a spurious 'module 0'. A module object must partition every ",
         "feature into a real module.", call. = FALSE)
  }

  # build module objects (a non-overlapping partition owns all of its own nodes,
  # so core.list == index.list)
  initial.index.vector <- as.numeric(initial.index.vector)
  final.index.vector   <- as.numeric(final.index.vector)
  build_mods <- function(src, iv){
    il <- split(seq_len(nrow(x)), iv)
    methods::new("module",
                 source = src,
                 data.dim = dim(x),
                 overlapping = FALSE,
                 index.vector = iv,
                 index.list = il,
                 name.list = split(rownames(x), iv),
                 core.list = il)
  }
  initial.mods <- build_mods("find_WGCNA_mods initial mods", initial.index.vector)
  final.mods   <- build_mods("find_WGCNA_mods sized mods",   final.index.vector)

  return(list(
    wgcna.adj = adj,
    initial.mods = initial.mods,
    final.mods = final.mods
  ))
}

#' Count the modules in an index vector, excluding the unassigned group (0)
#' @param index.vector an integer vector assigning each node to a module (0 = unassigned)

#' @return an integer, the number of distinct non-zero modules

#' @keywords internal
.n_modules <- function(index.vector){
  length(unique(index.vector[index.vector != 0]))
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

#' Helper to find_WGCNA_mods that enforces max.size by recursively splitting
#' oversized modules on their own sub-dendrograms
#'
#' Unlike node trading, a cutree split can always reach max.size (in the limit,
#' singletons) without moving a node into an unrelated module, so genuinely
#' connected features are never separated to satisfy the cap.
#' @param index.vector an integer vector of length p assigning each node to a module (0 = unassigned)
#' @param dis a p x p TOM distance matrix
#' @param max.size a numeric, the maximum number of nodes allowed in a module
#' @param min.size an integer, the minimum module size (informational; small pieces are reconciled by .merge_small_modules)
#' @param hclust.method a character, the agglomeration method passed to flashClust

#' @return an integer vector of length p with every non-zero module within max.size

#' @importFrom stats as.dist cutree
#' @keywords internal
.split_to_max_size <- function(index.vector,
                               dis,
                               max.size,
                               min.size,
                               hclust.method = "average") {
  if (is.infinite(max.size)) return(index.vector)
  repeat {
    sizes <- table(index.vector)
    big <- as.numeric(names(sizes[sizes > max.size]))
    big <- big[big != 0]                       # never split the unassigned group
    if (length(big) == 0) break

    m <- big[1]                                # split one module, then re-evaluate
    m.nodes <- which(index.vector == m)
    dend <- flashClust::flashClust(stats::as.dist(dis[m.nodes, m.nodes, drop = FALSE]),
                                   method = hclust.method)

    # smallest number of cuts that brings every piece within max.size
    k <- ceiling(length(m.nodes) / max.size)
    repeat {
      k <- min(k, length(m.nodes))
      sub <- stats::cutree(dend, k = k)
      if (max(table(sub)) <= max.size || k >= length(m.nodes)) break
      k <- k + 1
    }
    message("Splitting modules iteratively to satisfy max.size.")

    # first piece keeps the id m, the others get fresh ids beyond the current max
    offset <- max(index.vector)
    index.vector[m.nodes] <- ifelse(sub == 1L, m, offset + sub - 1L)
  }
  index.vector
}

#' Helper to find_WGCNA_mods that merges modules below min.size into their
#' nearest neighbouring module (by mean TOM distance) when there is room
#' @param index.vector an integer vector of length p assigning each node to a module (0 = unassigned)
#' @param dis a p x p TOM distance matrix
#' @param min.size an integer, the minimum number of nodes allowed in a module
#' @param max.size a numeric, the maximum module size a merge target may reach

#' @return an integer vector of length p with undersized modules merged where possible

#' @keywords internal
.merge_small_modules <- function(index.vector, dis, min.size, max.size) {
  repeat {
    sizes <- table(index.vector)
    small <- as.numeric(names(sizes[sizes < min.size]))
    small <- small[small != 0]
    if (length(small) == 0) break

    m <- small[which.min(sizes[as.character(small)])]   # smallest undersized module first
    m.nodes <- which(index.vector == m)
    others <- setdiff(unique(index.vector[index.vector != 0]), m)
    if (length(others) == 0) break

    # nearest other module by mean TOM distance that has room to absorb m
    mean.d <- vapply(others, function(o) mean(dis[m.nodes, which(index.vector == o)]), numeric(1))
    room   <- vapply(others, function(o) sum(index.vector == o) + length(m.nodes) <= max.size, logical(1))
    if (!any(room)) break                                # nowhere to merge without exceeding max.size

    target <- others[room][which.min(mean.d[room])]
    index.vector[m.nodes] <- target
  }
  index.vector
}

#' Helper to find_WGCNA_mods that renumbers module labels to a contiguous 1:N
#' sequence, closing the gaps and inflated ids left behind by iterative
#' splitting and small-module merging. The unassigned group (0) is left as-is.
#' @param index.vector an integer vector of length p assigning each node to a module (0 = unassigned)

#' @return an integer vector of length p whose non-zero modules are relabelled 1..N (in ascending order of their original ids)

#' @keywords internal
.relabel_sequential <- function(index.vector) {
  assigned <- index.vector != 0
  # map the distinct non-zero ids (ascending) onto 1..N, leaving 0 untouched
  labels <- sort(unique(index.vector[assigned]))
  index.vector[assigned] <- match(index.vector[assigned], labels)
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
  stalled <- FALSE
  while (any(index.vector == 0) && any(index.vector != 0) && !stalled) {
    unassigned <- which(index.vector == 0)
    assigned <- which(index.vector != 0)

    if (method == "adjacency") {
      # for each unassigned node, find its single best-connected assigned node
      score <- adj[unassigned, assigned, drop = FALSE]
      best.col <- max.col(score, ties.method = "first")
      best.score <- score[cbind(seq_along(unassigned), best.col)]
      best.module <- index.vector[assigned[best.col]]
    } else {
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
    if (any(recruit)) {
      index.vector[unassigned[recruit]] <- best.module[recruit]
      if (verbose) {
        message("Assigned ", sum(recruit), " node(s) at ", method, " threshold ", threshold)
      }
    } else if (threshold > 0) {
      # nothing clears the threshold; lower it
      threshold <- round(threshold - step, 2)
    } else {
      # threshold exhausted and nothing was recruited; stop
      stalled <- TRUE
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
#' @param t.x a samples x features numeric matrix (i.e. t(x))
#' @param cor.FN a character, the correlation function ('bicor' or 'cor')
#' @param cor.options a list of correlation options passed to mergeCloseModules
#' @param merging.cut a numeric, the eigengene dissimilarity threshold for merging
#' @param min.mods an integer, the minimum number of modules to retain
#' @param max.size a numeric, modules are only merged while all are below this size

#' @return an integer vector of length p with merged module assignments

#' @importFrom WGCNA mergeCloseModules
#' @keywords internal
.merge_modules <- function(index.vector, t.x, cor.FN, cor.options, merging.cut, min.mods, max.size) {
  message("Merging modules based on eigengene similarity...")
  n.mods <- .n_modules(index.vector)
  n.merges <- 0
  module.eigengenes <- NULL # reused across rounds to skip recomputing module eigengenes
  while (n.mods > min.mods && all(table(index.vector) < max.size)) {
    merged <- WGCNA::mergeCloseModules(
      exprData = t.x,
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
    n.mods <- .n_modules(index.vector)
    n.merges <- n.merges + 1
  }
  if (n.merges == 0) {
    message("No modules were small/similar enough to merge; using the initial modules.")
  } else {
    message("Merged modules in ", n.merges, " round(s), resulting in ", n.mods, " modules.")
  }
  return(index.vector)
}

#' Detect overlapping co-expression modules from a data matrix using Independent Component Analysis (ICA)
#'
#' ICA decomposes the features into \code{n.comp} independent components
#' (metagenes). Each feature is \emph{owned} by (is a core node of) its dominant
#' component -- the component in which it has the largest absolute loading -- so
#' the cores form a clean partition. A feature is additionally recruited as an
#' \emph{auxiliary} member of every other component in which its standardised
#' loading exceeds \code{membership.z}; this multi-membership is what makes the
#' returned modules overlap, so ICA yields divide-and-conquer-ready fuzzy modules
#' directly, without a separate growth step. \code{max.size} trims auxiliary
#' members (weakest loadings first) but never core members, so a component with
#' more than \code{max.size} cores signals that \code{n.comp} should be raised.
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param n.comp an integer, the number of independent components (modules) to extract
#' @param max.size an integer or NULL, the maximum number of nodes allowed in a module (NULL = unconstrained); trims auxiliary members only
#' @param membership.z a numeric, the standardised absolute-loading threshold above which a feature becomes an auxiliary member of a non-dominant component
#' @param ... additional arguments passed to fastICA::fastICA

#' @return a list with 'ica.loadings', the p x n.comp source (metagene) loading matrix, and 'mods', an overlapping module object whose core nodes are each feature's dominant component

#' @importFrom methods new

#' @export
find_ICA_mods <- function(x,
                          n.comp,
                          max.size = NULL,
                          membership.z = 2,
                          ...) {
  if (!requireNamespace("fastICA", quietly = TRUE)) {
    stop("Package fastICA is required. Install with: install.packages('fastICA')", call. = FALSE)
  }
  if (is.null(max.size)) max.size <- Inf

  # decompose features into independent components; S has one row per feature and
  # one column per component (the metagene loadings)
  ICA.results <- fastICA::fastICA(X = as.matrix(x), n.comp = n.comp, ...)
  S <- ICA.results$S
  rownames(S) <- rownames(x)
  abs.S <- abs(S)

  # ownership: each feature's core (dominant) component is its top absolute loading
  core.assign  <- max.col(abs.S, ties.method = "first")
  score.vector <- abs.S[cbind(seq_len(nrow(abs.S)), core.assign)]

  # multi-membership: a feature is an auxiliary member of any other component in
  # which its standardised loading clears membership.z (its own core is always in)
  z <- scale(abs.S)
  member <- z >= membership.z
  member[is.na(member)] <- FALSE
  member[cbind(seq_len(nrow(member)), core.assign)] <- TRUE
  message("ICA produced ", n.comp, " components.")

  # assemble per-component index (core + auxiliary) and core lists, trimming
  # auxiliary members by loading to respect max.size (cores are never trimmed)
  index.list <- vector("list", n.comp)
  core.list  <- vector("list", n.comp)
  for (cc in seq_len(n.comp)) {
    cores   <- which(core.assign == cc)
    members <- which(member[, cc])
    if (length(members) > max.size) {
      if (length(cores) > max.size) {
        warning("Component ", cc, " has ", length(cores),
                " core features (> max.size); raise n.comp to reduce module size.")
      }
      aux       <- setdiff(members, cores)
      keep.n    <- max(0, max.size - length(cores))
      aux.keep  <- aux[order(abs.S[aux, cc], decreasing = TRUE)][seq_len(min(keep.n, length(aux)))]
      members   <- c(cores, aux.keep)
    }
    index.list[[cc]] <- sort(members)
    core.list[[cc]]  <- sort(cores)
  }

  # drop components that own no feature (their core partition entry is empty)
  keep <- lengths(core.list) > 0
  index.list <- index.list[keep]
  core.list  <- core.list[keep]

  mods <- methods::new("module",
                       source = "find_ICA_mods",
                       data.dim = dim(x),
                       overlapping = TRUE,
                       score.vector = score.vector,
                       index.list = index.list,
                       name.list = lapply(index.list, function(m) rownames(x)[m]),
                       core.list = core.list)

  return(list(
    ica.loadings = S,
    mods = mods
  ))
}


###############################
### Expand to Fuzzy Modules ###
###############################

#' Expand non-overlapping modules to fuzzy (overlapping) modules by recruiting nodes correlated with each module's eigengene
#'
#' The nodes of the input partition become the \emph{core} (owned) nodes of each
#' fuzzy module; the recruited neighbours are auxiliary nodes, present only so the
#' graph learned on the module can condition on each core node's Markov blanket.
#' Growth stops at \code{ratio}/\code{max.size}, or earlier if \code{min.cor} is
#' supplied and no external feature clears that correlation (an approximate
#' blanket-closure criterion).
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param max.size an integer, the maximum number of nodes allowed in any fuzzy module
#' @param n.pc an integer, the number of principal components used to represent each module's eigengene
#' @param ratio a numeric, the maximum ratio of fuzzy module size to original module size
#' @param min.cor a numeric or NULL; if set, only features whose summed absolute correlation with the module eigengene exceeds this value are recruited (blanket-closure cutoff), so a module may grow by fewer than \code{ratio} implies

#' @return a module object

#' @importFrom stats prcomp cor
#' @importFrom methods new

#' @export
eigen_fuzzy_modules <- function(x, input.modules, max.size, n.pc = 2, ratio = 1.5, min.cor = NULL){
  # growth only defines ownership cleanly when the input is a partition
  if (input.modules@overlapping) {
    stop("input.modules must be a non-overlapping partition (its nodes become the core nodes of each fuzzy module).")
  }
  # check if any modules are too large
  if (any(lengths(input.modules@index.list) > max.size)) {
    stop("Some modules are too large, increase max.size.")
  }

  # create new index list (list of nodes in each fuzzy module by index)
  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod.nodes <- input.modules@index.list[[m]]
    # get number of required fuzzy nodes (capped by max.size and by nodes available)
    f.size <- length(mod.nodes) * ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- min(f.size - length(mod.nodes), nrow(x) - length(mod.nodes))
    if(n.fuzzy.nodes <= 0){ return(sort(mod.nodes)) }

    # get the module's principal components
    mod.PC <- stats::prcomp(t(x[mod.nodes, , drop = FALSE]), scale. = TRUE)

    # score each gene outside the module by its summed absolute correlation
    # with the module's first n.pc principal components
    n.pc.use <- min(n.pc, ncol(mod.PC$x))
    pc.cor <- stats::cor(t(x[-mod.nodes, , drop = FALSE]),
                         mod.PC$x[, seq_len(n.pc.use), drop = FALSE])
    eigen.cor <- rowSums(abs(pc.cor))

    corRank <- sort(eigen.cor, decreasing = TRUE)      # ranked absolute correlation
    corRank <- corRank[seq_len(n.fuzzy.nodes)]         # highest-ranked candidates
    if(!is.null(min.cor)){ corRank <- corRank[corRank >= min.cor] } # blanket-closure cutoff
    fuzzy.nodes <- which(rownames(x) %in% names(corRank))

    # return fuzzy module combining original (core) and recruited (auxiliary) nodes
    sort(c(mod.nodes, fuzzy.nodes))
  })

  # convert to module object and return; original nodes are the owned core nodes
  fuzzy.mods <- methods::new("module",
                             source = paste("eigen_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(m){rownames(x)[m]}),
                             core.list = input.modules@index.list
  )
  return(fuzzy.mods)

}

#' Expand non-overlapping modules to fuzzy (overlapping) modules using a thresholded WGCNA adjacency matrix
#'
#' The nodes of the input partition become the \emph{core} (owned) nodes of each
#' fuzzy module; the recruited neighbours are auxiliary nodes carried along so
#' the module can condition on each core node's Markov blanket. Growth stops at
#' \code{ratio}/\code{max.size}, or earlier if \code{min.adj} is supplied and no
#' external feature reaches that adjacency (an approximate blanket-closure cutoff).
#' @param x a numeric matrix with p features (rows) and n samples (columns)
#' @param adj a p x p numeric adjacency matrix (e.g. from WGCNA::adjacency) defining pairwise feature similarity
#' @param input.modules a module S4 object containing the non-overlapping module assignments to expand
#' @param max.size an integer, the maximum number of nodes allowed in any fuzzy module
#' @param ratio a numeric, the maximum ratio of fuzzy module size to original module size
#' @param min.adj a numeric or NULL; if set, only features whose maximum adjacency to a core node exceeds this value are recruited (blanket-closure cutoff)

#' @return a module object

#' @importFrom methods new

#' @export
adj_fuzzy_modules <- function(x, adj, input.modules, max.size, ratio, min.adj = NULL){
  # growth only defines ownership cleanly when the input is a partition
  if (input.modules@overlapping) {
    stop("input.modules must be a non-overlapping partition (its nodes become the core nodes of each fuzzy module).")
  }
  # check if any modules are too large
  if (any(lengths(input.modules@index.list) > max.size)) {
    stop("Some modules are too large, increase max.size.")
  }

  # create new index list (list of nodes in each fuzzy module by index)
  index.list <- lapply(seq_along(input.modules@index.list), function(m){
    mod.nodes <- input.modules@index.list[[m]]
    # get number of required fuzzy nodes (capped by max.size and by nodes available)
    f.size <- length(mod.nodes) * ratio
    if(f.size > max.size){f.size <- max.size}
    n.fuzzy.nodes <- min(f.size - length(mod.nodes), nrow(x) - length(mod.nodes))
    if(n.fuzzy.nodes <= 0){ return(sort(mod.nodes)) }

    # get matrix of nodes in module adj with nodes outside module
    in.out.adj <- adj[mod.nodes, -mod.nodes, drop = FALSE]

    # get max adj of node outside module with a node inside module
    out.max <- apply(in.out.adj, 2, max)

    # rank and select fuzzy nodes based on adjacency to the module
    adj.rank <- order(out.max, decreasing = TRUE)[seq_len(n.fuzzy.nodes)]
    if(!is.null(min.adj)){ adj.rank <- adj.rank[out.max[adj.rank] >= min.adj] } # closure cutoff
    fuzzy.nodes <- which(rownames(x) %in% colnames(in.out.adj)[adj.rank])

    # return fuzzy module combining original (core) and recruited (auxiliary) nodes
    sort(c(mod.nodes, fuzzy.nodes))
  })

  # convert to module object and return; original nodes are the owned core nodes
  fuzzy.mods <- methods::new("module",
                             source = paste("adj_fuzzy_modules", "generated from", input.modules@source),
                             data.dim = dim(x),
                             overlapping = TRUE,
                             index.list = index.list,
                             name.list = lapply(index.list, function(m){rownames(x)[m]}),
                             core.list = input.modules@index.list
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
  # membership indicators over the shared node universe; the matrix product then
  # counts, for each (i, j), the nodes shared between module i and module j
  membership <- function(idx){
    matrix(vapply(idx, function(s) as.integer(nodes1 %in% s), integer(length(nodes1))),
           nrow = length(nodes1))
  }
  overlap.matrix <- crossprod(membership(idx1), membership(idx2)) # t(A) %*% B
  dimnames(overlap.matrix) <- list(lab1, lab2)

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
  counts <- vapply(module.names, function(nm){
    # membership indicator over the graph's vertices
    member <- v.names %in% nm
    if(nrow(edges) == 0){
      return(c(within = 0, between = 0))
    }
    a.in <- member[edges[, 1]]
    b.in <- member[edges[, 2]]
    c(within  = sum(a.in & b.in),       # both endpoints inside the module
      between = sum(xor(a.in, b.in)))   # exactly one endpoint inside the module
  }, numeric(2))

  within.edges   <- counts["within", ]
  between.edges  <- counts["between", ]
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
    cm <- cor.list[[i]]
    vals <- cm[upper.tri(cm)]
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
