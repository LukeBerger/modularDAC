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
#' @param deep.split an integer 0-4 passed to \code{dynamicTreeCut::cutreeDynamic}, controlling how aggressively the dendrogram is split; 4 (the default) splits most aggressively, 0 least
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
    # under the same correlation and network settings the adjacency below is
    # built with; the two functions' defaults differ, which would otherwise fit
    # the power to a different matrix than the one it is applied to.
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

  # Retry across powers. pickSoftThreshold selects the smallest power whose
  # signed scale-free R^2 clears min.sft, which says nothing about whether the
  # resulting TOM resolves into modules, and on a weakly scale-free network that
  # crossing point is noise-dominated. When the natural cut returns fewer than
  # min.natural.mods modules, walk the remaining powers (best scale-free fit
  # first) and keep the first that does resolve.
  #
  # The degenerate case is otherwise silent: a cut that assigns every node to
  # the unassigned group (0) cannot be repaired downstream, since
  # .assign_unassigned has no module to recruit into and .split_to_max_size
  # never splits group 0, so the whole graph would be returned as a single
  # "module 0".
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
      fit <- trial   # keep the most recent attempt for the error message below
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

  # enforce max.size by recursively splitting oversized modules on their own
  # sub-dendrograms, then merging any piece left below min.size back into its
  # nearest neighbouring module
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

  # renumber the final modules to a contiguous 1:N sequence. Splitting
  # (.split_to_max_size) hands out inflated ids and small-module merging
  # (.merge_small_modules) deletes ids, so the labels arrive with gaps
  # (e.g. 1, 45, 101, 4). The unassigned group (0) is preserved, and the new
  # labels are not tied to the initial module ids.
  final.index.vector <- .relabel_sequential(final.index.vector)

  # The unassigned group (0) is not a module: split() below would otherwise emit
  # it as one, returning a failed detection as a single module covering the whole
  # graph. The retry above rules that out for the automatic path; this asserts
  # the contract for the others (an explicit `beta`, retry.betas = FALSE, or a
  # node that scores zero against every module in .assign_unassigned).
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

#' Helper to find_WGCNA_mods that enforces max.size by recursively splitting
#' oversized modules on their own sub-dendrograms
#'
#' Cutting a module's own sub-dendrogram always reaches max.size (in the limit,
#' singletons) without moving a node into an unrelated module, so no node is
#' evicted into a foreign module to satisfy the cap.
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
