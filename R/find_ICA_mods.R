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
