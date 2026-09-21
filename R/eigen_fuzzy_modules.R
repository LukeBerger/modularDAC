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
