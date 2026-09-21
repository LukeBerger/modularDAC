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
