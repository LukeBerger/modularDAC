# silence R CMD check NOTE for ggplot2 aes() non-standard evaluation in module_correlation
utils::globalVariables(c("module", "correlation"))

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
