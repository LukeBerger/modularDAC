#' Calculate the F1 score comparing a predicted graph against a true graph
#' @param g.true an igraph object, the true reference graph
#' @param g.pred an igraph object, the predicted graph to evaluate

#' @return a list of F1 score and its components

#' @importFrom igraph V as_edgelist

#' @export
calc_F1 <- function(g.true, g.pred) {
  # check that the graphs have the same nodes
  if (!all(igraph::V(g.true)$name %in% igraph::V(g.pred)$name)) {
    stop("True and Predicted graphs must contain the same nodes")
  }

  # extract edges as character vectors in a consistent format
  edges.true <- apply(igraph::as_edgelist(g.true), 1, function(x) paste(sort(x), collapse = " -- "))
  edges.pred <- apply(igraph::as_edgelist(g.pred), 1, function(x) paste(sort(x), collapse = " -- "))

  # count True Positives, False Positives, False Negatives
  tp <- sum(edges.pred %in% edges.true)
  fp <- sum(!(edges.pred %in% edges.true))
  fn <- sum(!(edges.true %in% edges.pred))

  # calc precision and recall
  precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))

  # calc F1 Score
  if (precision + recall == 0) {
    f1 <- 0 # to avoid zero divison error in edge case
  } else {
    f1 <- 2 * precision * recall / (precision + recall)
  }

  # return values as list
  return(list(
    F1 = f1,
    Precision = precision,
    Recall = recall,
    TruePos = tp,
    FalsePos = fp,
    FalseNeg = fn
  ))
}
