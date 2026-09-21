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
