# library(fastICA)
# library(DescTools)
# library(dplyr)


## Complex Module detection via multistage algorithm (based on ICA core modules)

## TAKES: data matrix (transformed eset exprs), number of components to look for,
## number of nodes in "core modules", intial threshold for determining neighbors
## minium threshold for determining neighbors

## RETURNS: list of module assignments

complex_ICA_modules <- function(x,n.comp, core.size = 5, starting.threshold = 0.5, min.threshold = 0){
  #########################
  ###Define core modules###
  #########################

  #preform ICA
  ICA.results <- fastICA::fastICA(X= x,
                    n.comp = n.comp,
                    verbose = 0
  )
  s.matrix <- abs(ICA.results$S)
  #get correlation matrix from data
  cor.matrix <- stats::cor(t(x))
  diag(cor.matrix) <- 0 # set diagonals as zero to ignore loops

  #get basic mods based on max component score of ICA
  basic.mods <- apply(s.matrix, 1, which.max) #modules named based on ICA column n

  #very stringent modules based on top 5 values
  complex.mods <- basic.mods #probably a better way to this, may not need to define basic mods
  top5 <- apply(s.matrix, 2, function(col){
    which(col %in% sort(col, decreasing = T)[1:core.size])
  })
  complex.mods[-top5] <- 0

  #########################
  ###Perform Assignments###
  #########################

  ### Assignments based on NEIGHBORS with decreasing threshold
  #while there are unassigned nodes
  last.mods <- complex.mods
  threshold = starting.threshold
  min.threshold = 0
  unassigned <- which(complex.mods == 0)
  while(length(unassigned) > 0){ #there might be a more efficient way to check this\
    unassigned <- which(complex.mods == 0)
    assigned <- which(complex.mods != 0)
    #for each unassigned node
    for(node in unassigned){
      #get neighbors based on current threshold
      nbs <- assigned[which(abs(cor.matrix[node, assigned]) > threshold)]
      #and assign this node based on the mode of their modules
      md <- DescTools::Mode(last.mods[nbs])[1]
      if(length(md) == 1){
        if(!is.na(md)){
          complex.mods[node] <- md
        }
      }
    }

    #if no changes
    if(all(last.mods == complex.mods)){
      #and not at min threshold
      if(threshold > min.threshold){
        #lower threshold
        threshold <- threshold - 0.05
      }else{
        #break if at min threshold and no longer changing
        break
      }
    }
    last.mods <- complex.mods
  }

  # convert to module object and return
  cICA.mods <- new("module",
                     source = "complexICA",
                     data.dim = dim(x),
                     overlapping = FALSE,
                     index.vector = complex.mods,
                     score.vector = apply(abs(ICA.results$S), 1, max),
                     index.list =split(1:nrow(x) , complex.mods),
                     name.list = split(rownames(x), complex.mods)
                     )

  return(cICA.mods)
}

##########
###TEST###
##########

source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/00.SimulateGraphs.R")
source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/modularDAC/R/01.DetectModules.R")

er <- make_modular_graph()
x <- sim_graph_data(er, n.samples = 100)
mod <- complex_ICA_modules(x, 3)

#############
###TESTING###
#############

#check feature number matches input data
if(length(mod@index.vector) != nrow(x)){
  stop(paste(mod@source, "produced a index vector with the incorrect number of features"))
}
#check that all feature names in modules come from the data
if(!all(unlist(mod@name.list) %in% rownames(x))){
  stop(paste(mod@source, "feature names do not match input data"))
}
#check that each module has the same number of feature indexes and names
if(
  !all(
    unlist(
      lapply(seq_along(mod@index.list), function(i) length(mod@index.list[[i]]) == length(mod@name.list[[i]]))
      )
    )
){stop(paste(mod@source, "produced differnet length index and name lists"))}






##################################################
###Alternative methods for defining core module###
##################################################

###Define core modules bases on some threshold of relative component score
# complex.mods[
#   which(!apply(s.matrix, 1, function(row) {
#     mxCol <- which.max(row)
#     #get assignments with a dominant s score and high overall values
#     (row[mxCol] > sum(row[-mxCol]*coreThreshold)) &&  sum(row) > 1
#   }
#   ))
# ] <- 0



#################################################
###Alternative methods for initial assignments###
#################################################

### Assignments based on avg covariance of the module with decreasing threshold
# #while there are unassigned nodes
# last.mods <- complex.mods
# threshold = 0.95 #will make proper args in function later
# min.threshold = 0
# unassigned <- which(complex.mods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each unassigned node
#   for(node in unassigned){
#     #get avg covariance with all current assigned nodes
#     covars <- abs(cor.matrix[node, ])
#     modCovars <- split(covars,last.mods)
#     modCovars["0"] <- NULL #remove unassigned
#     avgCovars <-lapply(modCovars, mean)
#
#     #assign based on max covariance and ICA score
#     maxCov <- max(unlist(avgCovars))
#     wchMax <- which.max(avgCovars)
#     if(maxCov > threshold ){
#       if(useICA){
#         if(which.max(s.matrix[node,]) == wchMax){
#           complex.mods[node]  <- wchMax #might have issues with name matching
#         }
#       }else{
#         complex.mods[node]  <- wchMax
#       }
#
#     }
#
#   }
#   #update based on new assignments
#   unassigned <- which(complex.mods == 0)
#
#   #if no changes
#   if(all(last.mods == complex.mods)){
#     #and not at min threshold
#     if(threshold > min.threshold){
#       #lower threshold
#       print("lowering threshold")
#       threshold <- threshold - 0.05
#     }else{
#       #break if at min threshold and no longer changing
#       break
#     }
#   }
#   last.mods <- complex.mods
#   print(threshold)
# }

### Assignments based on closest neighbor
#while there are unassigned nodes
# last.mods <- complex.mods
# unassigned <- which(complex.mods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each assigned node
#   assigned <- which(complex.mods != 0)
#   unassigned <- which(complex.mods == 0)
#   for(node in assigned){
#     #assign its "best" unassigned neighbor (based on covar) the same module
#     best <- unassigned[which.max(abs(cor.matrix[node, unassigned]))]
#     complex.mods[best] <- complex.mods[node]
#   }
#   #end iterations if stuck
#   if(all(last.mods == complex.mods)){break}
#   last.mods <- complex.mods
# }

### Assignments based on closeness to entire module
#while there are unassigned nodes
# last.mods <- complex.mods
# unassigned <- which(complex.mods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each module
#   unassigned <- which(complex.mods == 0)
#   for(module in 1:n.comp){ #what module goes last can steal nodes, need to store values then take max
#     subCov <- cor.matrix[which(complex.mods == module), unassigned]
#     #find the unassigned node closest to the overall module
#     complex.mods[unassigned[which.max(colMeans(abs(subCov)))]] <- module
#
#   }
#   #end iterations if stuck
#   if(all(last.mods == complex.mods)){break}
#   last.mods <- complex.mods
# }

###############################################
###Methods for correcting flawed assignments###
###############################################

### Reassign border nodes based on second degree neighbors
# threshold <- 0.2
# last.mods <- complex.mods
#
# #border nodes are nodes with neighbors of a different module
# borderNodes <- which(unlist(
#   lapply(complex.mods, function(node){
#     all(unique(complex.mods[which(abs(cor.matrix[node, ]) > threshold)]) == complex.mods[node])
#   })
# ))
# while(length(borderNodes > 0)){
#   for(node in borderNodes){
#     #get neighbors (based on threshold) and neighbors neighbors
#     n1 <- which(abs(cor.matrix[node, ]) > threshold)
#     n2 <- unique(unlist(lapply(n1, function(n){
#         which(abs(cor.matrix[n, ]) > threshold)
#       })))[-node]
#     #assign node to mode of n2 modules
#     md <- DescTools::Mode(last.mods[n2])[1]
#     if(length(md) == 1){
#       if(!is.na(md)){
#         complex.mods[node] <- md
#       }
#     }
#   }
#   if(all(last.mods == complex.mods)){break}
#   last.mods <- complex.mods
#   borderNodes <- which(unlist(
#     lapply(complex.mods, function(node){
#       all(unique(complex.mods[which(abs(cor.matrix[node, ]) > threshold)]) == complex.mods[node])
#     })
#   ))
# }

###Reassign border nodes based on covariance with final modules
# #border nodes are those whose closest neighbor is of a different module
# borderNodes <- which(complex.mods != complex.mods[apply(abs(cor.matrix),1, which.max)])
# last.mods <- complex.mods
# while(length(borderNodes > 0)){
#   for(node in borderNodes){
#     #get avg covariance with all modules
#     covars <- abs(cor.matrix[node, ])
#     modCovars <- split(covars,last.mods)
#     modCovars["0"] <- NULL #remove unassigned (if any still exist)
#     avgCovars <-lapply(modCovars, mean)
#
#     #assign based on max covariance
#     complex.mods[node] <- which.max(avgCovars)
#   }
#   if(all(last.mods == complex.mods)){break}
#   last.mods <- complex.mods
#   borderNodes <- which(complex.mods != complex.mods[apply(abs(cor.matrix),1, which.max)])
# }

