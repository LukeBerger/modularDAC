# library(fastICA)
# library(DescTools)
library(dplyr)


## Complex Module detection via multistage algorithm (based on ICA core modules)

## TAKES: data matrix (transformed eset exprs), number of components to look for,
## number of nodes in "core modules", intial threshold for determining neighbors
## minium threshold for determining neighbors

## RETURNS: list of module assignments

complexICA <- function(dataMat,nComp, coreSize = 5, intialThreshold = 0.5, minThreshold = 0){
  # set.seed(123) #setting seed removes randomness in fastaICA iterations

  #########################
  ###Define core modules###
  #########################

  #preform ICA
  ICAOut <- fastICA::fastICA(X= dataMat,
                    n.comp = nComp,
                    verbose = 0
  )
  sMat <- abs(ICAOut$S)
  #get covariance matrix from data
  covMat <- cov(t(dataMat))
  diag(covMat) <- 0 # set diagonals as zero to ignore loops

  #get basic mods based on max component score of ICA
  basicMods <- apply(sMat, 1, which.max) #modules named based on ICA column n

  #very stringent modules based on top 5 values
  complexMods <- basicMods #probably a better way to this, may not need to define basic mods
  top5 <- apply(sMat, 2, function(col){
    which(col %in% sort(col, decreasing = T)[1:coreSize])
  })
  complexMods[-top5] <- 0

  #########################
  ###Intial Asssignments###
  #########################

  ### Assignments based on NEIGHBORS with decreasing threshold
  #while there are unassigned nodes
  lastMods <- complexMods
  threshold = intialThreshold
  minThreshold = 0
  unassigned <- which(complexMods == 0)
  while(length(unassigned) > 0){ #there might be a more efficient way to check this\
    unassigned <- which(complexMods == 0)
    assigned <- which(complexMods != 0)
    #for each unassigned node
    for(node in unassigned){
      #get neighbors based on current threshold
      nbs <- assigned[which(abs(covMat[node, assigned]) > threshold)]
      #and assign this node based on the mode of their modules
      md <- DescTools::Mode(lastMods[nbs])[1]
      if(length(md) == 1){
        if(!is.na(md)){
          complexMods[node] <- md
        }
      }
    }

    #if no changes
    if(all(lastMods == complexMods)){
      #and not at min threshold
      if(threshold > minThreshold){
        #lower threshold
        threshold <- threshold - 0.05
      }else{
        #break if at min threshold and no longer changing
        break
      }
    }
    lastMods <- complexMods
  }

  ###################################
  ###Reassigned flawed assignments###
  ###################################

  #have yet to define a methodology that outperforms skipping this step

  return(complexMods)
}

##########
###TEST###
##########

# source("/restricted/projectnb/agedisease/personal/lberger/modular_graph_learning/ModularDAC/00.SimulateGraphs.R")
#
# er <- make_modular_graph()
# x <- sim_graph_data(er, n.samples = 100)
# mods <- complexICA(x, 3)



##################################################
###Alternative methods for defining core module###
##################################################

###Define core modules bases on some threshold of relative component score
# complexMods[
#   which(!apply(sMat, 1, function(row) {
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
# lastMods <- complexMods
# threshold = 0.95 #will make proper args in function later
# minThreshold = 0
# unassigned <- which(complexMods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each unassigned node
#   for(node in unassigned){
#     #get avg covariance with all current assigned nodes
#     covars <- abs(covMat[node, ])
#     modCovars <- split(covars,lastMods)
#     modCovars["0"] <- NULL #remove unassigned
#     avgCovars <-lapply(modCovars, mean)
#
#     #assign based on max covariance and ICA score
#     maxCov <- max(unlist(avgCovars))
#     wchMax <- which.max(avgCovars)
#     if(maxCov > threshold ){
#       if(useICA){
#         if(which.max(sMat[node,]) == wchMax){
#           complexMods[node]  <- wchMax #might have issues with name matching
#         }
#       }else{
#         complexMods[node]  <- wchMax
#       }
#
#     }
#
#   }
#   #update based on new assignments
#   unassigned <- which(complexMods == 0)
#
#   #if no changes
#   if(all(lastMods == complexMods)){
#     #and not at min threshold
#     if(threshold > minThreshold){
#       #lower threshold
#       print("lowering threshold")
#       threshold <- threshold - 0.05
#     }else{
#       #break if at min threshold and no longer changing
#       break
#     }
#   }
#   lastMods <- complexMods
#   print(threshold)
# }

### Assignments based on closest neighbor
#while there are unassigned nodes
# lastMods <- complexMods
# unassigned <- which(complexMods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each assigned node
#   assigned <- which(complexMods != 0)
#   unassigned <- which(complexMods == 0)
#   for(node in assigned){
#     #assign its "best" unassigned neighbor (based on covar) the same module
#     best <- unassigned[which.max(abs(covMat[node, unassigned]))]
#     complexMods[best] <- complexMods[node]
#   }
#   #end iterations if stuck
#   if(all(lastMods == complexMods)){break}
#   lastMods <- complexMods
# }

### Assignments based on closeness to entire module
#while there are unassigned nodes
# lastMods <- complexMods
# unassigned <- which(complexMods == 0)
# while(length(unassigned) > 0){ #there might be a more efficient way to check this
#   #for each module
#   unassigned <- which(complexMods == 0)
#   for(module in 1:nComp){ #what module goes last can steal nodes, need to store values then take max
#     subCov <- covMat[which(complexMods == module), unassigned]
#     #find the unassigned node closest to the overall module
#     complexMods[unassigned[which.max(colMeans(abs(subCov)))]] <- module
#
#   }
#   #end iterations if stuck
#   if(all(lastMods == complexMods)){break}
#   lastMods <- complexMods
# }

###############################################
###Methods for correcting flawed assignments###
###############################################

### Reassign border nodes based on second degree neighbors
# threshold <- 0.2
# lastMods <- complexMods
#
# #border nodes are nodes with neighbors of a different module
# borderNodes <- which(unlist(
#   lapply(complexMods, function(node){
#     all(unique(complexMods[which(abs(covMat[node, ]) > threshold)]) == complexMods[node])
#   })
# ))
# while(length(borderNodes > 0)){
#   for(node in borderNodes){
#     #get neighbors (based on threshold) and neighbors neighbors
#     n1 <- which(abs(covMat[node, ]) > threshold)
#     n2 <- unique(unlist(lapply(n1, function(n){
#         which(abs(covMat[n, ]) > threshold)
#       })))[-node]
#     #assign node to mode of n2 modules
#     md <- DescTools::Mode(lastMods[n2])[1]
#     if(length(md) == 1){
#       if(!is.na(md)){
#         complexMods[node] <- md
#       }
#     }
#   }
#   if(all(lastMods == complexMods)){break}
#   lastMods <- complexMods
#   borderNodes <- which(unlist(
#     lapply(complexMods, function(node){
#       all(unique(complexMods[which(abs(covMat[node, ]) > threshold)]) == complexMods[node])
#     })
#   ))
# }

###Reassign border nodes based on covariance with final modules
# #border nodes are those whose closest neighbor is of a different module
# borderNodes <- which(complexMods != complexMods[apply(abs(covMat),1, which.max)])
# lastMods <- complexMods
# while(length(borderNodes > 0)){
#   for(node in borderNodes){
#     #get avg covariance with all modules
#     covars <- abs(covMat[node, ])
#     modCovars <- split(covars,lastMods)
#     modCovars["0"] <- NULL #remove unassigned (if any still exist)
#     avgCovars <-lapply(modCovars, mean)
#
#     #assign based on max covariance
#     complexMods[node] <- which.max(avgCovars)
#   }
#   if(all(lastMods == complexMods)){break}
#   lastMods <- complexMods
#   borderNodes <- which(complexMods != complexMods[apply(abs(covMat),1, which.max)])
# }

