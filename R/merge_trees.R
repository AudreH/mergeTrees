hcToPath = function(hc.object){
  successives.steps = list()
  merge1 = hc.object$merge
  order1 = -hc.object$order # topolical order for a tree
  height1 = hc.object$height
  vect = list()

  vectors = hcToPath_cpp(successives.steps, merge1, order1, n = length(order1))

  path = do.call("cbind", vectors)
  path = path +1 # Cpp starts from 0, R starts from 1
  path = cbind(path, height1)

  return(list(path = path[nrow(path):1,], order = -order1))
}

#' Merge a set of hclust objet into a single tree
#'
#' @param hc.list a list with a least one hclust object to be merge in a single consensus tree
#' @param standardize a boolean indicating wether the heights of the different trees should be normalized before merged
#' @export
mergeTrees = function(hc.list, standardize = FALSE){
  # library(Rcpp)
  # sourceCpp("src/prune_splits.cpp")
  # sourceCpp("src/hcToPath.cpp")
  # source("R/agregation_to_hc.R")
  # source("R/generate_splits.R")
  #
  # sourceCpp("src/createMergeMatrix.cpp")
  #
  # hc_1 <- hclust(dist(iris[, 1:4], "euclidean"), method = "ward.D2")
  # hc_2 <- hclust(dist(iris[, 1:4], "euclidean"), method = "complete")
  # tree.list = list(hc_1, hc_2)
  # hc.list = tree.list # TODO A SUPPRIMER

  n = length(hc.list[[1]]$order) # tous les arbres doivent avoir le meme nombre d'element.
  p = length(hc.list)

  if(standardize){
    hc.list = lapply(hc.list, FUN = function(x){
      x$height = x$height/max(x$height) # tous compris entre 0 et 1
      return(x)
    })
  }

  #############################################
  # ----- Labels : ---------------------------
  #############################################

  # In case the trees have different labels: no merging possible (we do'nt know the corresponding labels between the trees)
  # TO DO : add a break here in case all the labels are not identical to those from the first tree.

  labels_list = lapply(hc.list, FUN = function(x){return(x$labels)})
  # orders_list = lapply(list.trees)

  # reference: les labels/order of the first tree in the list
  if(!is.null(labels_list[[1]])){
    labels.equal = lapply(labels_list, FUN = function(x) identical(x, labels_list[[1]]))

    hc.list = lapply(1:length(hc.list), FUN = function(x){
      if(labels.equal[[x]] == FALSE){
        hc.list[[x]]$labels = hc.list[[x]]$labels[order(match(hc.list[[x]]$labels, hc.list[[1]]$labels))]
        hc.list[[x]]$order = hc.list[[x]]$order[order(match(hc.list[[x]]$labels, hc.list[[1]]$order))]
      }else{
      }
      return(hc.list[[x]])
    })
  }

  #############################################
  # ----- Reconstitution paths : -------------
  #############################################

  DataSets = lapply(hc.list, FUN = function(hc.object) return(hcToPath(hc.object)))

  lSetRules  <- lapply(DataSets, function(path.hc) list(rules = path.hc$path,
                                                        lambda.rules = path.hc$path[,4],
                                                        order = path.hc$order))
  oRules <- orderRules(lSetRules)
  out_agregation <- pruneSplits(listSetRules = lSetRules, orderRules = as.matrix(oRules), n, p)

  index_rules = out_agregation$groupsIRule[-1]+1

  dimrule <- oRules[index_rules, 2]
  lambdaRules = unlist(lapply(1:length(dimrule), FUN = function(x){lSetRules[[oRules[index_rules[x],2]]]$lambda.rules[oRules[index_rules[x],1]]}))


  #############################################
  # ----- Reconstitution hclust : ------------
  #############################################


  mat_element = rbind(1:n, out_agregation$currentGroup)

  # Matrix Children and parents:
  mat_groupes = rbind(1:n, out_agregation$groupsParent, out_agregation$groupsChildCurrent, out_agregation$groupsIRule)

  # Reorder : correspondence between matrix of groups and matrix of elements.
  mat_element2 = mat_element[,order(mat_element[2,])]

  mat_aide = rbind(mat_element2, mat_groupes)
  mat_aide2 = mat_aide[-3,]

  l_element = 1 ; l_groupeAct = 2; l_parent = 3 ; l_enfant = 4 ;

  # Those who fathered 0 don't have any children
  mat_aide_save = mat_aide2
  mat_aide2[l_enfant, which(mat_aide2[l_enfant,]==0)] <- 3*n # IntegerMatrixNeeded
  mat_aide3 = mat_aide2[,order(mat_aide2[l_groupeAct,], decreasing = TRUE)]
  matrice_aide2 = mat_aide3[-5,]
  # save(matrice_aide2, file = "prune_res.RData")
  # matrice_aide3 = matrice_aide2[-4,]

  # Matrice Merge :
  # MatriceMerge = CreationMatriceMerge(n, matrice_aide2)
  # print("J'arrive jusqu'a la creation de la matrice merge!!!")
  mergeMatrix = createMergeMatrix(n, matrice_aide2)
  mergeMatrix[mergeMatrix<=n] <- -mergeMatrix[mergeMatrix<=n] # les elements commencent a 1 en R pas en Cpp
  mergeMatrix[mergeMatrix>n] <- mergeMatrix[mergeMatrix>n]-n

  # Ordre :
  mat_aide2 = mat_aide_save
  mat_aide2[l_enfant, which(mat_aide2[l_enfant,]==0)] <- NA
  mat_aide3 = mat_aide2[,order(mat_aide2[l_groupeAct,], decreasing = TRUE)]
  matrice_aide2 = mat_aide3[-5,]
  Order = OrdreIndividus(matrice_aide2)

  # Height :
  # dans le cas ou les individus sont degroupes artificiellement, completer les hauteurs de coupures.
  # TO DO : ne devrait pas intervenir ici..
  Height <- rev(lambdaRules)
  if(length(Height)!=(n-1)){Height = c(rep(0, n-1-length(Height)), Height)} # pas de raisons ici, mais bon...

  Cluster <- list(merge = mergeMatrix, height = Height, order = Order, labels = hc.list[[1]]$labels)
  Cluster <- unclass(Cluster)
  class(Cluster) <- "hclust"

  return(Cluster)
}

