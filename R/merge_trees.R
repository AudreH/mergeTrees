library(Rcpp)
sourceCpp("prune_splits.cpp")
sourceCpp("hcToPath.cpp")
source("generate_splits.R")
source("agregation_to_hc.R")

hcToPath = function(hc.object){
  successives.steps = list()
  merge1 = hc.object$merge
  order1 = -hc.object$order # topolical order for a tree
  height1 = hc.object$height
  vect = list()
  
  vectors = hcToPath_cpp(successives.steps, merge1, order1, n = length(order1))
  
  # for(i in 1:nrow(merge1)){
  #   # Cas ou les deux merge sont negatifs : deux singletons.
  #   if(!any(merge1[i,]>0)){
  #     # successives.steps[[i]] = c(-match(merge1[i,1]), )
  #     el1 = match(merge1[i,1], -order1)
  #     el2 = match(merge1[i,2], -order1)
  #     position = ifelse(el1 > el2, 1, 0)
  #     if(position == 1){# le max est el1, le min el2
  #       successives.steps[[i]] = c(el2, el1)
  #       vect[[i]]= c(el2, el2, el1, height1[i]) # el2 est le split
  #     }else{
  #       successives.steps[[i]] = c(el1, el2)
  #       vect[[i]] = c(el1, el1, el2, height1[i])
  #     }
  #   }else if(merge1[i,1] < 0 & merge1[i,2] > 0){
  #     # merge1[i,2] est un cluster, l'autre est un singleton
  #     el_singleton = match(merge1[i,1], -order1) # le singleton
  #     el_cluster = successives.steps[[merge1[i,2]]] # elements du groupe deja forme
  #     position = ifelse(max(el_cluster) < el_singleton, 1, 0) 
  #     if(position == 1){ # ca veut dire que : le cluster est situe avant le singleton 
  #       successives.steps[[i]] = c(el_cluster, el_singleton)
  #       vect[[i]] = c(min(el_cluster), max(el_cluster), el_singleton, height1[i]) 
  #     }else{ # le cluster est situe apres le singleton
  #       successives.steps[[i]] = c(el_singleton, el_cluster)
  #       vect[[i]] = c(el_singleton, el_singleton, max(el_cluster), height1[i])
  #     }
  #   }else if(merge1[i,1] > 0 & merge1[i,2] < 0){
  #     # merge[i,1] est un cluster, l'autre est un singleton. 
  #     # Meme code que precedent.
  #     el_singleton = match(merge1[i,2], -order1) # le singleton
  #     el_cluster = successives.steps[[merge1[i,1]]] # elements du groupe deja forme
  #     position = ifelse(max(el_cluster) < el_singleton, 1, 0) 
  #     if(position == 1){ # ca veut dire que : le cluster est situe avant le singleton 
  #       successives.steps[[i]] = c(el_cluster, el_singleton)
  #       vect[[i]] = c(min(el_cluster), max(el_cluster), el_singleton, height1[i]) 
  #     }else{ # le cluster est situe apres le singleton
  #       successives.steps[[i]] = c(el_singleton, el_cluster)
  #       vect[[i]] = c(el_singleton, el_singleton, max(el_cluster), height1[i])
  #     }
  #   }else if(merge1[i,1] > 0 & merge1[i,2] > 0){
  #     # les deux sont des clusters a regrouper
  #     el_cluster1 = successives.steps[[merge1[i,1]]] 
  #     el_cluster2 = successives.steps[[merge1[i,2]]]
  #     position = ifelse(max(el_cluster1) < min(el_cluster2), 1, 0)
  #     if(position == 1){ # ca veut dire que le cluster 1 est place avant le cluster 2
  #       successives.steps[[i]] = c(el_cluster1, el_cluster2)
  #       vect[[i]] = c(min(el_cluster1), max(el_cluster1), max(el_cluster2), height1[i])
  #     }else{
  #       successives.steps[[i]] = c(el_cluster2, el_cluster1)
  #       vect[[i]] = c(min(el_cluster2), max(el_cluster2), max(el_cluster1), height1[i])
  #     }
  #   }
  # }
  
  path = do.call("cbind", vectors)
  path = path +1 # Cpp starts from 0, R starts from 1
  path = cbind(path, height1)
  # path = do.call("rbind", vect1)
  # path[, c(1,2,3)] = path[,c(1,2,3)]+1 # changement d'indice entre R et cpp
  return(list(path = path[nrow(path):1,], order = -order1))
}

merge.trees = function(hc.list, standardize = FALSE){
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
  
  labels_list = lapply(list.trees, FUN = function(x){return(x$labels)})
  # orders_list = lapply(list.trees)
  
  # reference: les labels/order of the first tree in the list
  if(!is.null(labels_list[[1]])){
    labels.equal = lapply(labels_list, FUN = function(x) identical(x, labels_list[[1]]))
    
    list.trees = lapply(1:length(list.trees), FUN = function(x){
      if(labels.equal[[x]] == FALSE){
        list.trees[[x]]$labels = list.trees[[x]]$labels[order(match(list.trees[[x]]$labels, list.trees[[1]]$labels))]
        list.trees[[x]]$order = list.trees[[x]]$order[order(match(list.trees[[x]]$labels, list.trees[[1]]$order))]
      }else{
      }
      return(list.trees[[x]])
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
  mat_aide2[l_enfant, which(mat_aide2[l_enfant,]==0)] <- NA
  mat_aide3 = mat_aide2[,order(mat_aide2[l_groupeAct,], decreasing = TRUE)]
  matrice_aide2 = mat_aide3[-5,]
  
  # Matrice Merge : 
  MatriceMerge = CreationMatriceMerge(n, matrice_aide2)
  
  # Ordre :
  Order = OrdreIndividus(matrice_aide2)
  
  # Height :
  # dans le cas ou les individus sont degroupes artificiellement, completer les hauteurs de coupures.
  # TO DO : ne devrait pas intervenir ici.. 
  Height <- rev(lambdaRules)
  if(length(Height)!=(n-1)){Height = c(rep(0, n-1-length(Height)), Height)} # pas de raisons ici, mais bon... 
  
  Cluster <- list(merge = MatriceMerge, height = Height, order = Order, labels = list.trees[[1]]$labels)
  Cluster <- unclass(Cluster)
  class(Cluster) <- "hclust"
  
  return(Cluster)
}

