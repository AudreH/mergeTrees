hcToPath = function(hc.object){
  successives.steps = list()
  merge1 = hc.object$merge
  order1 = -hc.object$order # topolical order for a tree
  height1 = hc.object$height
  vect = list()

  # vectors = hcToPath_cpp(successives.steps, merge1, order1, n = length(order1))

  for(i in 1:nrow(merge1)){
    # Cas ou les deux merge sont negatifs : deux singletons.
    if(!any(merge1[i,]>0)){
      # successives.steps[[i]] = c(-match(merge1[i,1]), )
      el1 = match(merge1[i,1], -order1)
      el2 = match(merge1[i,2], -order1)
      position = ifelse(el1 > el2, 1, 0)
      if(position == 1){# le max est el1, le min el2
        successives.steps[[i]] = c(el2, el1)
        vect[[i]]= c(el2, el2, el1, height1[i]) # el2 est le split
      }else{
        successives.steps[[i]] = c(el1, el2)
        vect[[i]] = c(el1, el1, el2, height1[i])
      }
    }else if(merge1[i,1] < 0 & merge1[i,2] > 0){
      # merge1[i,2] est un cluster, l'autre est un singleton
      el_singleton = match(merge1[i,1], -order1) # le singleton
      el_cluster = successives.steps[[merge1[i,2]]] # elements du groupe deja forme
      position = ifelse(max(el_cluster) < el_singleton, 1, 0)
      if(position == 1){ # ca veut dire que : le cluster est situe avant le singleton
        successives.steps[[i]] = c(el_cluster, el_singleton)
        vect[[i]] = c(min(el_cluster), max(el_cluster), el_singleton, height1[i])
      }else{ # le cluster est situe apres le singleton
        successives.steps[[i]] = c(el_singleton, el_cluster)
        vect[[i]] = c(el_singleton, el_singleton, max(el_cluster), height1[i])
      }
    }else if(merge1[i,1] > 0 & merge1[i,2] < 0){
      # merge[i,1] est un cluster, l'autre est un singleton.
      # Meme code que precedent.
      el_singleton = match(merge1[i,2], -order1) # le singleton
      el_cluster = successives.steps[[merge1[i,1]]] # elements du groupe deja forme
      position = ifelse(max(el_cluster) < el_singleton, 1, 0)
      if(position == 1){ # ca veut dire que : le cluster est situe avant le singleton
        successives.steps[[i]] = c(el_cluster, el_singleton)
        vect[[i]] = c(min(el_cluster), max(el_cluster), el_singleton, height1[i])
      }else{ # le cluster est situe apres le singleton
        successives.steps[[i]] = c(el_singleton, el_cluster)
        vect[[i]] = c(el_singleton, el_singleton, max(el_cluster), height1[i])
      }
    }else if(merge1[i,1] > 0 & merge1[i,2] > 0){
      # les deux sont des clusters a regrouper
      el_cluster1 = successives.steps[[merge1[i,1]]]
      el_cluster2 = successives.steps[[merge1[i,2]]]
      position = ifelse(max(el_cluster1) < min(el_cluster2), 1, 0)
      if(position == 1){ # ca veut dire que le cluster 1 est place avant le cluster 2
        successives.steps[[i]] = c(el_cluster1, el_cluster2)
        vect[[i]] = c(min(el_cluster1), max(el_cluster1), max(el_cluster2), height1[i])
      }else{
        successives.steps[[i]] = c(el_cluster2, el_cluster1)
        vect[[i]] = c(min(el_cluster2), max(el_cluster2), max(el_cluster1), height1[i])
      }
    }
  }

  path = do.call("rbind", vect1)
  path[, c(1,2,3)] = path[,c(1,2,3)]+1 # changement d'indice entre R et cpp
  return(list(path = path[nrow(path):1,], order = -order1))
}
