#####################################################-
# Transformation en matrice merge type hclust$merge #
#####################################################-

CreationMatriceMerge <- function(n, matrice_aide){
  l_element = 1 ; l_groupeAct = 2; l_parent = 3 ; l_enfant = 4 ;
    
  merge_mat = matrix(NA, ncol = 2, nrow = n-1)
  iterations = 1
  elementSuivant = n+1
  
  while(iterations<(n-1)){
    elementCourant = matrice_aide[l_element,1]
    groupe_parent = matrice_aide[l_parent, 1]
    
    which_merge = which(matrice_aide[l_groupeAct,]==groupe_parent)
    merge_element = matrice_aide[,which_merge]
    
    groupe_act_old = matrice_aide[l_groupeAct,1]
    groupe_act_new = matrice_aide[l_groupeAct, which_merge]
    
    matrice_aide[l_enfant,which(matrice_aide[l_enfant,]==groupe_act_old)] = groupe_act_new
    merge_mat[iterations, 1] = elementCourant
    merge_mat[iterations, 2] = merge_element[l_element]
    matrice_aide[l_element, which_merge] = elementSuivant
    matrice_aide[l_groupeAct, which_merge] = merge_element[l_groupeAct]
    matrice_aide[l_parent, which_merge] = merge_element[l_parent]
    matrice_aide[l_enfant, which_merge] = merge_element[l_enfant]
    matrice_aide = matrice_aide[,-1]
    
    elementSuivant = elementSuivant+1
    iterations = iterations+1
  }
  
  # fusion des deux dernieres colonnes (il est cense n'en rester que deux)
  merge_mat[iterations,] = c(matrice_aide[l_element,1], matrice_aide[l_element,2]) 
  
  merge_mat[merge_mat<=n] <- -merge_mat[merge_mat<=n] # les elements commencent a 1 en R pas en Cpp
  merge_mat[merge_mat>n] <- merge_mat[merge_mat>n]-n
  
  return(merge_mat)
}

#########################################
# Recuperation de l'ordre des individus #
#########################################

OrdreIndividus <- function(mat_aide){
  
  l_element = 1 ; l_groupeAct = 2; l_parent = 3 ; l_enfant = 4 ;
  
  list_groupes_parents = by(mat_aide[l_groupeAct,], INDICES = mat_aide[l_parent,], FUN = function(x) return(x))
  
  # Ordonner les groupes parents et leurs enfants
  first_order = c(0,list_groupes_parents$`0`)# la base, engendre tous les groupes/individus
  # deux fois 0, normal. c'est voulu.
  
  # si tout le monde degroupe de la base : liste de longueur 1, l'ordre est deja etabli.
  
  if(length(list_groupes_parents)>1){
    for(i in 2:length(list_groupes_parents)){
      enfantsDuParent = list_groupes_parents[[i]]
      parent = as.numeric(names(list_groupes_parents)[i])
      indice_replace = which(first_order==parent)
      if(length(indice_replace)>0){
        first_order <- c(first_order[1:(indice_replace)], enfantsDuParent,
                         first_order[(indice_replace+1):length(first_order)]) 
      }else{
        first_order = c(first_order, parent, enfantsDuParent)
      }
    }
  }
  
  Order = unique(first_order[!is.na(first_order)])
  Order = rev(mat_aide[l_element,])[Order+1]
  
  return(Order)
}

# Sortie :
# - Order : l'ordre des individus, tel qu'il n'y aura pas de coupure dans hclust
