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
