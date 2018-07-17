#####################################################-
# Transformation en matrice merge type hclust$merge #
#####################################################-

# NOTE 17/07/18: cette fonction a ete recodee en cpp. Il y a beaucoup de choses qui sont inutiles
#   et beaucoup de lignes ne servent a rien. Si jamais elle est reprise : enlever ce qui ne sert a rien...

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

    groupe_act_old = matrice_aide[l_groupeAct,1] # inutile
    groupe_act_new = matrice_aide[l_groupeAct, which_merge] # inutile

    matrice_aide[l_enfant,which(matrice_aide[l_enfant,]==groupe_act_old)] = groupe_act_new # inutile
    merge_mat[iterations, 1] = elementCourant
    merge_mat[iterations, 2] = merge_element[l_element]
    matrice_aide[l_element, which_merge] = elementSuivant
    matrice_aide[l_groupeAct, which_merge] = merge_element[l_groupeAct] #inutile
    matrice_aide[l_parent, which_merge] = merge_element[l_parent] # inutile
    matrice_aide[l_enfant, which_merge] = merge_element[l_enfant] # inutile
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
