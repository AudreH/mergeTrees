#' Create a consensus tree.
#'
#' The function assume all of the data.frame or matrices have similar dimensions: for directClustering methods, the datasets are merged using "cbind", it is expected that the matrices have the same number of rows.
#' The class of the objects in the list and the symmetry of the matrices, are determined by the first object in the list.
#'
#' @section Warning: There is no check step for the rownames/labels for the data.frame and trees. Please make sure the distance matrices or data.frame in the list are so the rows are ordered in the same ways.
#'
#' @param lis a list with at least one object to be merge in a single consensus tree. Can be a list of data.frame, distances matrices, matrices or hclust objects.
#' If the matrices are symmetrical, they will be considered as distance matrices unless the method is "directClustering". The function is expecting that the dataframe have the same number of rows,
#' that the distance matrices have the same dimensions and that the trees have the same number of leaves. It is expected that the
#' @param method the method to use for building the consensus tree. Can be "directClustering", "distanceAverage" or "mergeTrees". If NULL, then the method will be chosen according to the type of objects in lis.
#' @param tree_st a boolean indicating wether the trees need to be standardized before merging (only apply to mergeTrees method)
#' @param dist_st a boolean indicating wether the distance matrices need to be standardized before averaging them or performing clustering on them. Apply to methods mergeTrees and "distanceAverage."
#' @param data_st a boolean indicating wether the data need to be standardized (column wise). Used if lis is a list of aSymmetrical dataSets or if "directClustering" is the method.
#' @param data_center a boolean. Used if lis is a list of aSymmetrical dataSets or if "directClustering" is the method.
#' @param distance_method distance used to create distance matrices. Default is "euclidean". Not used if objects are trees or distance matrices.
#' @param linkage_method linkage criteria used to create the trees. Default is "ward.D2". Not used if objects are trees.
#' @param verbose boolean. Display (or not) what the method is doing. Default is FALSE.
#' @return A list of class hclust, being the consensus tree, with the following components: height, merge, method, order, and, if any, labels. For more information about these components, please see hclust help page.
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{dist}}, \code{\link{mergeTrees}}, \code{\link[ape]{consensus}}
#' @author Audrey Hulot, \email{audrey.hulot@@inra.fr}, Julien Chiquet, Guillem Rigaill
#' @export


consensusTree = function(lis,
                         method = NULL,
                         data_st = FALSE,
                         dist_st = FALSE,
                         tree_st = FALSE,
                         data_center = TRUE,
                         distance_method = "euclidean",
                         linkage_method = "ward.D2",
                         verbose = FALSE){

  # method = NULL

  # browser()
  ##########################################
  # --- Check the objects of the list ---  #
  ##########################################
  classObjects = class(lis[[1]])

  if(any(lapply(lis, class)!=classObjects)) stop("Objects in the list are not all of the same class.")

  if(!classObjects%in%c("data.frame", "hclust", "matrix", "dist")) stop("Only support data.frame, hclust, dist or matrix objects.")

  if(classObjects == "data.frame") lis = lapply(lis, as.matrix) # transform into matrix.

  ##########################################
  # --- Chosing the method to use -------  #
  ##########################################

  if(!is.null(method)){
    if(length(grep("direct", method, ignore.case = TRUE))>0) method = "directClustering"
    if(length(grep("average", method, ignore.case = TRUE))>0 | length(grep("distance", method, ignore.case = TRUE))>0) method = "distanceAverage"
    if(length(grep("tree", method, ignore.case = TRUE))>0 | length(grep("merge", method, ignore.case = TRUE))>0) method = "mergeTrees"

    # if(verbose) cat(paste0("Method: ", method, "\n"))
  }

  # if method is not provided by the user, then the method is determined by the class of the objects
  if(is.null(method)){
    if(classObjects == "data.frame" | classObjects == "matrix"){
      if(!isSymmetric(lis[[1]])) method = "directClustering"
      else method = "distanceAverage"
    }
    if(classObjects == "dist")  method = "distanceAverage"
    if(classObjects == "hclust") method = "mergeTrees"

    if(verbose) cat(paste0("Chosing: ", method, "\n"))

  }else{ # check compatibility between class of objects and method
    if((classObjects == "hclust" | classObjects == "dist") & method == "directClustering" ) stop("Can't perform directClustering on distance matrices or hclust objects")
    if(classObjects == "hclust" & method == "distanceAverage") stop("Can't perform averageClustering on hclust objects")
  }


  ##########################################
  # --- Direct Clustering ---------------- #
  ##########################################

  if(method == "directClustering" & classObjects %in% c("data.frame", "matrix")){
    # Verbose
    if(verbose) cat("Performing Direct Clustering \n")
    if(verbose) cat(paste0("* Center the data: ", data_center, "\n"))
    if(verbose) cat(paste0("* Scale the data: ", data_st, "\n"))

    # Standardization
    lis = lapply(lis, scale, center = data_center, scale = data_st)

    # Average distance clustering
    hc = hclust(dist(Reduce("cbind", lis), method = distance_method), method = linkage_method)
    hc$method = "directClustering"

  }else if(method == "directClustering" & classObjects%in%c("dist", "hclust")){
    stop("Can't perform direct clustering on distance or hclust object")
  }

  ##########################################
  # --- Average Distance ----------------- #
  ##########################################

  if(method == "distanceAverage" & classObjects %in% c("data.frame", "matrix", "dist")){
    # Verbose (and scaling data if needed)
    if(verbose) cat("Performing Average Distance Clustering \n")
    if(classObjects!="dist"){
      if(!isSymmetric(lis[[1]])){
        if(verbose)  cat(paste0("* Center data: ", data_center, "\n")) ; cat(paste0("* Standardize data: ", data_st, "\n"))
      }
    }
    if(verbose) cat(paste0("* Standardize distance: ", dist_st, "\n"))

    # Check classes and transform data if needed
    # if(classObjects=="data.frame") lis = lapply(lis, matrix) # not needed if line 40
    if(classObjects!="dist"){
      if(!isSymmetric(lis[[1]]) & (data_st | data_center)) lapply(lis, scale, center = data_center, scale = data_st)

      lis = lapply(lis, dist, method = distance_method)
    }

    # Standardization
    if(dist_st) lis = lapply(lis, FUN = function(mat) mat/max(mat))

    # Average distance clustering
    hc = hclust(as.dist(Reduce("+", lis)/length(lis)), method = linkage_method)
    hc$method = "distanceAverage"

  }else if(method == "distanceAverage" & classObjects=="hclust"){
    stop("Can't perform averageDistance clustering on trees")
  }

  ##########################################
  # --- Merge Trees ---------------------- #
  ##########################################

  if(method == "mergeTrees"){ # classObject can be any of data.frame, matrix, dist, hclust.
    if(verbose) cat("Performing merge tree method \n")
    if(!classObjects%in%c("dist", "hclust")){
      if(!isSymmetric(lis[[1]])){
        if(verbose)  cat(paste0("* Center data: ", data_center, "\n")) ; cat(paste0("* Standardize data: ", data_st, "\n"))
      }
    }
    if(verbose & classObjects!="hclust") cat(paste0("* Standardize distance: ", dist_st, "\n"))
    if(verbose) cat(paste0("* Standardize trees: ", tree_st, "\n"))

    if(classObjects!="hclust"){
      if(classObjects!="dist"){
        if(!isSymmetric(lis[[1]]) & (data_st | data_center))lapply(lis, scale, center = data_center, scale = data_st)
        lis = lapply(lis, dist, method = distance_method)
      }
      if(dist_st) lis = lapply(lis, FUN = function(mat) mat/max(mat))
      lis = lapply(lis, FUN = function(dist_mat) hclust(as.dist(dist_mat), method = linkage_method))
    }

    hc = mergeTrees(lis, standardize = tree_st)
  }

  return(hc)

}



