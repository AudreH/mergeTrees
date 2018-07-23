


## randomly generate one set of split rules and lambda
## a rule is
## lambda.rules:  from largest to smallest
## order: for data-points
## rules: a matrix of fusion, a line (i, j, k) correspond to fusion of group [i, j] with group ]j, k]
rsplitRules <- function(n){
  lGroup <- list()
  for(i in 1:n){
    lGroup[[i]] <- c(i, i, i)
  }

  mat <- matrix(NA, nrow=n-1, ncol=3)
  for(i in (n-1):1){
    nfuse <- sample.int(n=i, size=1)
    lGroup[[nfuse]][2] <- lGroup[[nfuse]][3]
    lGroup[[nfuse]][3] <- lGroup[[nfuse+1]][3]
    mat[i, ] <- lGroup[[nfuse]]
    lGroup[[nfuse+1]] <- NULL
  }

  return(list(rules=mat, lambda.rules=sort(runif(n-1), decreasing=TRUE), order=sample.int(n)))

}

## randomly generate several set of split rules
## simple iterative call to the previous function
rseveralSplitRules <- function(nb, n){
  return(lapply(rep(n, nb), FUN=rsplitRules))
}



## order splits coming from different set of split rules
## it returns a matrix
## first column  = rule number
## second column = which set of split rules
# orderRules <- function(listRules){
#   n <- length(listRules[[1]]$order)
#   p <- length(listRules)
#   allRulesIndex <- cbind(rep(1:(n-1), p), rep(1:p, each=n-1))
#   return(allRulesIndex[order(sapply(listRules, FUN=function(x) x$lambda.rules), decreasing=TRUE), ])
# }


orderRules <- function(listRules){

  n <- length(listRules[[1]]$order)

  p <- length(listRules)

  allRulesIndex <- cbind(rep(1:(n-1), p), rep(1:p, each = n - 1))

  # vect1 = unlist(lapply(1:length(listRules), FUN = function(x){return(1:length(listRules[[x]]$lambda.rules))}))
  #
  # vect2 = unlist(lapply(1:length(listRules), FUN = function(x){return(rep(x, length(listRules[[x]]$lambda.rules)))}))
  # allRulesIndex <- cbind(vect1, vect2)

  o <- order(sapply(listRules, function(l) l$height), decreasing = TRUE)

  allRulesIndex[o,]

}

