hcToPath = function(hc.object, n = length(hc.object$order)){
  # hc.object = hc_1
  merge1 = hc.object$merge
  # order1 = -hc.object$order # topolical order for a tree # COMMENT ON 19/07/18, UNCOMMENT IF NOT CONCLUSIVE
  order1 = hc.object$order # 19/07/18 : REMOVE IF NOT CONCLUSIVE
  # vectors = hcToPath_cpp(successives.steps, merge1, order1, n = length(order1)) # UNCOMMENT


  # AJOUT 19/07/18 : attempt to speed up code. Remove if not ok
  match_order = match(1:max(order1), order1)
  vectors = hcToPath_cpp(merge1, match_order, n)

  path = do.call("cbind", vectors)
  # path = path +1 # Cpp starts from 0, R starts from 1 // 19/07/18, UNCOMMENT IF NOT CONCLUSIVE
  path = cbind(path, hc.object$height)

  return(list(path = path[nrow(path):1,], order = order1)) # 19/07/18 : add a (-) if not conclusive
}
