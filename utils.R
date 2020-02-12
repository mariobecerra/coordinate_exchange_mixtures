# Coordinate exchange for mixtures
# Based on section 6.3.5 of Goos and Jones
# and Piepel, Cooley and Jones (2005)

library(Rcpp)
sourceCpp("utils.cpp")


compute_cox_direction = function(x, comp, n_points = 11){
  # Call C++ function
  cox_direction = computeCoxDirection(x, comp, n_points, verbose = 0)
  
  return(cox_direction)
}


plot_cox_direction = function(x_in, comp = NULL, n_points = 3){
  # x_in: vector of length 3 that sums up to 1
  # comp: which component should be plotted. Can be a scalar or a vector
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(ggtern)
  
  if(length(x_in) != 3) stop("x_in must be of length 3")
  if(sum(x_in) != 1) stop("x_in must sum up to 1")
  
  # If component is null
  if(!is.null(comp) & length(comp) == 1){
    out = compute_cox_direction(x_in, comp, n_points) %>% 
      as_tibble() %>% 
      set_names(c("c1", "c2", "c3")) %>% 
      ggplot(aes(c1, c2, c3)) +
      coord_tern() + 
      geom_path(linetype = "dashed") + 
      theme_minimal() +
      geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3]))
  } else{
    if(is.null(comp)) comp = 1:length(x_in)
    
    cox_dirs = lapply(1:length(comp), function(i){
      compute_cox_direction(x_in, i, 3) %>% 
        as_tibble() %>% 
        set_names(c("c1", "c2", "c3")) %>% 
        mutate(comp = i)
    }) %>% 
      bind_rows()
 
    out = cox_dirs %>% 
      ggplot(aes(c1, c2, c3)) +
      coord_tern() + 
      geom_path(linetype = "dashed", aes(group = comp)) + 
      theme_minimal() +
      geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3]))
  }
  
  return(out)
}



create_random_initial_design_gaussian = function(n_runs, q, seed = NULL){
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)
  
  if(!is.null(seed)) set.seed(seed)
  
  for(i in 1:nrow(X)){
    rands = runif(q)
    
    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }
  
  return(X)
}

mixture_coord_ex_gaussian = function(X, order = 1, n_cox_points = 100, max_it = 50, plot_designs = F, verbose = 1){
  
  n_runs = nrow(X)
  q = ncol(X)
  
  # Check that rows in X sum to 1
  row_sums = apply(X, 1, sum)
  
  if(sum(abs(row_sums - rep(1, n_runs)) < 1e-10) != n_runs){
    stop("Rows in X must sum 1")
  }
  
  
  # Coordinate exchanges:
  X_result = mixtureCoordinateExchangeGaussian(X, order, n_cox_points, max_it, verbose)
  
  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    d_eff_orig = X_result$d_eff_orig,
    d_eff = X_result$d_eff,
    n_iter = X_result$n_iter
  )
  
  if(plot_designs) {
    if(q == 3) plot_result_gaussian(out_list)
    else warning("Could not plot results because q != 3")
  }
  return(out_list)
  
}



plot_result_gaussian = function(res_alg){
  # res_alg: output of a call to mixture_coord_ex_gaussian() function
  
  ggtern::grid.arrange(
    res_alg$X_orig %>% 
      as_tibble() %>% 
      set_names(c("c1", "c2", "c3")) %>% 
      ggplot(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      coord_tern() + 
      theme_minimal() +
      ggtitle(label = "",
              subtitle = paste0("log D-efficiency = ", round(res_alg$d_eff_orig, 3)))
    ,
    res_alg$X %>% 
      as_tibble() %>% 
      set_names(c("c1", "c2", "c3")) %>% 
      ggplot(aes(c1, c2, c3)) +
      coord_tern() + 
      geom_point(shape = "x", size = 4) +
      theme_minimal() +
      ggtitle(label = "",
              subtitle = paste0("log D-efficiency = ", round(res_alg$d_eff, 3)))
    ,
    ncol = 2
  )
  
  
}





# optimize_cox_direction = function(
#   X_in, k, i, log_d_eff_best,
#   n_cox_points = 100,
#   order = 1,
#   method = NULL,
#   fn = NULL, gr = NULL, 
#   verbose = 0){
#   
#   X = X_in
#   
#   # Convergence parameter flag
#   optim_conv = 0
#   
#   if(is.null(fn)){
#     
#     # If fn is not given, then use a discrete approximation to Cox's direction and find optimal values
#     
#     # get points from Cox's direction
#     cox_direction = compute_cox_direction(X[k,], i, n_cox_points)
#     # X = findBestCoxDirOld(cox_direction, X, k, order, log_d_eff_best)
#     X = findBestCoxDir(cox_direction, X, k, order, log_d_eff_best)
#     log_d_eff_best = getScheffeLogDEfficiency(X, order);
#     
#     
#   } # end if
#   else{
#     if(is.null(method)){
#       warning("Optimizing method missing. L-BFGS-B will be used.")
#       method = "L-BFGS-B"
#     } 
#     
#     trace = verbose > 1
#     
#     optim_res = optim(
#       X[k,i],
#       fn = fn_cox_scheffe,
#       # gr = gr_cox_scheffe,
#       X = X, k = k, i = i, order = order,
#       method = method,
#       lower = 0.0,
#       upper = 1.0,
#       control = list(maxit = 10000, trace = trace))
#     
#     # optim_res = optim(
#     #   0.5,
#     #   fn = fn_cox_scheffe,
#     #   gr = gr_cox_scheffe,
#     #   X = X, k = k, i = i, order = order,
#     #   method = "L-BFGS")
#     # 
#     # optim_res = optim(
#     #   0.5, 
#     #   fn = fn_cox_scheffe, 
#     #   # gr = gr_cox_scheffe,
#     #   gr = NULL,
#     #   X = X, k = k, i = i, order = order,
#     #   method = "Brent", 
#     #   lower = 0.0,
#     #   upper = 1.0)
#     
#     optim_conv = optim_res$convergence
#     
#     if(optim_conv != 0) {
#       warning('Function "optim" failed in iteration: ',
#               'k = ', k,
#               ', i = ', i,
#               '.\nMessage: "', optim_res$message,
#               '"\nUsing finite approximation for this iteration instead.')
#     }
#     x = optim_res$par
#     x_row = X[k,]
#     
#     delta = x - x_row[i]
#     # recompute proportions in cox dir
#     for(j in setdiff(1:q, i)){
#       # In case it's a corner case, i.e., x[i] = 1
#       # The corner case might be wrong. Check later.
#       if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
#       else{
#         res = x_row[j] - delta*x_row[j]/(1 - x_row[i])
#       }
#       x_row[j] = res
#     } # end j
#     x_row[i] = x
#     
#     X[k,] = x_row
#     log_d_eff_best = -optim_res$value
#   }
#   
#   return(list(
#     X = X, 
#     log_d_eff_best = log_d_eff_best,
#     optim_conv = optim_conv))
#   
# }
# 
# 
# 
# 
# fn_cox_scheffe = function(x, X, k, i, order){
#   
#   x_row = X[k,]
#   
#   delta = x - x_row[i]
#   # recompute proportions in cox dir
#   for(j in setdiff(1:q, i)){
#     # In case it's a corner case, i.e., x[i] = 1
#     # This may be wrong. Check later.
#     if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
#     else{
#       res = x_row[j] - delta*x_row[j]/(1 - x_row[i])
#     }
#     x_row[j] = res
#   } # end j
#   x_row[i] = x
#   
#   Y = X
#   Y[k,] = x_row
#   
#   utility_funct = get_scheffe_log_D_efficiency(Y, order = order)
#   return(-as.numeric(utility_funct))
# }



























