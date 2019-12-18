# Coordinate exchange for mixtures
# Based on section 6.3.5 of Goos and Jones
# and Piepel, Cooley and Jones (2005)

library(Rcpp)
sourceCpp("utils.cpp")

create_random_initial_design = function(n_runs, q, seed = NULL){
  X = matrix(rep(NA_real_, n_runs*q), nrow = n_runs)
  
  if(!is.null(seed)) set.seed(seed)
  
  for(i in 1:nrow(X)){
    rands = runif(q)
    
    # rows in X sum to 1:
    X[i,] = rands/sum(rands)
  }
  
  return(X)
}


compute_cox_direction = function(x, comp, n_points = 11){
  
  i = comp
  q = length(x)
  
  
  # points in Cox's direction
  seq_points = seq(from = 0, to = 1, length.out = n_points)
  
  diffs = seq_points - x[i]
  ix1 = which.min(abs(diffs))
  
  cox_direction_aux = seq_points - diffs[ix1]
  cox_direction_aux2 = cox_direction_aux[-ix1]
  betw_0_1 = cox_direction_aux2 >= 0 & cox_direction_aux2 <= 1
  cox_direction_aux3 = c(0, cox_direction_aux2[betw_0_1], 1)
  
  cox_direction = matrix(rep(NA_real_, q*length(cox_direction_aux3)), ncol = q)
  cox_direction[,i] = cox_direction_aux3
  deltas = cox_direction_aux3 - x[i]
  
  for(n in seq_along(cox_direction_aux3)){
    # recompute proportions:
    for(j in setdiff(1:q, i)){
      # In case it's a corner case, i.e., x[i] = 1
      if(abs(1 - x[i]) < 1e-16) res = (1 - cox_direction[n, i])/(q-1)
      else{
        res = x[j] - deltas[n]*x[j]/(1 - x[i])
      }
      cox_direction[n, j] = res
    } # end j
    
    if(any(cox_direction[n, ] < -1e-10 | cox_direction[n, ] > 1 + 1e10)) {
      stop("Error while computing Cox direction. ",
           "Value out of bounds.\n", 
           "Cox direction computed:\n\tc(", 
           paste(cox_direction[n, ], collapse = ", "), ")")
    }
  }
  cox_direction = unique(cox_direction)
  
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
    
    # out = cox_dirs %>% 
    #   ggtern(aes(c1, c2, c3)) + 
    #   geom_path(linetype = "dashed", aes(group = comp)) + 
    #   theme_minimal() +
    #   geom_point(data = tibble(x = x_in[1], y = x_in[2], z = x_in[3]))
    
    out = cox_dirs %>% 
      ggplot(aes(c1, c2, c3)) +
      coord_tern() + 
      geom_path(linetype = "dashed", aes(group = comp)) + 
      theme_minimal() +
      geom_point(data = tibble(c1 = x_in[1], c2 = x_in[2], c3 = x_in[3]))
  }
  
  return(out)
}


get_scheffe_order_2 = function(X){
  q = ncol(X)
  n = nrow(X)
  n_col_X_m = q + (q-1)*q/2
  X_m = matrix(rep(NA_real_, n_col_X_m*n), nrow = n)
  X_m[,1:q] = X
  
  k = q
  for(i in 1:(q-1)){
    for(j in (i+1):q){
      # cat("i = ", i, ", j = ", j, "\n")
      k = k+1
      X_m[,k] = X[,i]*X[,j]
    }  
  }
  return(X_m)
}


get_scheffe_order_3 = function(X){
  q = ncol(X)
  n = nrow(X)
  X_m = get_scheffe_order_2(X)
  for(i in 1:(q-2)){
    for(j in (i+1):(q-1)){
      for(k in (j+1):q){
        # cat("i = ", i, ", j = ", j, "\n")
        # not the most efficient, but works
        X_m = cbind(X_m, X[,i]*X[,j]*X[,k])
      }  
    }
  }
  return(X_m)
}


# get_scheffe = function(X, order = 1){
#   stopifnot(order %in% 1:3)
#   if(order == 1)  X_m = X
#   else{
#     if(order == 2){
#       
#       X_m = get_scheffe_order_2(X)
#       
#     } else {
#       # order = 3
#       X_m = get_scheffe_order_3(X)
#     }
#   }
#   return(X_m)
# }


get_scheffe = function(X, order = 1){
  X_m = getScheffe(X, order)
  return(X_m)
}


get_scheffe_log_D_efficiency = function(X, order = 1){
  X_m = get_scheffe(X, order = order)
  # det_X_m = det(t(X_m) %*% X_m)
  log_det_X_m = determinant(t(X_m) %*% X_m)$modulus
  log_D_eff = log_det_X_m/ncol(X_m) - log(nrow(X_m))
  # log_D_eff = det_X_m^(1/ncol(X_m))/nrow(X_m)
  return(as.numeric(log_D_eff))
}




coord_ex_mixt = function(n_runs = 10, q = 3, n_cox_points = 100, order = 1, max_it = 50, seed = NULL, X = NULL, plot_designs = F, verbose = 1, method = "Brent"){
  
  if(is.null(X)){
    # If no initial design is provided, create a random starting design
    X = create_random_initial_design(n_runs, q, seed)
    
  } 
  
  n_runs = nrow(X)
  q = ncol(X)
  
  # Check that rows in X sum to 1
  row_sums = apply(X, 1, sum)
  
  if(sum(abs(row_sums - rep(1, n_runs)) < 1e-10) != n_runs){
    stop("Rows in X must sum 1")
  }
  
  available_methods = c("finite_approximation", "Brent", "L-BFGS-B")
  
  pmatch_vec = sapply(available_methods, 
         function(x) pmatch(method, x))
  
  method = names(which(!is.na(pmatch_vec)))
  
  if(all(is.null(pmatch_vec))){
    # if no method is recognized
    stop('Optimization method unknown. Must be one of the following:\n\t',
         paste(available_methods, collapse = ", "))
  } else{
    if(method == "finite_approximation"){
      fn_optim = NULL
      gr_optim = NULL
      method_optim = NULL
    } else{
      fn_optim = fn_cox_scheffe
      # gr_optim = gr_cox_scheffe
      gr_optim = NULL
      method_optim = method
    }
  }
  
  cat('Optimizing using "', method, '".\n', sep = "")
  
  
  if(!is.null(seed)) set.seed(seed)
  
  # Design matrix
  # Scheffé model of first order
  # X_m = get_scheffe(X, order = 1)
  # Scheffé model of second order
  # X_m = get_scheffe(X, order = 2)
  # Scheffé model of third order
  # X_m = get_scheffe(X, order = 3)
  
  X_orig = X
  
  # log D-efficiency of current design
  log_d_eff_orig = get_scheffe_log_D_efficiency(X, order = order)
  log_d_eff_best = log_d_eff_orig
  log_d_eff_aux = Inf
  
  # Coordinate exchanges:
  it = 0
  while(it < max_it){
    it = it + 1
    if(verbose >= 1) cat("Iter: ", it, ", log D-efficiency: ", log_d_eff_best, "\n")
    # If there was no improvement in this iteration
    if(abs(log_d_eff_aux - log_d_eff_best) < 1e-16) break
    
    log_d_eff_aux = log_d_eff_best
    
    for(k in 1:nrow(X)){
      if(verbose >= 2) cat("k = ", k, "\n")
      for(i in 1:q){
        if(verbose >= 2) cat("i = ", i, "\n")
        
        # The algorithm then evaluates the optimality criterion for a certain number of different designs, say k, obtained by replacing the original coordinate with k equidistant points on the Cox-effect direction line between the lower and upper limit.
        
        optim_list = optimize_cox_direction(
          X_in = X, 
          k = k, 
          i = i, 
          order = order,
          log_d_eff_best = log_d_eff_best,
          fn = fn_optim, gr = gr_optim, method = method,
          n_cox_points = n_cox_points,
          verbose = verbose)
        
        if(optim_list$optim_conv != 0){
          
          cat('Using finite approximation with 1000 points for iteration',
              'k =', k,
              ', i =', i,
              '.\n')
          
          optim_list = optimize_cox_direction(
            X_in = X, 
            k = k, 
            i = i, 
            order = order,
            log_d_eff_best = log_d_eff_best,
            fn = NULL, gr = NULL, method = NULL,
            n_cox_points = 1000,
            verbose = verbose)
        }
        
        
        X = optim_list$X
        log_d_eff_best = optim_list$log_d_eff_best
        
        if(verbose >= 5) {
          print(X)
          cat("Log D-eff:", log_d_eff_best, "\n")
        }
        
      } # end for q
      
    } # end for nrow
    
    if(verbose >= 3) print(X)
    
  } # end while
  
  if(verbose >= 1){
    cat("\n")
    cat("Original log D-efficiency: ", log_d_eff_orig)
    cat("\n")
    cat("Final log D-efficiency: ", log_d_eff_best)
    cat("\n")
    cat("Number of iterations: ", it)
    cat("\n")
  } 
  
  
  out_list = list(
    X_orig = X_orig,
    X = X,
    d_eff_orig = log_d_eff_orig,
    d_eff = log_d_eff_best,
    n_iter = it
  )
  
  if(plot_designs) {
    
    # ggtern::grid.arrange(
    #   out_list$X_orig %>% 
    #     as_tibble() %>% 
    #     set_names(c("c1", "c2", "c3")) %>% 
    #     ggplot(aes(c1, c2, c3)) +
    #     geom_point(shape = "x", size = 2) +
    #     coord_tern() + 
    #     theme_minimal() +
    #     ggtitle(label = "",
    #             subtitle = paste0("log D-efficiency = ", round(out_list$d_eff_orig, 3)))
    #   ,
    #   out_list$X %>% 
    #     as_tibble() %>% 
    #     set_names(c("c1", "c2", "c3")) %>% 
    #     ggplot(aes(c1, c2, c3)) +
    #     coord_tern() + 
    #     geom_point(shape = "x", size = 2) +
    #     theme_minimal() +
    #     ggtitle(label = "",
    #             subtitle = paste0("log D-efficiency = ", round(out_list$d_eff, 3)))
    #   ,
    #   ncol = 2
    # )
    plot_result(out_list)
  }
  
  return(out_list)
  
}



plot_result = function(res_alg){
  # res_alg: output of a call to coord_ex_mixt() function
  
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





optimize_cox_direction = function(
  X_in, k, i, log_d_eff_best,
  n_cox_points = 100,
  order = 1,
  method = NULL,
  fn = NULL, gr = NULL, 
  verbose = 0){
  
  X = X_in
  
  # Convergence parameter flag
  optim_conv = NULL
  
  if(is.null(fn)){
    
    # If fn is not given, then use a discrete approximation to Cox's direction and find optimal values
    
    # get points from Cox's direction
    cox_direction = compute_cox_direction(X[k,], i, n_cox_points)
    
    # compute optimality criterion for points in Cox's direction
    for(j in 1:nrow(cox_direction)){
      if(verbose >= 2) cat("j = ", j, "\n")
      
      X_new = X
      # replace point of Cox direction in original matrix
      X_new[k,] = cox_direction[j,]
      
      # Get D-efficiency of new design
      log_d_eff_j = get_scheffe_log_D_efficiency(X_new, order = order)
      if(verbose >= 2) cat("\t", log_d_eff_j, "\n")
      
      # If new D-efficiency is better, then keep the new one
      if(log_d_eff_j > log_d_eff_best) {
        log_d_eff_best = log_d_eff_j
        X = X_new
      }
    } # end for Cox
    
    
  } # end if
  else{
    if(is.null(method)){
      warning("Optimizing method missing. L-BFGS-B will be used.")
      method = "L-BFGS-B"
    } 
    
    trace = verbose > 1
    
    optim_res = optim(
      X[k,i],
      fn = fn_cox_scheffe,
      # gr = gr_cox_scheffe,
      X = X, k = k, i = i, order = order,
      method = method,
      lower = 0.0,
      upper = 1.0,
      control = list(maxit = 10000, trace = trace))
    
    # optim_res = optim(
    #   0.5,
    #   fn = fn_cox_scheffe,
    #   gr = gr_cox_scheffe,
    #   X = X, k = k, i = i, order = order,
    #   method = "L-BFGS")
    # 
    # optim_res = optim(
    #   0.5, 
    #   fn = fn_cox_scheffe, 
    #   # gr = gr_cox_scheffe,
    #   gr = NULL,
    #   X = X, k = k, i = i, order = order,
    #   method = "Brent", 
    #   lower = 0.0,
    #   upper = 1.0)
    
    optim_conv = optim_res$convergence
    
    if(optim_conv != 0) {
      warning('Function "optim" failed in iteration: ',
              'k = ', k,
              ', i = ', i,
              '.\nMessage: "', optim_res$message,
              '"\nUsing finite approximation for this iteration instead.')
    }
    x = optim_res$par
    x_row = X[k,]
    
    delta = x - x_row[i]
    # recompute proportions in cox dir
    for(j in setdiff(1:q, i)){
      # In case it's a corner case, i.e., x[i] = 1
      # The corner case might be wrong. Check later.
      if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
      else{
        res = x_row[j] - delta*x_row[j]/(1 - x_row[i])
      }
      x_row[j] = res
    } # end j
    x_row[i] = x
    
    X[k,] = x_row
    log_d_eff_best = -optim_res$value
  }
  
  return(list(
    X = X, 
    log_d_eff_best = log_d_eff_best,
    optim_conv = optim_conv))
  
  
}







fn_cox_scheffe = function(x, X, k, i, order){
  
  x_row = X[k,]
  
  delta = x - x_row[i]
  # recompute proportions in cox dir
  for(j in setdiff(1:q, i)){
    # In case it's a corner case, i.e., x[i] = 1
    # This may be wrong. Check later.
    if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
    else{
      res = x_row[j] - delta*x_row[j]/(1 - x_row[i])
    }
    x_row[j] = res
  } # end j
  x_row[i] = x
  
  Y = X
  Y[k,] = x_row
  
  utility_funct = get_scheffe_log_D_efficiency(Y, order = order)
  return(-as.numeric(utility_funct))
}








find_ith_row_inv_A = function(A, i, tol = 1e-32){
  n = nrow(A)
  eye = diag(n)[,i]
  a_i = solve(A, eye, tol = tol)
  return(a_i)
}

find_inv_A_ij = function(A, i, j, tol = 1e-32){
  ith_row = find_ith_row_inv_A(A, 1)
  return(ith_row[j])
}



gr_cox_scheffe = function(x, X, k, i, order){
  x_row = X[k,]
  
  delta = x - x_row[i]
  # recompute proportions in cox dir
  for(j in setdiff(1:q, i)){
    # In case it's a corner case, i.e., x[i] = 1
    # This may be wrong. Check later.
    if(abs(1 - x_row[i]) < 1e-16) res = (1 - x_row[i])/(q-1)
    else{
      res = x_row[j] - delta*x_row[j]/(1 - x_row[i])
    }
    x_row[j] = res
  } # end j
  x_row[i] = x
  
  Y = X
  Y[k,] = x_row
  
  Y = get_scheffe(Y, order = order)
  
  YtY = t(Y) %*% Y
  
  ji_element = find_inv_A_ij(YtY, i, k, tol = 1e-32)
  return(ji_element/q)
}




