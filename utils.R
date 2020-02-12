# Coordinate exchange for mixtures for Gaussian and multinomial logit models. Based on:
# 1) section 6.3.5 of Goos and Jones; 
# 2) Piepel, Cooley and Jones (2005); and
# 3) Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)

library(Rcpp)

sourceCpp("utils.cpp")


##########################################################
##########################################################
## General functions
##########################################################
##########################################################

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



##########################################################
##########################################################
## Functions for Gaussian model
##########################################################
##########################################################

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




##########################################################
##########################################################
## Functions for MNL model
##########################################################
##########################################################

create_random_initial_MNL_design = function(q, J, S, seed = NULL){
  # Creates a random initial design of the specified dimensions
  # Returns a 3-dimensional array of dimensions (q, J, S)
  X = array(rep(NA_real_, q*J*S), dim = c(q, J, S))
  
  if(!is.null(seed)) set.seed(seed)
  
  for(j in 1:J){
    for(s in 1:S){
      rands = runif(q)
      # ingerdients must sum up to 1
      X[,j, s] = rands/sum(rands)
    }
  }
  
  return(X)
}




create_random_beta = function(q){
  # Creates a random parameter vector.
  # Returns a list in which the first element of the list is a numerical vector with the parameters and the second element is a matrix with the indices as in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  
  beta_1 = rnorm(q)
  
  beta_2 = rnorm(q*(q-1)/2)
  
  beta_3 = rnorm(q*(q-1)*(q-2)/6)
  beta = c(beta_1, beta_2, beta_3)
  
  n = (q*q*q+5*q)/6
  
  # Just in case my formula is wrong:
  if(length(beta) != n) stop("Error in dimensions")
  
  beta_ix = matrix(rep(NA_integer_, 3*n), ncol = 3)
  colnames(beta_ix) = c("i", "k", "l")
  
  beta_ix[1:q, "i"] = 1:q
  counter = q
  
  for(i in 1:(q-1)){
    for(k in (i+1):(q)){
      counter = counter + 1
      beta_ix[counter, "i"] = i
      beta_ix[counter, "k"] = k
    }
  }
  
  for(i in 1:(q-2)){
    for(k in (i+1):(q-1)){
      for(l in (k+1):(q)){
        counter = counter + 1
        beta_ix[counter, "i"] = i
        beta_ix[counter, "k"] = k
        beta_ix[counter, "l"] = l
      }
    }
  }
  
  return(list(beta = beta, beta_ix = beta_ix))
}



mixture_coord_ex_mnl = function(X, beta, n_cox_points = 100, max_it = 50, plot_designs = F, verbose = 1){
  # Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model.
  # X: 3 dimensional array with dimensions (q, J, S) where:
  #    q is the number of ingredient proportions
  #    J is the number of alternatives within a choice set
  #    S is the number of choice sets
  # beta: vector of parameters. Should be of length (q^3 + 5*q)/6
  # n_cox_points: Number of points to use in the discretization of Cox direction
  # max_it: Maximum number of iteration that the coordinate exchange algorithm will do
  # plot_designs: If TRUE, shows a plot of the initial and the final design. Only works if q is 3.
  # verbose: level of verbosity. 6 levels, in which level prints the previous plus additional things:
  #    1: Print the log D efficiency in each iteration and a final summary
  #    2: Print the values of k, s, i, and log D efficiency in each subiteration
  #    3: Print the resulting X after each iteration, i.e., after each complete pass on the data
  #    4: Print log D efficiency for each point in the Cox direction discretization
  #    5: Print the resulting X and information matrix after each subiteration
  #    6: Print the resulting X or each point in the Cox direction discretization
  # Returns alist with the following objects:
  #    X_orig: The original design. Array with dimensions (q, J, S).
  #    X: The optimized design. Array with dimensions (q, J, S).
  #    d_eff_orig: log D-efficiency of the original design.
  #    d_eff: log D-efficiency of the optimized design.
  #    n_iter: Number of iterations performed.
  
  
  # Some input checks
  dim_X = dim(X)
  
  if(length(dim_X) != 3) stop("X must be a 3 dimensional array.")
  if(!is.vector(beta)) stop("beta is not a vector. It must be a numerical or integer vector.")
  if(!(is.numeric(beta) | !is.integer(beta))) stop("beta is not numeric or integer. It must be a numerical or integer vector.")
  
  q = dim_X[1]
  m = (q*q*q + 5*q)/6
  
  if(m != length(beta)) stop("Incompatible length in beta and q: beta must be of length (q^3 + 5*q)/6")
  
  # Call to C++ function
  # Note: In the future use a C++ implementation of Brent's method like the following
  # https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.html
  # https://github.com/fditraglia/RcppBrent
  # It has the implementation of the algorithm in Chapter 6 of Brent's book
  # (Ch 6: Global Minimization Given an Upper Bound on the Second Derivative)
  X_result = mixtureCoordinateExchangeMNL(
    X_orig = X, 
    beta = beta, 
    n_cox_points = n_cox_points, 
    max_it = max_it, 
    verbose = verbose
  )
  
  out_list = list(
    X_orig = X_result$X_orig,
    X = X_result$X,
    d_eff_orig = X_result$d_eff_orig,
    d_eff = X_result$d_eff,
    n_iter = X_result$n_iter
  )
  
  if(plot_designs) {
    if(q == 3) mnl_plot_result(out_list)
    else warning("Could not plot results because q != 3")
  }
  return(out_list)
}



mnl_plot_result = function(res_alg){
  # res_alg: output of a call to mixture_coord_ex_mnl() function.
  # It must be a design of 3 ingredients.
  
  dim_X = dim(res_alg$X_orig)
  q = dim_X[1]
  S = dim_X[3]
  
  if(q != 3) stop("Design must be of 3 ingredients.")
  
  library(ggtern)
  
  
  # Convert 3 dimensional arrays into matrices by vertically binding them
  X_orig_mat = t(res_alg$X_orig[,,1])
  for(s in 2:S){
    X_orig_mat = rbind(X_orig_mat, t(res_alg$X_orig[,,s]))
  }
  
  X_final_mat = t(res_alg$X_orig[,,1])
  for(s in 2:S){
    X_final_mat = rbind(X_final_mat, t(res_alg$X[,,s]))
  }
  
  # PLot matrices
  ggtern::grid.arrange(
    X_orig_mat %>% 
      as_tibble() %>% 
      set_names(c("c1", "c2", "c3")) %>% 
      ggplot(aes(c1, c2, c3)) +
      geom_point(shape = "x", size = 4) +
      coord_tern() + 
      theme_minimal() +
      ggtitle(label = "",
              subtitle = paste0("log D-efficiency = ", round(res_alg$d_eff_orig, 3)))
    ,
    X_final_mat %>% 
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



