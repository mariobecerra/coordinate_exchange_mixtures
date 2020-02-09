
mnl_mixture_coord_ex = function(X, beta, n_cox_points = 100, max_it = 50, plot_designs = F, verbose = 1){
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
  # res_alg: output of a call to mnl_mixture_coord_ex() function.
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


