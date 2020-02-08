
mnl_mixture_coord_ex = function(X, beta, n_cox_points = 100, max_it = 50, plot_designs = F, verbose = 1){
  
  # Add some input checks
  
  q = dim(X)[1]
  
  X_result = mixtureCoordinateExchangeMNL(
    X_orig = X, 
    beta = beta, 
    n_cox_points = n_cox_points, 
    max_it = max_it, 
    verbose = verbose
  )
  
  if(verbose >= 1){
    cat("\n")
    cat("Original log D-efficiency: ", X_result$d_eff_orig)
    cat("\n")
    cat("Final log D-efficiency: ", X_result$d_eff)
    cat("\n")
    cat("Number of iterations: ", X_result$n_iter)
    cat("\n")
  } 
  
  
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
  library(ggtern)
  
  # res_alg: output of a call to mnl_mixture_coord_ex() function
  
  # q = dim(res_alg$X_orig)[1]
  S = dim(res_alg$X_orig)[3]
  
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