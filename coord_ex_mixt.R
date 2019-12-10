# Coordinate exchange for mixtures
# Based on section 6.3.5 of Goos and Jones
# and Piepel, Cooley and Jones (2005)


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
      cox_direction[n, j] = x[j] - deltas[n]*x[j]/(1 - x[i] + 1e-16)
    }
  }
  
  return(unique(cox_direction))
}



plot_cox_direction = function(x_in, comp = NULL){
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
    out = compute_cox_direction(x_in, comp, 3) %>% 
      as_tibble() %>% 
      set_names(c("x", "y", "z")) %>% 
      ggplot(aes(x, y, z)) +
      coord_tern() + 
      geom_path(linetype = "dashed") + 
      theme_bw() +
      geom_point(data = tibble(x = x_in[1], y = x_in[2], z = x_in[3]))
  } else{
    if(is.null(comp)) comp = 1:length(x_in)
    
    cox_dirs = lapply(1:length(comp), function(i){
      compute_cox_direction(x_in, i, 3) %>% 
        as_tibble() %>% 
        set_names(c("x", "y", "z")) %>% 
        mutate(comp = i)
    }) %>% 
      bind_rows()
    
    # out = cox_dirs %>% 
    #   ggtern(aes(x, y, z)) + 
    #   geom_path(linetype = "dashed", aes(group = comp)) + 
    #   theme_bw() +
    #   geom_point(data = tibble(x = x_in[1], y = x_in[2], z = x_in[3]))
    
    out = cox_dirs %>% 
      ggplot(aes(x, y, z)) +
      coord_tern() + 
      geom_path(linetype = "dashed", aes(group = comp)) + 
      theme_bw() +
      geom_point(data = tibble(x = x_in[1], y = x_in[2], z = x_in[3]))
  }
  
  return(out)
}

library(tidyverse)
library(ggtern)

x_in = c(0.4725440, 0.2627838, 0.2646722)
plot_cox_direction(x_in, 1)
plot_cox_direction(x_in, 1:3)


# plot several Cox's directions for random designs
random_designs = create_random_initial_design(8, q = 3)

cox_direction_plots = lapply(1:nrow(random_designs), function(i){
  return(plot_cox_direction(random_designs[i,], 1:3))
})

# Don't understand why the plots are not right when arranged in a grid with gridExtra
# I need to use ggtern's own grid.arrange
# https://stackoverflow.com/questions/42825983/issue-when-using-grid-arrange-with-ggtern-generated-plots
# do.call("grid.arrange", c(cox_direction_plots, ncol = 4)) # This works with ggtern::grid.arrange
ggtern::grid.arrange(grobs = cox_direction_plots, ncol = 4)



get_scheffe = function(X, order = 1){
  stopifnot(order %in% 1:3)
  if(order == 1)  X_m = X
  else{
    if(order == 2){
      # TODO
      X_m = X
    } else {
      # order = 3
      # TODO
      X_m = X
    }
  }
  return(X_m)
}


get_scheffe_D_efficiency = function(X, order = 1){
  X_m = get_scheffe(X, order = order)
  det_X_m = det(t(X_m) %*% X_m)
  D_eff = det_X_m^(1/ncol(X_m))/nrow(X_m)
  return(D_eff)
}



coord_ex_mixt = function(n_runs, q = 3, n_cox_points = 100, seed = NULL, X = NULL){
  
  if(is.null(X)){
    # If no initial design is provided, create a random starting design
    X = create_random_initial_design(n_runs, q, seed)
    
  } else{
    # Check that rows in X sum to 1
    row_sums = apply(X, 1, sum)
    if(sum(abs(row_sums - rep(1, n_runs)) > 1e-16) != n_runs){
      stop("Rows in X must sum 1")
    }
  }
  
  if(!is.null(seed)) set.seed(seed)
  
  # Design matrix
  # Scheffé model of first order
  # X_m = get_scheffe(X, order = 1)
  # Scheffé model of second order
  # X_m = get_scheffe(X, order = 2)
  # Scheffé model of third order
  # X_m = get_scheffe(X, order = 3)
  
  X_orig = X
  
  # D-efficiency of current design
  d_eff_0 = get_scheffe_D_efficiency(X)
  d_eff_best = d_eff_0
  
  # Coordinate exchanges:
  
  for(k in 1:nrow(X)){
    # cat("k = ", k, "\n")
    for(i in 1:q){
      # cat("i = ", i, "\n")
      
      # The algorithm then evaluates the optimality criterion for a certain number of different designs, say k, obtained by replacing the original coordinate with k equidistant points on the Cox-effect direction line between the lower and upper limit.
      
      # get points from Cox's direction
      cox_direction = compute_cox_direction(X[k,], i, n_cox_points)
      
      # compute optimality criterion for points in Cox's direction
      for(j in 1:nrow(cox_direction)){
        # cat("j = ", j, "\n")
        
        X_new = X
        # replace point of Cox direction in original matrix
        X_new[k,] = cox_direction[j,]
        
        # Get D-efficiency of new design
        d_eff_j = get_scheffe_D_efficiency(X_new, order = 1)
        
        # If new D-efficiency is better, then keep the new one
        if(d_eff_j > d_eff_best) {
          d_eff_best = d_eff_j
          X = X_new
        }
      } # end for Cox
      
    } # end for q
    
  } # end for nrow
  
  return(list(
    X_orig = X_orig,
    X = X,
    d_eff_orig = d_eff_0,
    d_eff = d_eff_best
  ))
  
}




res_alg = coord_ex_mixt(10, q = 3, n_cox_points)

ggtern::grid.arrange(
  res_alg$X_orig %>% 
    as_tibble() %>% 
    set_names(c("x", "y", "z")) %>% 
    ggplot(aes(x, y, z)) +
    geom_point() +
    coord_tern() + 
    theme_minimal() ,
  res_alg$X %>% 
    as_tibble() %>% 
    set_names(c("x", "y", "z")) %>% 
    ggplot(aes(x, y, z)) +
    coord_tern() + 
    geom_point() +
    theme_minimal() ,
  ncol = 2
)

