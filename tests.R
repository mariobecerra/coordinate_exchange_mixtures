library(tidyverse)
library(ggtern)

source("utils.R")


##########################################################
##########################################################
## Tests functions common to both models
##########################################################
##########################################################

x_in = c(0.4725440, 0.2627838, 0.2646722)
plot_cox_direction(x_in, 1)
plot_cox_direction(x_in, 1:3)

plot_cox_direction(c(1, 0, 0), 1:3)
plot_cox_direction(c(0, 1, 0), 1:3)
plot_cox_direction(c(0, 0, 1), 1:3)




##########################################################
##########################################################
## Gaussian model
##########################################################
##########################################################

# plot several Cox's directions for random designs
# Don't understand why the plots are not right when arranged in a grid with gridExtra
# I need to use ggtern's own grid.arrange
# https://stackoverflow.com/questions/42825983/issue-when-using-grid-arrange-with-ggtern-generated-plots
# do.call("grid.arrange", c(cox_direction_plots, ncol = 4)) # This works with ggtern::grid.arrange
random_designs = create_random_initial_design_gaussian(8, q = 3)

cox_direction_plots = lapply(1:nrow(random_designs), function(i){
  return(plot_cox_direction(random_designs[i,], 1:3))
})

ggtern::grid.arrange(grobs = cox_direction_plots, ncol = 4)


# First degree

X_1 = create_random_initial_design_gaussian(9, 3, seed = 10)

res_alg_order_1_1 = mixture_coord_ex_gaussian(
  X_1, 
  order = 1, 
  plot_designs = T,
  n_cox_points = 200)


# Second degree
res_alg_order_2_1 = mixture_coord_ex_gaussian(
  X_1, 
  order = 2, 
  plot_designs = T,
  n_cox_points = 1000)




# Third degree
res_alg_order_3_1 = mixture_coord_ex_gaussian(
  X_1, 
  order = 3, 
  plot_designs = T,
  n_cox_points = 1000)








X_2 = create_random_initial_design_gaussian(15, 3, seed = 10)

res_alg_order_1_2 = mixture_coord_ex_gaussian(
  X_2, 
  order = 1, 
  plot_designs = T,
  n_cox_points = 200)


# Second degree
res_alg_order_2_2 = mixture_coord_ex_gaussian(
  X_2, 
  order = 2, 
  plot_designs = T,
  n_cox_points = 1000)




# Third degree
res_alg_order_3_2 = mixture_coord_ex_gaussian(
  X_2, 
  order = 3, 
  plot_designs = T,
  n_cox_points = 1000)





##########################################################
##########################################################
## MNL model
##########################################################
##########################################################


q = 3
J = 5
S = 3
X2 = create_random_initial_MNL_design(q, J, S, seed = 4)
beta2 = rep(0, (q*q*q + 5*q)/6)



X2_opt2 = mixture_coord_ex_mnl(
  X = X2, 
  beta = beta2, 
  n_cox_points = 5, 
  max_it = 2, 
  verbose = 5,
  plot_designs = T
)










q = 3
J = 5
S = 4
X3 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta3 = rep(0, (q*q*q + 5*q)/6)
beta3_2 = create_random_beta(q)

X3_opt = mixture_coord_ex_mnl(
  X = X3, 
  beta = beta3, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)


X3_2_opt = mixture_coord_ex_mnl(
  X = X3, 
  beta = beta3_2$beta, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)







q = 3
J = 5
S = 4
X4 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta4 = rep(0, (q*q*q + 5*q)/6)
beta4_2 = create_random_beta(q)

X4_opt = mixture_coord_ex_mnl(
  X = X4, 
  beta = beta4, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)


X4_2_opt = mixture_coord_ex_mnl(
  X = X4, 
  beta = beta4_2$beta, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)







q = 5
J = 15
S = 10
X5 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta5 = rep(0, (q*q*q + 5*q)/6)
beta5_2 = create_random_beta(q)

X5_opt = mixture_coord_ex_mnl(
  X = X5, 
  beta = beta5, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)


Sys.time()
X5_2_opt = mixture_coord_ex_mnl(
  X = X5, 
  beta = beta5_2$beta, 
  n_cox_points = 1000, 
  max_it = 5, 
  verbose = 1, 
  plot_designs = T
)
Sys.time()







