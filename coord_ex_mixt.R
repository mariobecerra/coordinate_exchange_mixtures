library(tidyverse)
library(ggtern)

source("utils.R")

x_in = c(0.4725440, 0.2627838, 0.2646722)
plot_cox_direction(x_in, 1)
plot_cox_direction(x_in, 1:3)

plot_cox_direction(c(1, 0, 0), 1:3)
plot_cox_direction(c(0, 1, 0), 1:3)
plot_cox_direction(c(0, 0, 1), 1:3)


# plot several Cox's directions for random designs
# Don't understand why the plots are not right when arranged in a grid with gridExtra
# I need to use ggtern's own grid.arrange
# https://stackoverflow.com/questions/42825983/issue-when-using-grid-arrange-with-ggtern-generated-plots
# do.call("grid.arrange", c(cox_direction_plots, ncol = 4)) # This works with ggtern::grid.arrange
random_designs = create_random_initial_design(8, q = 3)

cox_direction_plots = lapply(1:nrow(random_designs), function(i){
  return(plot_cox_direction(random_designs[i,], 1:3))
})

ggtern::grid.arrange(grobs = cox_direction_plots, ncol = 4)



# First degree
res_alg_order_1 = coord_ex_mixt(9, q = 3, n_cox_points = 1000, order = 1, plot_designs = T)


# Second degree
res_alg_order_2 = coord_ex_mixt(9, q = 3, n_cox_points = 1000, order = 2, plot_designs = T)


# Third degree
res_alg_order_3_1 = coord_ex_mixt(12, q = 3, n_cox_points = 100, order = 3, plot_designs = T)
res_alg_order_3_2 = coord_ex_mixt(12, q = 3, n_cox_points = 1000, order = 3, plot_designs = T)











