# Optimal D-optimality for multinomial logit Scheffé model for mixtures

library(Rcpp)
library(tidyverse)
library(ggtern)

source("utils_MNL.R")
sourceCpp("utils_MNL.cpp")




q = 3
J = 5
S = 3
X2 = create_random_initial_MNL_design(q, J, S, seed = 4)
beta2 = rep(0, (q*q*q + 5*q)/6)

(Xs2_1 = getXs(X2, 1))
(Xs2_2 = getXs(X2, 2))
(Xs2_3 = getXs(X2, 3))

(ps2_1 = getPs(X2, beta2, 1, Xs2_1))
(ps2_2 = getPs(X2, beta2, 2, Xs2_2))
(ps2_3 = getPs(X2, beta2, 3, Xs2_3))

middle2_1 = create_diag_mat(ps2_1) - ps2_1 %*% t(ps2_1)
middle2_2 = create_diag_mat(ps2_2) - ps2_2 %*% t(ps2_2)
middle2_3 = create_diag_mat(ps2_3) - ps2_3 %*% t(ps2_3)

(I2_1 = t(Xs2_1) %*% middle2_1 %*% Xs2_1)
(I2_2 = t(Xs2_2) %*% middle2_2 %*% Xs2_2)
(I2_3 = t(Xs2_3) %*% middle2_3 %*% Xs2_3)

I2_1 + I2_2 + I2_3

getInformationMatrix(X2, beta2)
(I2_1 + I2_2 + I2_3) == getInformationMatrix(X2, beta2)

getLogDEfficiency(X2, beta2, 1)
determinant(I2_1 + I2_2 + I2_3)

findBestCoxDir(computeCoxDirection(X2[, 1,1], 1, 10), X2, beta2, 1, 1, 100000, 1)

X2_opt = mixtureCoordinateExchangeMixtureMNL(
  X_orig = X2, 
  beta = beta2, 
  n_cox_points = 5, 
  max_it = 2, 
  verbose = 5
)

  
X2_opt2 = mnl_mixture_coord_ex(
  X = X2, 
  beta = beta2, 
  n_cox_points = 5, 
  max_it = 2, 
  verbose = 5
)










q = 3
J = 5
S = 4
X3 = create_random_initial_MNL_design(q, J, S, seed = 3)
beta3 = rep(0, (q*q*q + 5*q)/6)
beta3_2 = create_random_beta(q)

X3_opt = mnl_mixture_coord_ex(
  X = X3, 
  beta = beta3, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)


X3_2_opt = mnl_mixture_coord_ex(
  X = X3, 
  beta = beta3_2$beta, 
  n_cox_points = 100, 
  max_it = 10, 
  verbose = 1, 
  plot_designs = T
)











