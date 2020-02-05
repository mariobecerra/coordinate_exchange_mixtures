# Optimal D-optimality for multinomial logit Scheffé model for mixtures

library(Rcpp)
library(tidyverse)

sourceCpp("utils_MNL.cpp")



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


q = 3
J = 4
S = 2
X = create_random_initial_MNL_design(q, J, S, seed = 4)
beta = create_random_beta(q)

sapply(1:4, function(i) getUjs(X, beta$beta, beta$beta_ix, i, 1))
getUs(X, beta$beta, 1) %>% as.numeric()

getInformationMatrix(X, beta$beta)
getLogDEfficiency(X, beta$beta)

aaa = findBestCoxDir(computeCoxDirection(X[, 1,1], 1, 10), X, beta$beta, 1, 1, 100000)

bbbb = mixtureCoordinateExchangeMixtureMNL(
  X, beta$beta, 100, 2, 5
)




