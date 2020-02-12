// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


////////////////////////////////////////
// Auxiliary functions
////////////////////////////////////////


arma::vec removeElement(vec x, int ix){
  // Auxiliary function.
  // Given a vector x, returns a vector out without the element with index ix.
  // Note: this is 0-based indexing like in C++ and not 1-based indexing like in R
  int n = x.n_elem;
  vec out;
  // If ix is not the first or last element
  if(ix != 0 & ix != n-1){
    out = join_vert(x.subvec(0, ix-1), x.subvec(ix + 1, n-1)); 
  } else{
    if(ix == 0){
      // If ix is the first element
      out = x.subvec(1, n-1); 
    } else{
      // If ix is the last element
      out = x.subvec(0, n-2); 
    }
  }
  return(out);
}



// [[Rcpp::export]]
arma::mat computeCoxDirection(arma::vec& x, int comp, int n_points, int verbose){
  // Function that returns a discretization of the Cox direction for vector x in component comp.
  // Returns a matrix of dimension (cox_dir_n_elems, q), where cox_dir_n_elems is roughly equal to n_points and q is the length of x.
  // Input: 
  //     x: q dimensional vector of proportions. Must sum up to 1. 
  //        There's no input check of this because this operation is done many times in the coordinate exchange algorithm and because the x vector is fed by another function.
  //     comp: component in which the Cox direction is computed. Must be an integer between 1 and q.
  //     n_points: Number of points to use in the discretization.
  //     verbose: integer that expresses the level of verbosity. Mainly used in other functions and too much useful by itself.
  
  if(verbose >= 2) Rcout << "Computing Cox direction" << std::endl;
  
  int i = comp - 1;
  int q = x.n_elem;
  
  //   points in Cox's direction
  vec seq_points = linspace<vec>(0, 1, n_points);
  
  vec diffs = seq_points - x(i);
  
  int ix1 = index_min(abs(diffs));
  
  // Auxiliary vector 1
  vec cox_direction_aux_1 = seq_points - diffs(ix1);
  
  // Auxiliary vector 2
  // copy all elements of cox_direction_aux_1 except the one with index ix1
  // the equivalent in R: cox_direction_aux_2 = cox_direction_aux_1[-ix1]
  vec cox_direction_aux_2 = removeElement(cox_direction_aux_1, ix1);
  
  uvec bigger_zero = find(cox_direction_aux_2 >= 0); // indices of elements bigger than 0
  uvec less_one = find(cox_direction_aux_2 <= 1); // indices of elements less than 1
  uvec betw_0_1 = intersect(bigger_zero, less_one); // indices of elements between 0 and 1
  
  // Auxiliary vector 3
  vec cox_direction_aux_3 = join_vert(zeros(1), cox_direction_aux_2.elem(betw_0_1), ones(1));
  
  // Number of elements in final vector
  int cox_dir_n_elems = cox_direction_aux_3.n_elem;
  
  // Initialize final matrix with random numbers
  arma::mat cox_direction = arma::mat(cox_dir_n_elems, q, fill::randu);
  
  // Fill corresponding column
  cox_direction.col(i) = cox_direction_aux_3;
  
  // Vector of differences
  vec deltas = cox_direction_aux_3 - x(i);
  
  // Fill matrix
  for(int n = 0; n < cox_dir_n_elems; n++){
    
    // recompute proportions:
    vec setDiff_aux = linspace<vec>(0, q-1, q);
    vec setDiff = removeElement(setDiff_aux, i);
    int j;
    double result;
    
    // Iterate over ingredient proportions
    for(int j_aux = 0; j_aux < setDiff.n_elem; j_aux++){
      j = setDiff(j_aux);
      if(abs(1 - x(i)) < 1e-16) {
        // In case x(i) is numerically 1
        result = (1 - cox_direction(n, i))/(q-1);
      } else{
        // In case x(i) is not numerically 1
        result = x(j) - deltas(n)*x(j)/(1 - x(i));
      }
      cox_direction(n, j) = result;
      j++;
    } // end for j
    
    // Check that the computed directions are not smaller than 0 or bigger than 1 because of numerical error
    if(any(cox_direction.row(n) < -1e-10 || cox_direction.row(n) > 1 + 1e10)) {
      Rcout << cox_direction << std::endl;
      stop("Error while computing Cox direction. Value out of bounds.\n");
    }
  } // end for n
  return(cox_direction);
}



////////////////////////////////////////
// Gaussian model
////////////////////////////////////////

arma::mat getScheffeGaussianOrder2(arma::mat& X){
  int q = X.n_cols;
  int n = X.n_rows;
  int n_col_X_m = q + (q-1)*q/2;
  arma::mat X_m(n, n_col_X_m);
  
  // Copy X matrix into first q columns of X_m 
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      X_m(i,j) = X(i,j);
    }
  }
  
  // Fill rest of matrix element-wise
  int k = q-1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      k = k+1;
      for(int row = 0; row < n; row++){
        X_m(row,k) = X(row,i)*X(row,j);
      }
    }
  }
  
  return X_m;
}



arma::mat getScheffeGaussianOrder3(arma::mat& X){
  int q = X.n_cols;
  int n = X.n_rows;
  
  // This probably slows down the whole thing
  arma::mat X_ord_2 = getScheffeGaussianOrder2(X);
  
  int n_col_X_ord_2 = X_ord_2.n_cols;
  
  int n_col_X_m = X_ord_2.n_cols;
  // compute number of columns in X_m
  // There's a formula to find this number, but I'll work it out later.
  for(int i = 0; i < (q-2); i++){
    for(int j = (i+1); j < (q-1); j++){
      for(int k = (j+1); k < q; k++){
        n_col_X_m++;
      }  
    }
  }
  
  arma::mat X_m(n, n_col_X_m);
  
  // Copy X_ord_2 matrix into first q columns of X_m
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n_col_X_ord_2; j++){
      X_m(i,j) = X_ord_2(i,j);
    }
  }
  
  // // Fill rest of matrix column-wise
  int l = q - 1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      l++;
      for(int row = 0; row < n; row++){
        X_m(row,l) = X(row,i)*X(row,j);
      }
    }
  }
  
  l = X_ord_2.n_cols - 1;
  for(int i = 0; i < (q - 2); i++){
    for(int j = i + 1; j < (q - 1); j++){
      for(int k = j + 1; k < q; k++){
        l++;
        for(int row = 0; row < n; row++){
          X_m(row,l) = X(row,i)*X(row,j)*X(row,k);
        }
      }  
    }
  }
  
  
  return X_m;
}





arma::mat getScheffeGaussian(arma::mat& X, int order){
  
  // not the most elegant, but works
  bool flag = (order != 1 & order != 2 & order != 3);
  
  if(flag){
    stop("Inadmissible value for order. Must be 1, 2 or 3");
  }
  
  arma::mat X_m;
  
  if(order == 1)  X_m = X;
  else{
    if(order == 2){
      X_m = getScheffeGaussianOrder2(X);
    } else {
      X_m = getScheffeGaussianOrder3(X);
    }
  }
  // arma::mat_out = as<mat>(X_m);
  // return as<mat>(X_m);
  return X_m;
}




double getLogDEfficiencyGaussian(arma::mat& X, int order){
  arma::mat X_m = getScheffeGaussian(X, order);
  arma::mat X_mT = trans(X_m);
  arma::mat I = X_mT * X_m; // Information matrix
  double log_D_eff;
  
  // Attempt to do a Cholesky decomposition on the information matrix
  arma::mat L;
  double log_det_I;
  try{
    L = chol(I);
    // Compute the determinant of information matrix using the decomposition
    log_det_I = 2*sum(log(L.diag()));
    log_D_eff = log_det_I/X_m.n_cols - log(X_m.n_rows);
  }
  catch(const std::runtime_error& e){
    // If Cholesky decomposition fails, it is likely because information matrix
    // was not numerically positive definite.
    // If this happens, it is probably because a numerical inestability.
    // The function then returns the log D efficiency as a big negative number, this 
    // way the algorithm does nothing in this iteration because  the algorithm thinks
    // there was no improvement when swapping the proportions.
    Rcout << "Information matrix:\n" << I << "\n";
    Rcout << "X:\n" << X << "\n";
    Rcout << "Error in Cholesky decomposition with message: " << e.what() << std::endl;
    log_D_eff = -10000;
    Rcout << "Returning log_D_eff = " << log_D_eff << std::endl;
  }
  
  return log_D_eff;
}




arma::mat findBestCoxDirGaussian(arma::mat& cox_dir, arma::mat& X_in, int k, int order, double log_d_eff_best) {
  arma::mat X = X_in;
  int n_col_X = X.n_cols;
  arma::vec x_k(n_col_X);
  
  double log_d_eff_j;
  int n_cox_points = cox_dir.n_rows;
  for(int j = 0; j < n_cox_points; j++){
    
    // In Rcpp: x_k = X(k-1,_);
    for(int elem = 0; elem < n_col_X; elem++){
      x_k(elem) = X(k-1, elem);  
    }
    
    // In Rcpp: X(k-1,_) = cox_dir(j, _);
    for(int elem = 0; elem < n_col_X; elem++){
      X(k-1, elem) = cox_dir(j, elem);
    }
    
    log_d_eff_j = getLogDEfficiencyGaussian(X, order);
    
    // If new D-efficiency is better, then keep the new one.
    if(log_d_eff_j > log_d_eff_best) {
      // This design has a better D-efficiency, so we keep the design and update the best value
      log_d_eff_best = log_d_eff_j;
    } else{
      // This design does not have a better D-efficiency, so we return to the old design.
      // In Rcpp: X(k-1,_) = x_k;
      for(int elem = 0; elem < n_col_X; elem++){
        X(k-1,elem) = x_k(elem);
      }
      
    }
  }
  
  return X;
}




// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeGaussian(arma::mat X_orig, int order, int n_cox_points, int max_it, int verbose){
  // Performs the coordinate exchange algorithm for a Multinomial Logit Scheffé model.
  // Based on a special cubic Scheffé model as described in Ruseckaite, et al - Bayesian D-optimal choice designs for mixtures (2017)
  // X: armadillo matrix with dimensions (n, q) where:
  //    q is the number of ingredient proportions
  //    n is the number of runs
  // n_cox_points: Number of points to use in the discretization of Cox direction.
  // max_it: Maximum number of iteration that the coordinate exchange algorithm will do.
  // verbose: level of verbosity. 6 levels, in which level prints the previous plus additional things:
  //    1: Print the log D efficiency in each iteration and a final summary
  //    2: Print the values of k, s, i, and log D efficiency in each subiteration
  //    3: Print the resulting X after each iteration, i.e., after each complete pass on the data
  //    4: Print log D efficiency for each point in the Cox direction discretization
  //    5: Print the resulting X and information matrix after each subiteration
  //    6: Print the resulting X or each point in the Cox direction discretization
  // Returns an Rcpp::List object with the following objects:
  //    X_orig: The original design. Armadillo matrix with dimensions (n, q).
  //    X: The optimized design. Armadillo matrix with dimensions (n, q).
  //    d_eff_orig: log D-efficiency of the original design.
  //    d_eff: log D-efficiency of the optimized design.
  //    n_iter: Number of iterations performed.
  
  // Does not do input checks because the R wrapper function does them.
  
  // Create new matrix, otherwise it is modified in R too
  arma::mat X = X_orig;
  
  int n_runs = X.n_rows;
  int q = X.n_cols;
  
  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);
  
  // Vector of ingredient proportions
  arma::vec x(q);
  
  
  double log_d_eff_orig = getLogDEfficiencyGaussian(X, order);
  double log_d_eff_best = log_d_eff_orig;
  double log_d_eff_aux = -1e308; // -Inf
  
  // Coordinate exchanges
  int it = 0;
  while(it < max_it){
    it = it + 1;
    if(verbose >= 1) Rcout << "Iter: " << it << ", log D-efficiency: " << log_d_eff_best << std::endl;
    
    // If there was no improvement in this iteration
    if(abs(log_d_eff_aux - log_d_eff_best) < 1e-16) break;
    
    log_d_eff_aux = log_d_eff_best;
    
    for(int k = 1; k <= n_runs; k++){
      for(int i = 0; i < q; i++){
        
        if(verbose >= 2) Rcout << "\nIter: " << it <<  ", k = " << k << ", i = " << i << std::endl;
        
        // populate x vector with corresponding ingredient proportions
        for(int l = 0; l < q; l++){
          x(l) = X(k-1, l);
        }
        
        cox_dir = computeCoxDirection(x, i+1, n_cox_points, verbose);
        X = findBestCoxDirGaussian(cox_dir, X, k, order, log_d_eff_best);
        log_d_eff_best = getLogDEfficiencyGaussian(X, order);
        
        if(verbose >= 2) Rcout << "Log D-eff: " << log_d_eff_best << std::endl;
        
        if(verbose >= 5){  
          Rcout << "X =\n" << X << std::endl;
        }
        
      } // end for i
      
    } // end for k
    
    if(verbose >= 3) Rcout << "X =\n" << X << std::endl;
    
    if(verbose >= 2) Rcout << std::endl << std::endl;
    
  } // end while
  
  if(verbose >= 1){
    Rcout << std::endl;
    Rcout << "Original log D-efficiency: " << log_d_eff_orig;
    Rcout << std::endl;
    Rcout << "Final log D-efficiency: " << log_d_eff_best;
    Rcout << std::endl;
    Rcout << "Number of iterations: " << it;
    Rcout << std::endl;
  }
  
  // return object
  return Rcpp::List::create(
    _["X_orig"] = X_orig,
    _["X"] = X,
    _["d_eff_orig"] = log_d_eff_orig,
    _["d_eff"] = log_d_eff_best,
    _["n_iter"] = it
  );
  
} // end function



////////////////////////////////////////
// Multinomial logit (MNL) model
////////////////////////////////////////







