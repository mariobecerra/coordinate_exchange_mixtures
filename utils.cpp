// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/



// [[Rcpp::export]]
NumericMatrix getScheffeOrder2(NumericMatrix X){
  int q = X.ncol();
  int n = X.nrow();
  int n_col_X_m = q + (q-1)*q/2;
  NumericMatrix X_m(n, n_col_X_m);
  
  // Copy X matrix into first q columns of X_m 
  for(int i = 0; i < n; i++){
    for(int j = 0; j < q; j++){
      X_m(i,j) = X(i,j);
    }
  }
  
  // Fill rest of matrix column-wise
  int k = q-1;
  for(int i = 0; i < (q-1); i++){
    for(int j = i+1; j < q; j++){
      k = k+1;
      X_m(_,k) = X(_,i)*X(_,j);
    }
  }
  
  // element-wise
  // int k = q-1;
  // for(int i = 0; i < (q-1); i++){
  //   for(int j = i+1; j < q; j++){
  //     k = k+1;
  //     for(int l = 0; l < n; l++){
  //       X_m(l,k) = X(l,i)*X(l,j);
  //     }
  //   }  
  // }
  
  return X_m;
}



// [[Rcpp::export]]
NumericMatrix getScheffeOrder3(NumericMatrix X){
  int q = X.ncol();
  int n = X.nrow();
  
  NumericMatrix X_ord_2 = getScheffeOrder2(X);
  
  int n_col_X_ord_2 = X_ord_2.ncol();
  
  int n_col_X_m = X_ord_2.ncol();
  // compute number of columns in X_m
  // There's a formula to find this number, but I'll work it out later.
  for(int i = 0; i < (q-2); i++){
    for(int j = (i+1); j < (q-1); j++){
      for(int k = (j+1); k < q; k++){
        n_col_X_m++;
      }  
    }
  }
  
  NumericMatrix X_m(n, n_col_X_m);
  
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
      X_m(_,l) = X(_,i)*X(_,j);
    }
  }
  
  l = X_ord_2.ncol() - 1;
  for(int i = 0; i < (q - 2); i++){
    for(int j = i + 1; j < (q - 1); j++){
      for(int k = j + 1; k < q; k++){
        l++;
        X_m(_,l) = X(_,i)*X(_,j)*X(_,k);
      }  
    }
  }

  
  return X_m;
}








// [[Rcpp::export]]
arma::mat getScheffe(NumericMatrix X, int order){
  
  // not the most elegant, but works
  bool flag = (order != 1 & order != 2 & order != 3);
  
  if(flag){
    stop("Inadmissible value for order. Must be 1, 2 or 3");
  }
  
  NumericMatrix X_m;
  
  if(order == 1)  X_m = X;
  else{
    if(order == 2){
      X_m = getScheffeOrder2(X);
    } else {
      X_m = getScheffeOrder3(X);
    }
  }
  // arma::mat_out = as<mat>(X_m);
  return as<mat>(X_m);
}




// [[Rcpp::export]]
NumericMatrix transposeNumeric(const NumericMatrix & x) {
  return transpose(x);
}


// [[Rcpp::export]]
cx_double getScheffeLogDEfficiency(NumericMatrix X, int order){
  arma::mat X_m = getScheffe(X, order);
  arma::mat X_mT = trans(X_m);
  arma::mat XtX = X_mT * X_m;
  cx_double log_det_X_m = arma::log_det(XtX);
  cx_double denom = cx_double(X_m.n_cols);
  cx_double rhs = cx_double(log(X_m.n_rows));
  cx_double log_D_eff = log_det_X_m/denom - rhs;
  return log_D_eff;
}



// // [[Rcpp::export]]
// NumericMatrix findBestCoxDir(NumericMatrix cox_dir, NumericMatrix X, int i, int k, int j, int order) {
//   NumericMatrix X_new;
//   NumericVector cox_dir_j;
//   double log_d_eff_j;
//   int n_cox_points = cox_dir.nrow();
//   for(int j = 0; j < n_cox_points; j++){
//     X_new = clone(X);
//     cox_dir_j = cox_dir(j, _);
//     X_new(k,_) = cox_dir_j;
// 
//     log_d_eff_j = get_scheffe_log_D_efficiency(X_new, order = order);
//   }
// 
// 
// 
// 
//   return X;
// }



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
