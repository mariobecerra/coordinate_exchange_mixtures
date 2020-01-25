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
double getScheffeLogDEfficiency(NumericMatrix X, int order){
  arma::mat X_m = getScheffe(X, order);
  arma::mat X_mT = trans(X_m);
  arma::mat XtX = X_mT * X_m;
  cx_double log_det_X_m = arma::log_det(XtX);
  double log_D_eff;
  
  // If there is an imaginary part, then the output is -Inf
  if(log_det_X_m.imag() != 0){
    log_D_eff = -10000000.0;
  } else{
    log_D_eff = log_det_X_m.real()/X_m.n_cols - log(X_m.n_rows);
  }
  
  return log_D_eff;
}



// [[Rcpp::export]]
NumericMatrix findBestCoxDirOld(NumericMatrix cox_dir, NumericMatrix X_in, int k, int order, double log_d_eff_best) {
  NumericMatrix X = clone(X_in);
  NumericMatrix X_new = clone(X);
  NumericVector cox_dir_j;
  double log_d_eff_j;
  int n_cox_points = cox_dir.nrow();
  for(int j = 0; j < n_cox_points; j++){
    X_new = clone(X);
    cox_dir_j = cox_dir(j, _);
    X_new(k-1,_) = cox_dir_j;

    log_d_eff_j = getScheffeLogDEfficiency(X_new, order);

    // If new D-efficiency is better, then keep the new one
    if(log_d_eff_j > log_d_eff_best) {
      log_d_eff_best = log_d_eff_j;
      X = clone(X_new);
    }
  }

  return X;
}



// [[Rcpp::export]]
NumericMatrix findBestCoxDir(NumericMatrix cox_dir, NumericMatrix X_in, int k, int order, double log_d_eff_best) {
  NumericMatrix X = clone(X_in);
  NumericVector x_k(X.ncol());
  
  // numeric vector for Cox direction in j-th row
  NumericVector cox_dir_j(cox_dir.ncol());
  double log_d_eff_j;
  int n_cox_points = cox_dir.nrow();
  for(int j = 0; j < n_cox_points; j++){
    x_k = X(k-1,_);
    X(k-1,_) = cox_dir(j, _);
    
    log_d_eff_j = getScheffeLogDEfficiency(X, order);
    
    // If new D-efficiency is better, then keep the new one.
    // If it's not, keep the old one.
    if(log_d_eff_j > log_d_eff_best) {
      log_d_eff_best = log_d_eff_j;
    } else{
      X(k-1,_) = x_k;
    }
  }
  
  return X;
}



// [[Rcpp::export]]
arma::vec prueba_vec(int n){
  vec out(n);
  for(int i = 0; i < n; i++){
    out[i] = 4;
  }
  return(out);
}


// [[Rcpp::export]]
NumericVector prueba_vec2(int n){
  vec out = prueba_vec(n);
  
  return(NumericVector(out.begin(), out.end()));
}


// [[Rcpp::export]]
arma::vec seq0(int init, int end, int n_points){
  vec seq_points = linspace<vec>(init, end, n_points);
  return(seq_points);
}



// [[Rcpp::export]]
arma::vec zerosCpp(int n){
  vec out = zeros(n);
  return(out);
}



// [[Rcpp::export]]
arma::vec removeElement(vec x, int ix){
  int n = x.n_elem;
  vec out;
  if(ix != 0 & ix != n-1){
    out = join_vert(x.subvec(0, ix-1), x.subvec(ix + 1, n-1)); 
  } else{
    if(ix == 0){
      out = x.subvec(1, n-1); 
    } else{
      out = x.subvec(0, n-2); 
    }
  }
  
  return(out);
}



// // [[Rcpp::export]]
// Rcpp::NumericVector armaSetDiff(arma::uvec& x, arma::uvec& y){
//   // https://stackoverflow.com/questions/29724083
//   x = arma::unique(x);
//   y = arma::unique(y);
//   
//   for (size_t j = 0; j < y.n_elem; j++) {
//     arma::uvec q1 = arma::find(x == y[j]);
//     if (!q1.empty()) {
//       x.shed_row(q1(0));
//     }
//   }
//   
//   Rcpp::NumericVector x2 = Rcpp::wrap(x);
//   x2.attr("dim") = R_NilValue;
//   return x2;
// }



// [[Rcpp::export]]
arma::mat computeCoxDirection(NumericVector x, int comp, int n_points){
   
  int i = comp - 1;
  int q = x.length();
   
  //   # points in Cox's direction
  //   seq_points = seq(from = 0, to = 1, length.out = n_points)
  vec seq_points = linspace<vec>(0, 1, n_points);
  
  //   diffs = seq_points - x[i]
  vec diffs = seq_points - x(i);
  
  // ix1 = which.min(abs(diffs))
  int ix1 = index_min(abs(diffs));
  
  //   cox_direction_aux = seq_points - diffs[ix1]
  //   cox_direction_aux2 = cox_direction_aux[-ix1]
  //   betw_0_1 = cox_direction_aux2 >= 0 & cox_direction_aux2 <= 1
  //   cox_direction_aux3 = c(0, cox_direction_aux2[betw_0_1], 1)
  vec cox_direction_aux = seq_points - diffs(ix1);

  // vec cox_direction_aux2(cox_direction_aux.n_elem - 1);
  // // copy all elements of cox_direction_aux except the one with index ix1
  // // the equivalent in R: cox_direction_aux2 = cox_direction_aux[-ix1]
  // // could not find a simpler way to do it with Armadillo
  // int aux_int = 0;
  // for(int j = 0; j < cox_direction_aux.n_elem - 1; j++){
  //   if(j != ix1) cox_direction_aux2(j) = cox_direction_aux(aux_int);
  //   aux_int++;
  // }
  vec cox_direction_aux2 = removeElement(cox_direction_aux, ix1);
  
  uvec bigger_zero = find(cox_direction_aux2 >= 0);
  uvec less_one = find(cox_direction_aux2 <= 1);
  uvec betw_0_1 = intersect(bigger_zero, less_one);
  vec cox_direction_aux3 = join_vert(zeros(1), cox_direction_aux2.elem(betw_0_1), ones(1));
  
  // cox_direction = matrix(rep(NA_real_, q*length(cox_direction_aux3)), ncol = q)
  arma::mat cox_direction = arma::mat(cox_direction_aux3.n_elem, q, fill::randu);
  


  // cox_direction[,i] = cox_direction_aux3
  cox_direction.col(i) = cox_direction_aux3;
  
  
  // deltas = cox_direction_aux3 - x[i]
  vec deltas = cox_direction_aux3 - x(i);

    for(int n = 0; n < cox_direction_aux3.n_elem; n++){
      // recompute proportions:
      vec setDiff_aux = linspace<vec>(0, q-1, q);
      vec setDiff = removeElement(setDiff_aux, i);
      int j;
      double res;
      for(int j_aux = 0; j_aux < setDiff.n_elem; j_aux++){
        j = setDiff(j_aux);
        // In case it's a corner case, i.e., x[i] = 1
        if(abs(1 - x(i)) < 1e-16) res = (1 - cox_direction(n, i))/(q-1);
        else{
          res = x(j) - deltas(n)*x(j)/(1 - x(i));
        }
        cox_direction(n, j) = res;
        j++;
      } // end j

      // if(any(cox_direction[n, ] < -1e-10 | cox_direction[n, ] > 1 + 1e10)) {
      //   stop("Error while computing Cox direction. ",
      //        "Value out of bounds.\n",
      //        "Cox direction computed:\n\tc(",
      //        paste(cox_direction[n, ], collapse = ", "), ")")
      // }
    }
    // cox_direction = unique(cox_direction)


//   for(n in seq_along(cox_direction_aux3)){
//     # recompute proportions:
//     for(j in setdiff(1:q, i)){
//       # In case it's a corner case, i.e., x[i] = 1
//       if(abs(1 - x[i]) < 1e-16) res = (1 - cox_direction[n, i])/(q-1)
//       else{
//         res = x[j] - deltas[n]*x[j]/(1 - x[i])
//       }
//       cox_direction[n, j] = res
//     } # end j
//     
//     if(any(cox_direction[n, ] < -1e-10 | cox_direction[n, ] > 1 + 1e10)) {
//       stop("Error while computing Cox direction. ",
//            "Value out of bounds.\n", 
//            "Cox direction computed:\n\tc(", 
//            paste(cox_direction[n, ], collapse = ", "), ")")
//     }
//   }
//   cox_direction = unique(cox_direction)
//   
  return(cox_direction);
}




// compute_cox_direction = function(x, comp, n_points = 11){
//   
//   i = comp
//   q = length(x)
//   
//   
//   # points in Cox's direction
//   seq_points = seq(from = 0, to = 1, length.out = n_points)
//   
//   diffs = seq_points - x[i]
//   ix1 = which.min(abs(diffs))
//   
//   cox_direction_aux = seq_points - diffs[ix1]
//   cox_direction_aux2 = cox_direction_aux[-ix1]
//   betw_0_1 = cox_direction_aux2 >= 0 & cox_direction_aux2 <= 1
//   cox_direction_aux3 = c(0, cox_direction_aux2[betw_0_1], 1)
//   
//   cox_direction = matrix(rep(NA_real_, q*length(cox_direction_aux3)), ncol = q)
//   cox_direction[,i] = cox_direction_aux3
//   deltas = cox_direction_aux3 - x[i]
//   
//   for(n in seq_along(cox_direction_aux3)){
//     # recompute proportions:
//     for(j in setdiff(1:q, i)){
//       # In case it's a corner case, i.e., x[i] = 1
//       if(abs(1 - x[i]) < 1e-16) res = (1 - cox_direction[n, i])/(q-1)
//       else{
//         res = x[j] - deltas[n]*x[j]/(1 - x[i])
//       }
//       cox_direction[n, j] = res
//     } # end j
//     
//     if(any(cox_direction[n, ] < -1e-10 | cox_direction[n, ] > 1 + 1e10)) {
//       stop("Error while computing Cox direction. ",
//            "Value out of bounds.\n", 
//            "Cox direction computed:\n\tc(", 
//            paste(cox_direction[n, ], collapse = ", "), ")")
//     }
//   }
//   cox_direction = unique(cox_direction)
//   
//   return(cox_direction)
// }





// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// timesTwo(42)
// */
