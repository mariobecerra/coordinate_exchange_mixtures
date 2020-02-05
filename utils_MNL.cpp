// #include <Rcpp.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double getUjs(arma::cube X, arma::vec beta, arma::mat beta_ix, int j, int s){
  // j and s are 1-indexed
  
  int q = X.n_rows;
  if((q*q*q+5*q)/6 != beta.n_elem) stop("Incompatible q in beta and X");
  
  // int J = X.n_cols;
  // int S = X.n_elem_slice;
  double Ujs;
  
  double first_term = 0;
  double second_term = 0;
  double third_term = 0;
  
  
  arma::uvec ids_na;
  arma::uvec ids_i;
  arma::uvec ids_k;
  arma::uvec ids_l;
  arma::uvec ids_iter;
  int id;
  
  
  // the first q items in vector correspond to beta_i
  for(int i = 0; i < q-1; i++){
    // subtract 1 from j and s because it is 1-indexed
    first_term = first_term + (beta(i) - beta(q-1))*X(i, j-1, s-1);
  }
  
  // second part
  double beta_ik;
  for(int i = 0; i < q-1; i++){
    for(int k = i+1; k < q; k++){
      ids_na = intersect(find_finite(beta_ix.col(1)), find_nonfinite(beta_ix.col(2)));
      ids_i = find(beta_ix.col(0) == i+1);
      ids_k = find(beta_ix.col(1) == k+1);
      // Gives the correct index in the beta vector
      ids_iter = intersect(ids_na, intersect(ids_i, ids_k));
      // extract first (and only element) of this vector
      // this is a trick so that I don't have to cast a vector to integer
      id = ids_iter(0);
      beta_ik = beta(id);
      // subtract 1 from j and s because it is 1-indexed
      second_term = second_term + beta_ik*X(i, j-1, s-1)*X(k, j-1, s-1); 
    }
  }
  
  
  // third part
  double beta_ikl;
  for(int i = 0; i < q-2; i++){
    for(int k = i+1; k < q-1; k++){
      for(int l = k+1; l < q; l++){
        ids_na = intersect(find_finite(beta_ix.col(1)), find_finite(beta_ix.col(2)));
        ids_i = find(beta_ix.col(0) == i+1);
        ids_k = find(beta_ix.col(1) == k+1);
        ids_l = find(beta_ix.col(2) == l+1);
        // Gives the correct index in the beta vector
        ids_iter = intersect(ids_na, intersect(ids_i, intersect(ids_k, ids_l)));
        // extract first (and only element) of this vector
        // this is a trick so that I don't have to cast a vector to integer
        id = ids_iter(0);
        beta_ikl = beta(id);
        // subtract 1 from j and s because it is 1-indexed
        third_term = third_term + beta_ikl*X(i, j-1, s-1)*X(k, j-1, s-1)*X(l, j-1, s-1);
      }
    }
  }
  
  Ujs = first_term + second_term + third_term;
  return Ujs;
}


// [[Rcpp::export]]
arma::mat getXs(arma::cube X, int s){
  
  int q = X.n_rows;
  int J = X.n_cols;
  int m = (q*q*q + 5*q)/6;
  
  // arma::mat Xs(J, m);
  arma::mat Xs(J, m-1, fill::zeros);
  
  int col_counter = 0;
  
  for(int i = 0; i < q-1; i++){
    for(int j = 0; j < J; j++){
      // subtract 1 from s because it is 1-indexed
      Xs(j, col_counter) = X(i, j, s-1);
    }
    col_counter++;
  }
  
  
  // second part
  for(int i = 0; i < q-1; i++){
    for(int k = i+1; k < q; k++){
      for(int j = 0; j < J; j++){
        Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1);
      }
      col_counter++;
    }
  }



  // third part
  for(int i = 0; i < q-2; i++){
    for(int k = i+1; k < q-1; k++){
      for(int l = k+1; l < q; l++){
        for(int j = 0; j < J; j++){
          Xs(j, col_counter) = X(i, j, s-1)*X(k, j, s-1)*X(l, j, s-1);
        }
        col_counter++;
      }
    }
  }
  
  // Xs.print();
  
  return Xs;
  
}


// [[Rcpp::export]]
arma::vec getUs(arma::cube X, arma::vec beta, int s, arma::mat Xs){
  int J = X.n_cols;
  int q = X.n_rows;
  
  // arma::mat Xs = getXs(X, s);
  int m = Xs.n_cols + 1;
  if(m != beta.n_elem) stop("Incompatible q in beta and X");
  
  arma::vec beta2(m-1);
  // compute beta_i_star = beta_i - beta_q
  for(int i = 0; i < q-1; i++){
    beta2(i) = beta(i) - beta(q-1);
    // cout << "beta_star_" << i << ": " << beta2(i) << "\n";
  }
  
  for(int i = q-1; i < m-1; i++){
    beta2(i) = beta(i+1);
    // cout << "beta2_" << i << ": " << beta2(i) << "\n";
  }
  
  
  arma::vec Us(J);
  
  Us = Xs*beta2;
  return Us;
}


// [[Rcpp::export]]
arma::vec getPs(arma::cube X, arma::vec beta, int s, arma::mat Xs){
  int J = X.n_cols;

  arma::vec Us(J);
  Us = getUs(X, beta, s, Xs);
  
  arma::vec exp_Ujs(J);
  arma::vec P(J);
  
  exp_Ujs = exp(Us);
  
  double sum_exp_Ujs = sum(exp_Ujs);
  
  
  P = exp_Ujs/sum_exp_Ujs;
  
  if(abs(sum(P) - 1) > 1e-10) warning("Sum may not be numerically equal to 1.");

  return P;
}





// [[Rcpp::export]]
arma::mat getInformationMatrix(arma::cube X, arma::vec beta){
  int J = X.n_cols;
  int S = X.n_elem/(X.n_cols*X.n_rows);
  int m = beta.n_elem;
  
  arma::mat Xs(J, m-1);
  arma::mat I(m-1, m-1, fill::zeros);
  arma::mat identity(J, J, fill::eye);
  arma::mat ps_ps_t(J, J);
  arma::mat middle(J, J);
  arma::vec ps;

  for(int s = 1; s <= S; s++){
    Xs = getXs(X, s);
    ps = getPs(X, beta, s, Xs);;
    ps_ps_t = ps*ps.t();
    middle = ps_ps_t;
    middle.diag() = ps_ps_t.diag() - ps;
    I = I - (Xs.t())*middle*Xs;
  }
  
  return I;

}


// [[Rcpp::export]]
double getLogDEfficiency(arma::cube X, arma::vec beta){
  
  int m = beta.n_elem;
  arma::mat I(m-1, m-1, fill::zeros);
  
  I = getInformationMatrix(X, beta);
  
  cx_double log_det_I = arma::log_det(I);
  double log_D_eff;
  
  // If there is an imaginary part, then the output is -Inf
  if(log_det_I.imag() != 0){
    log_D_eff = -10000000.0;
  } else{
    // Don't know if I should scale or just return log determinant
    log_D_eff = -log_det_I.real()/I.n_cols - log(I.n_rows);
  }
  
  return log_D_eff;
}




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





// [[Rcpp::export]]
arma::mat computeCoxDirection(arma::vec x, int comp, int n_points){
  
  int i = comp - 1;
  int q = x.n_elem;
  
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
  
  int cox_dir_n_elems = cox_direction_aux3.n_elem;
  // if repeated elements
  // if(cox_direction_aux3(1) == 0) cox_direction_aux3 = removeElement(cox_direction_aux3, 1);
  // if(cox_direction_aux3(cox_dir_n_elems-2) == 1) cox_direction_aux3 = removeElement(cox_direction_aux3, cox_dir_n_elems-2);
  
  // cox_direction = matrix(rep(NA_real_, q*length(cox_direction_aux3)), ncol = q)
  arma::mat cox_direction = arma::mat(cox_dir_n_elems, q, fill::randu);
  
  
  
  // cox_direction[,i] = cox_direction_aux3
  cox_direction.col(i) = cox_direction_aux3;
  
  
  // deltas = cox_direction_aux3 - x[i]
  vec deltas = cox_direction_aux3 - x(i);
  
  for(int n = 0; n < cox_dir_n_elems; n++){
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
      // if(res < -1e-10 | )
      cox_direction(n, j) = res;
      j++;
    } // end j
    
    if(any(cox_direction.row(n) < -1e-10 || cox_direction.row(n) > 1 + 1e10)) {
      cox_direction.print();
      stop("Error while computing Cox direction. Value out of bounds.\n");
      
    }
  }
  // cox_direction = unique(cox_direction)
  return(cox_direction);
}





// [[Rcpp::export]]
arma::cube findBestCoxDir(arma::mat cox_dir, arma::cube X_in, arma::vec beta, int k, int s, double log_d_eff_best) {

  // k: Cox direction index (1 to q)
  // s: choice set (1 to S)
  
  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_in;
  
  // // Create new cube
  // arma::cube X(size(X_in));
  // 
  // // copy elements of cube
  // X = X_in;
  
  
  // int J = X.n_cols;
  int q = X.n_rows;
  // int S = X.n_elem/(X.n_cols*X.n_rows);

  // NumericVector x_k(X.ncol());
  arma::vec x_k(q);



  
  double log_d_eff_j;
  int n_cox_points = cox_dir.n_rows;
  
  for(int j = 0; j < n_cox_points; j++){
  
    // same as:
    // X(arma::span::all, arma::span(k-1), arma::span(s-1)) = cox_dir.row(j);
    // if the previous worked
    // Basically, here we replace the j-th row of the Cox direction matrix
    // in the corresponding slice and row of X
    for(int l = 0; l < q; l++){
      X(l, k-1, s-1) = cox_dir(j, l);
    }
    
    log_d_eff_j = getLogDEfficiency(X, beta);
    
    x_k = X(arma::span::all, arma::span(k-1), arma::span(s-1));
    // x_k.print();
    cout << "j = " << j << "\n";
    cout << "d ef = "  << log_d_eff_j << "\n\n";

    // If new D-efficiency is better, then keep the new one.
    // If it's not, keep the old one.
    if(log_d_eff_j > log_d_eff_best) {
      log_d_eff_best = log_d_eff_j;
    } else{
      
      // return to old state
      for(int l = 0; l < q; l++){
        X(l, k-1, s-1) = x_k(l);
      }
      
    }
  }

  return X;
}






// [[Rcpp::export]]
Rcpp::List mixtureCoordinateExchangeMixtureMNL(arma::cube X_orig, arma::vec beta, int n_cox_points, int max_it, int verbose = 1){

  // Should add some code to check the input dimensions

  // Create new cube, otherwise it is modified in R too
  arma::cube X = X_orig;
  
  int J = X.n_cols;
  int q = X.n_rows;
  int S = X.n_elem/(X.n_cols*X.n_rows);
  
  // Create matrix with appropriate dimensions for Cox direction in each iteration
  arma::mat cox_dir(n_cox_points, q);
  
  // Vector of ingredient proportions
  arma::vec x(q);


  double log_d_eff_orig = getLogDEfficiency(X, beta);
  double log_d_eff_best = log_d_eff_orig;
  double log_d_eff_aux = 1e308; // Max double

  // Coordinate exchanges
  int it = 0;
  while(it < max_it){
    it = it + 1;
    if(verbose >= 1) cout << "Iter: " << it << ", log D-efficiency: " << log_d_eff_best, "\n";

    // If there was no improvement in this iteration
    if(abs(log_d_eff_aux - log_d_eff_best) < 1e-16) break;

    log_d_eff_aux = log_d_eff_best;

    for(int k = 1; k <= J; k++){
      if(verbose >= 2) cout << "k = " << k << "\n";

      for(int s = 1; s <= S; s++){
        if(verbose >= 2) cout << "\ts = " << s << "\n";

        for(int i = 0; i < q; i++){
          if(verbose >= 2) cout << "\t\ti = " << i << "\n";
          
          // populate x vector with corresponding ingredient proportions
          for(int l = 0; l < q; l++){
            x(l) = X(l, k-1, s-1);
          }
          
          cox_dir = computeCoxDirection(x, i+1, n_cox_points);
          X = findBestCoxDir(cox_dir, X, beta, k, s, log_d_eff_best);
          log_d_eff_best = getLogDEfficiency(X, beta);

          if(verbose >= 5) {
            // X.print(); // check later because output may not be coordinated
            cout << "\t\t\tLog D-eff:" << log_d_eff_best << "\n";
          }

        } // end for i

      } // end for s

    } // end for k

    // if(verbose >= 3) X.print();

  } // end while

  if(verbose >= 1){
    cout << "\n";
    cout << "Original log D-efficiency: " << log_d_eff_orig;
    cout << "\n";
    cout << "Final log D-efficiency: " << log_d_eff_best;
    cout << "\n";
    cout << "Number of iterations: " << it;
    cout << "\n";
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









