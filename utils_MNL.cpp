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
arma::mat prueba(arma::mat X){
  X.print();
  return X;
}



