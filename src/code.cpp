// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP exp_uv(const Eigen::Map<Eigen::MatrixXd> u,
            Eigen::Map<Eigen::MatrixXd> v, 
            int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = (u * v.transpose()).array().exp();
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP grad_v(const Eigen::Map<Eigen::MatrixXd> x,
            Eigen::Map<Eigen::MatrixXd> u, 
            Eigen::Map<Eigen::MatrixXd> v,
            Eigen::Map<Eigen::MatrixXd> uv_exp,
            int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = ((x - uv_exp).transpose()) * u;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP grad_u(const Eigen::Map<Eigen::MatrixXd> x,
            Eigen::Map<Eigen::MatrixXd> u, 
            Eigen::Map<Eigen::MatrixXd> v,
            Eigen::Map<Eigen::MatrixXd> uv_exp,
            Eigen::Map<Eigen::MatrixXd> w,
            Eigen::Map<Eigen::MatrixXd> j,
            Eigen::Map<Eigen::MatrixXd> one,
            double lambda,
            int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = ((x - uv_exp) * v) - 2*lambda*(((((j * w.transpose()) * u).diagonal()).sum())*one - (w.transpose() * u) - (w * u));
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP lik_c(const Eigen::Map<Eigen::MatrixXd> x,
                 Eigen::Map<Eigen::MatrixXd> u, 
                 Eigen::Map<Eigen::MatrixXd> v,
                 int n_cores){
  
  Eigen::setNbThreads(n_cores);
  double L = ((u * v.transpose()) * x).sum() - ((u * v.transpose()).array().exp()).sum();
  return Rcpp::wrap(L);
}

// [[Rcpp::export]]
SEXP penal_c(const Eigen::Map<Eigen::MatrixXd> u, 
                 Eigen::Map<Eigen::MatrixXd> w,
                 Eigen::Map<Eigen::MatrixXd> j,
                 int n_cores){
  
  Eigen::setNbThreads(n_cores);
  double P = 2*((w.transpose()*(((u.cwiseProduct(u)) * j) - u*(u.transpose()))).diagonal()).sum();
  return Rcpp::wrap(P);
}