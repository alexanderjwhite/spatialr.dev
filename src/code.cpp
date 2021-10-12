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
SEXP grad_v(const Eigen::SparseMatrix<double> x,
            Eigen::Map<Eigen::MatrixXd> u, 
            Eigen::Map<Eigen::MatrixXd> v,
            Eigen::Map<Eigen::MatrixXd> v_penal,
            Eigen::Map<Eigen::MatrixXd> uv_exp,
            int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = ((x - uv_exp).transpose()) * u - v_penal;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP grad_u(const Eigen::SparseMatrix<double> x,
            Eigen::Map<Eigen::MatrixXd> u, 
            Eigen::Map<Eigen::MatrixXd> v,
            Eigen::Map<Eigen::MatrixXd> u_penal,
            Eigen::Map<Eigen::MatrixXd> uv_exp,
            Eigen::SparseMatrix<double> w,
            Eigen::Map<Eigen::MatrixXd> j2,
            Eigen::Map<Eigen::MatrixXd> one,
            double lambda,
            int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = ((x - uv_exp) * v) - lambda*((2 * w.transpose() * j2.transpose()).cwiseProduct(u) + (2 * w.transpose() * j2.transpose()).cwiseProduct(u) - w.transpose() * u - w * u) - u_penal;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP lik_c(const Eigen::SparseMatrix<double> x,
                 Eigen::Map<Eigen::MatrixXd> u, 
                 Eigen::Map<Eigen::MatrixXd> v,
                 Eigen::Map<Eigen::MatrixXd> uv_exp,
                 Eigen::Map<Eigen::MatrixXd> j1,
                 int n_cores){
  
  Eigen::setNbThreads(n_cores);
  double L = ((x.transpose() * u * v.transpose()).diagonal()).sum() - ((uv_exp * j1).diagonal()).sum() ;
  return Rcpp::wrap(L);
}

// [[Rcpp::export]]
SEXP penal_c(const Eigen::Map<Eigen::MatrixXd> u,
             Eigen::SparseMatrix<double> w,
             Eigen::Map<Eigen::MatrixXd> j2,
                 int n_cores){
  
  Eigen::setNbThreads(n_cores);
  double P = (w.transpose() * ((u.cwiseProduct(u)) * j2 + j2.transpose() * (u.transpose().cwiseProduct(u.transpose())) - 2 * u * u.transpose())).diagonal().sum();
  return Rcpp::wrap(P);
}

// [[Rcpp::export]]
SEXP reg_c(const Eigen::Map<Eigen::MatrixXd> comp, 
           double eta,
             int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd R = eta*(4*comp*comp.transpose()*comp - 4*comp);
  return Rcpp::wrap(R);
}

// [[Rcpp::export]]
SEXP norm_c(const Eigen::Map<Eigen::MatrixXd> comp, 
           double eta,
           const Eigen::Map<Eigen::MatrixXd> ident, 
           int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd CTC = comp.transpose() * comp - ident;
  double N = eta*((CTC.transpose() * CTC).diagonal().sum());
  return Rcpp::wrap(N);
}



