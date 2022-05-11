// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List fct_c_opt_adam(arma::mat gradients,
                    List state){
  if (state.size() == 0){
    state["beta1"] = 0.9;
    state["beta2"] = 0.999;
    state["epsilon"] = 1e-8;
    state["iteration"] = 1;
    state["m"] = 0*gradients;
    state["v"] = 0*gradients;
    state["alpha"] = 1e-3;
  }
  
  double beta1 = as<double>(state["beta1"]);
  double beta2 = as<double>(state["beta2"]);
  double epsilon = as<double>(state["epsilon"]);
  int iteration = as<int>(state["iteration"]);
  arma::mat m = as<arma::mat>(state["m"]);
  arma::mat v = as<arma::mat>(state["v"]);
  double alpha = as<double>(state["alpha"]);
  
  m = beta1 * m + (1.0 - beta1)*gradients;
  v = beta2 * v + (1.0 - beta2)*(gradients % gradients);
  arma::mat mhat = m / (1.0 - pow(beta1, iteration));
  arma::mat vhat = v / (1.0 - pow(beta2, iteration));
  state["updates"] = alpha * mhat / (sqrt(vhat) + epsilon);
  
  state["iteration"] = iteration + 1;
  
  return(state);
}

// [[Rcpp::export]]
double pdist(arma::mat x, arma::sp_mat w, arma::mat index){
  double result = 0;
  arma::mat diff;
  index = index - 1;
  for (int i = 0; i < index.n_rows; i++){
    for (int j = 0; j < index.n_cols; j++){
      diff = x.col(i) - x.col(index(i,j));
      result = result + w(i,index(i,j)) * accu(diff % diff);
    }
  }
  return(result);
}

// [[Rcpp::export]]
List fct_c_optimize(arma::sp_mat x,
                    arma::mat u,
                    arma::mat v,
                    arma::sp_mat w,
                    arma::mat index,
                    double lambda,
                    double epsilon,
                    int maxiter,
                    bool display_progress=true){
  
  // Progress p(maxiter, display_progress);
  bool converged = false;
  double objective = -std::numeric_limits<double>::max();
  double objective_prev;
  double lik;
  double diff;
  double osp;
  double u_diff;
  double v_diff;
  List results;
  List u_state;
  List v_state;
  arma::mat u_prev;
  arma::mat v_prev;
  arma::mat uv_t;
  arma::mat uv_exp;
  arma::mat u_gradient;
  arma::mat v_gradient;
  arma::mat lik_store(maxiter, 1, fill::zeros);
  arma::mat obj_store(maxiter, 1, fill::zeros);
  arma::mat udiff_store(maxiter, 1, fill::zeros);
  arma::mat vdiff_store(maxiter, 1, fill::zeros);
  arma::mat osp_store(maxiter, 1, fill::zeros);
  List u_grad_desc;
  List v_grad_desc;
  arma::mat u_update;
  arma::mat v_update;
  arma::mat j;
  j.ones(u.n_cols, x.n_cols);
  int i = 1;
  
  while (!converged) {
    
    
    
    u_prev = u;
    v_prev = v;
    objective_prev = objective;
    uv_t = u * v.t();
    uv_exp = exp(uv_t);
    
    
    
    if (i % 1 == 0){
      Rcpp::checkUserInterrupt();
    }
    
    u_gradient = (x - uv_exp) * v;
    v_gradient = ((x - uv_exp).t() * u)  - lambda*((2 * (w + w.t()) * j.t()) % v - (w + w.t()) * v);
    
    
    u_grad_desc = fct_c_opt_adam(u_gradient, u_state);
    v_grad_desc = fct_c_opt_adam(v_gradient, v_state);
    
    u_state = u_grad_desc;
    v_state = v_grad_desc;
    
    u_update = as<arma::mat>(u_grad_desc["updates"]);
    v_update = as<arma::mat>(v_grad_desc["updates"]);
    
    u = u_prev + u_update;
    v = v_prev + v_update;
    
    lik = accu(uv_t % x) - accu(uv_exp) ;
    
    osp = lambda*pdist(v.t(), w, index);
    objective = lik-osp;
    diff = abs(objective - objective_prev)/abs(objective_prev);
    u_diff = accu(abs(u_prev - u))/accu(abs(u_prev));
    v_diff = accu(abs(v_prev - v))/accu(abs(v_prev));
    
    lik_store[i] = lik;
    osp_store[i] = osp;
    obj_store[i] = objective;
    udiff_store[i] = u_diff;
    vdiff_store[i] = v_diff;
    
    if (i % 1 == 0){
      Rcout << "Iteration: " << i << " | Objective: " << diff << "\n";
    }
    
    if((i >= maxiter) || ((diff < epsilon) && (u_diff < epsilon) && (v_diff < epsilon))){
      converged = true;
    }
    
    i = i + 1;
    
  }
  results["u"] = u;
  results["v"] = v;
  results["likelihood"] = lik_store;
  results["osp"] = osp_store;
  results["objective"] = obj_store;
  results["udiff"] = udiff_store;
  results["vdiff"] = vdiff_store;
  return(results);
  
}



