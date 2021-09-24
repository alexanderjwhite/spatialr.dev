#' Internal Optimization Function
#'
#' @inheritParams spatial_clust
#' @param u matrix; u matrix
#' @param v matrix; v matrix
#' @param u_penal numeric; penalization if u is normalized
#' @param v_penal numeric; penalization if v is normalized
#'
#' @return list of results
#' @export
#' 
#' @useDynLib spatialr.dev, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#'
#' @examples 
fct_optimize <- function(x, u, v, w, lambda, optimizer, epsilon, max_iter, eta, norm_comp = norm_comp, verbose, fast, cores){
  
  optimize <- optimizer
  if(fast){
    f_u_grad <- grad_u
    f_v_grad <- grad_v
    f_exp_uv <- exp_uv
    f_lik <- lik_c
    f_penal <- penal_c
    f_reg <- reg_c
    f_norm <- norm_c
  } else {
    f_u_grad <- fct_u_grad
    f_v_grad <- fct_v_grad
    f_exp_uv <- fct_exp_uv
    f_lik <- fct_lik
    f_penal <- fct_penal
    f_reg <- fct_reg
    f_norm <- fct_norm
  }
   
  
  if(verbose){print("Optimizing...")}
  iter <- 0
  objective <- -Inf
  obj_track <- NULL
  lik_track <- NULL
  penal_track <- NULL
  converged <- FALSE
  u_state <- NULL
  v_state <- NULL
  u_penal <- matrix(0, nrow = nrow(u), ncol = ncol(u))
  v_penal <- matrix(0, nrow = nrow(v), ncol = ncol(v))
  comp_normalize <- 0
  
  j1 <- matrix(1, nrow = nrow(v), ncol = nrow(u))
  j2 <- matrix(1, nrow = ncol(u), ncol = nrow(u))
  one <- matrix(1, nrow = nrow(u), ncol = ncol(u))
  ident <- diag(rep(1, ncol(u)))
  
  while(!converged){

    
    u_prev <- u
    v_prev <- v
    uv_exp <- f_exp_uv(u, v, cores)
    
    iter <- iter + 1
    
    if(norm_comp == "u"){
      u_penal <- f_reg(u, eta, cores)
      comp_normalize <- f_norm(u, eta, ident, cores)
    } else if(norm_comp == "v"){
      v_penal <- f_reg(v, eta, cores)
      comp_normalize <- f_norm(v, eta, ident, cores)
    } 
    
    u_gradient <- f_u_grad(x, u, v, u_penal, uv_exp, w, j2, one, lambda, cores)
    v_gradient <- f_v_grad(x, u, v, v_penal, uv_exp, cores)
    
    u_grad_desc <- optimize(u_gradient, u_state)
    v_grad_desc <- optimize(v_gradient, v_state)
    
    u_state <- u_grad_desc[["state"]]
    v_state <- v_grad_desc[["state"]]
    
    u_update <- u_grad_desc[["updates"]]
    v_update <- v_grad_desc[["updates"]]
    
    u <- u_prev + u_update
    v <- v_prev + v_update
    
    objective_prev <- objective
    
    lik <- f_lik(x, u, v, uv_exp, j1, cores)
    lik_track <- c(lik_track, lik)
    
    penal <- f_penal(u, w, j2, cores)
    penal_track <- c(penal_track,penal)
    
    r <- lik/penal
    
    objective <- lik-lambda*penal-comp_normalize
    obj_track <- c(obj_track,objective)
    
    diff <- abs(objective - objective_prev)/abs(objective_prev)
    
    if(is.nan(diff)){diff <- Inf}
    if(verbose & ((iter %% 10) == 0)){print(paste0("Objective at iteration ", iter, ": ", round(objective, digits = 10), "|lik:", round(lik, digits = 2), "|penal:", round(penal, digits = 2)))}
    if(verbose & ((iter %% 10) == 0)){fct_push_me(paste("lambda:",lambda,"|iter:", iter,"|diff:",round(diff, digits = 10)))}
    if((diff < epsilon) | (iter >= max_iter)){converged <- TRUE}

  }
  if(verbose & iter >= max_iter){print("Maximum iterations reached, consider increasing.")}
  if(verbose){print("Done")}
  return(list(u=u, v=v, penal = penal, lik = lik, r = r, lambda = lambda, objective = objective, convergence = diff, iter = iter))
}
