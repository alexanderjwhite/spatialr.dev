#' Internal Optimization Function
#'
#' @inheritParams spatial_clust
#'
#' @return list of results
#' @export
#' 
#' @useDynLib spatialr.dev, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#'
#' @examples 
fct_optimize <- function(x, u, v, w, lambda, optimizer, epsilon, max_iter, verbose, fast, cores){
  
  if(optimizer == "fct_amsgrad"){optimize <- fct_amsgrad}
  if(fast){
    f_u_grad <- grad_u
    f_v_grad <- grad_v
    f_ojb <- objective_c
    f_exp_uv <- exp_uv
  } else {
    f_u_grad <- fct_u_grad
    f_v_grad <- fct_v_grad
    f_ojb <- fct_objective
    f_exp_uv <- fct_exp_uv
  }
   
  
  if(verbose){print("Optimizing...")}
  iter <- 0
  objective <- -Inf
  converged <- FALSE
  u_state <- NULL
  v_state <- NULL
  
  j <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  u_iter <- u
  v_iter <- v
  one <- matrix(1, nrow = nrow(u), ncol = ncol(u))
  
  while(!converged){
    iter <- iter + 1
    u_prev <- u_iter
    v_prev <- v_iter
    uv_exp <- f_exp_uv(u, v, cores)
    u_gradient <- f_u_grad(x, u, v, uv_exp, w, j, one, lambda, cores)
    v_gradient <- f_v_grad(x, u, v, uv_exp, cores)
    
    u_grad_desc <- optimize(u_gradient, u_state)
    v_grad_desc <- optimize(v_gradient, v_state)
    
    u_state <- u_grad_desc[["state"]]
    v_state <- v_grad_desc[["state"]]
    
    u_update <- u_grad_desc[["updates"]]
    v_update <- v_grad_desc[["updates"]]
    
    u <- u_prev + u_update
    v <- v_prev + v_update
    
    
    objective_prev <- objective
    objective <- f_ojb(x, u, v, cores)
    diff <- abs(objective - objective_prev)/abs(objective_prev)
    
    if(is.nan(diff)){diff <- Inf}
    if(verbose & ((iter %% 10) == 0)){print(paste0("Objective convergence at iteration ", iter, ": ", round(diff, digits = 10)))}
    if((diff < epsilon) | (iter >= max_iter)){converged <- TRUE}

  }
  if(verbose & iter >= max_iter){print("Maximum iterations reached, consider increasing.")}
  if(verbose){print("Done")}
  return(list(u=u, v=v, convergence = diff, iter = iter))
}
