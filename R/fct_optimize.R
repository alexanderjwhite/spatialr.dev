#' Internal Optimization Function
#'
#' @inheritParams spatial_clust
#' @param u matrix; u matrix
#' @param v matrix; v matrix
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
  
  optimize <- optimizer
  if(fast){
    f_u_grad <- grad_u
    f_v_grad <- grad_v
    f_exp_uv <- exp_uv
    f_lik <- lik_c
    f_penal <- penal_c
  } else {
    f_u_grad <- fct_u_grad
    f_v_grad <- fct_v_grad
    f_exp_uv <- fct_exp_uv
    f_lik <- fct_lik
    f_penal <- fct_penal
  }
   
  
  if(verbose){print("Optimizing...")}
  iter <- 0
  objective <- -Inf
  converged <- FALSE
  u_state <- NULL
  v_state <- NULL
  
  j1 <- matrix(1, nrow = nrow(v), ncol = nrow(u))
  j2 <- matrix(1, nrow = ncol(u), ncol = nrow(u))
  u_iter <- u
  v_iter <- v
  one <- matrix(1, nrow = nrow(u), ncol = ncol(u))
  
  while(!converged){
    iter <- iter + 1
    u_prev <- u_iter
    v_prev <- v_iter
    uv_exp <- f_exp_uv(u, v, cores)
    u_gradient <- f_u_grad(x, u, v, uv_exp, w, j2, one, lambda, cores)
    v_gradient <- f_v_grad(x, u, v, uv_exp, cores)
    
    u_grad_desc <- optimize(u_gradient, u_state)
    v_grad_desc <- optimize(v_gradient, v_state)
    
    u_state <- u_grad_desc[["state"]]
    v_state <- v_grad_desc[["state"]]
    
    u_update <- u_grad_desc[["updates"]]
    v_update <- v_grad_desc[["updates"]]
    
    u <- u_prev + u_update
    v <- v_prev + v_update
    
    q <- diag(apply(v, 2, function(y){norm(y, type="2")}))
    q_inv <- solve(q)

    u <- u %*% q_inv
    v <- v %*% q
    
    
    objective_prev <- objective
    lik <- f_lik(x, u, v, uv_exp, j1, cores)
    penal <- f_penal(u, w, j2, cores)
    r <- lik/penal
    objective <- lik-lambda*penal
    diff <- abs(objective - objective_prev)/abs(objective_prev)
    
    if(is.nan(diff)){diff <- Inf}
    if(verbose & ((iter %% 10) == 0)){print(paste0("Objective convergence | r at iteration ", iter, ": ", round(diff, digits = 10), "|lik:", round(lik, digits = 2), "|penal:", round(penal, digits = 2)))}
    if((diff < epsilon) | (iter >= max_iter)){converged <- TRUE}

  }
  if(verbose & iter >= max_iter){print("Maximum iterations reached, consider increasing.")}
  if(verbose){print("Done")}
  return(list(u=u, v=v, penal = penal, lik = lik, r = r, lambda = lambda, objective = objective, convergence = diff, iter = iter))
}
