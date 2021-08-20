#' Internal Optimization Function
#'
#' @param x data matrix
#' @param u initial u matrix
#' @param v initial v matrix
#' @param w distance matrix
#' @param lambda numeric; penalization parameter
#' @param optimizer internal function; only fct_amsgrad for now
#' @param epsilon numeric; convergence criterion
#' @param max_iter integer; maximum number of iterations
#' @param verbose TRUE or FALSE; print to screen?
#'
#' @return list of results
#' @export
#'
#' @examples fct_optimize(x = matrix(rnorm(100), nrow = 10, ncol = 10), 
#' u = matrix(rnorm(20), nrow = 10, ncol = 2), 
#' v = matrix(rnorm(20), nrow = 10, ncol = 2), 
#' w = matrix(rnorm(100), nrow = 10, ncol = 10), 
#' lambda = 1, optimizer = "fct_amsgrad", epsilon = 1e-8, max_iter = 10, verbose = TRUE)
fct_optimize <- function(x, u, v, w, lambda, optimizer = "fct_amsgrad", epsilon = 1e-8, max_iter = 1e3, verbose = TRUE){
  
  if(optimizer == "fct_amsgrad"){optimize <- fct_amsgrad}
  
  if(verbose){print("Optimizing...")}
  iter <- 0
  objective <- -Inf
  converged <- FALSE
  u_state <- NULL
  v_state <- NULL
  
  j <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  u_iter <- u
  v_iter <- v
  
  while(!converged){
    iter <- iter + 1
    u_prev <- u_iter
    v_prev <- v_iter
    u_gradient <- fct_u_grad(x, u, v, w, j, lambda)
    v_gradient <- fct_v_grad(x, u, v)
    
    u_grad_desc <- optimize(u_gradient, u_state)
    v_grad_desc <- optimize(v_gradient, v_state)
    
    u_state <- u_grad_desc[["state"]]
    v_state <- v_grad_desc[["state"]]
    
    u_update <- u_grad_desc[["updates"]]
    v_update <- v_grad_desc[["updates"]]
    
    u <- u_prev + u_update
    v <- v_prev + v_update
    
    objective_prev <- objective
    objective <- sum((u %*% t(v)) %*% t(x)) - sum(exp(u %*% t(v)))
    diff <- abs(objective - objective_prev)/abs(objective_prev)
    
    if(is.nan(diff)){diff <- Inf}
    if(verbose & ((iter %% 100) == 0)){print(paste0("Objective convergence at iteration ", iter, ": ", round(diff, digits = 10)))}
    if((diff < epsilon) | (iter >= max_iter)){converged <- TRUE}

  }
  if(verbose & iter >= max_iter){print("Maximum iterations reached, consider increasing.")}
  if(verbose){print("Done")}
  return(list(u=u, v=v, convergence = diff, iter = iter))
}
