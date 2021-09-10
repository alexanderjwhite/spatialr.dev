#' Internal Objective Function
#'
#' @inheritParams fct_optimize
#'
#' @return numeric; objective value
#' @export
#'
#' @examples
fct_objective <- function(x, u, v, w, j, lambda, cores){
  
  lik <- sum((u %*% t(v)) %*% x) - sum(exp(u %*% t(v)))
  penal <- 2*lambda*sum(diag(t(w)%*%((u^2)%*%j - u%*%t(u))))
  
  return(lik-penal)
}