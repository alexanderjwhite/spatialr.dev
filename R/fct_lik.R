#' Compute the likelihood portion of the objective function
#'
#' @inheritParams fct_optimize
#' @param j1 matrix; matrix of 1s with dimension p by n
#'
#' @return likelihood value
#' @export
#'
#' @examples
fct_lik <- function(x, u, v, uv_exp, j1, cores){
  sum(diag(t(x)%*%u%*%t(v))) - sum(diag(uv_exp%*%j1))
  #return(sum((u %*% t(v))%*%t(x)) - sum(exp(u %*% t(v))))
}