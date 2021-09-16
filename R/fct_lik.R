#' Compute the likelihood portion of the objective function
#'
#' @inheritParams fct_optimize
#'
#' @return likelihood value
#' @export
#'
#' @examples
fct_lik <- function(x, u, v,cores){
  return(sum((u %*% t(v))%*%t(x)) - sum(exp(u %*% t(v))))
}