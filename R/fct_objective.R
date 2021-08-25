#' Internal Objective Function
#'
#' @inheritParams fct_optimize
#'
#' @return numeric; objective value
#' @export
#'
#' @examples
fct_objective <- function(x, u, v, cores){
  return(sum((u %*% t(v)) %*% t(x)) - sum(exp(u %*% t(v))))
}