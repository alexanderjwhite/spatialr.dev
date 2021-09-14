#' Compute the likelihood portion of the objective function
#'
#' @inheritParams spatial_clust
#'
#' @return likelihood value
#' @export
#'
#' @examples
fct_lik <- function(x, u, v,cores){
  return(sum((u %*% t(v)) %*% x) - sum(exp(u %*% t(v))))
}