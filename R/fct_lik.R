#' Compute the likelihood portion of the objective function
#'
#' @inheritParams spatial_clust
#'
#' @return likelihood value
#' @export
#'
#' @examples
fct_lik <- function(x, u, v,cores){
  return(sum(x%*%(u %*% t(v))) - sum(exp(u %*% t(v))))
}