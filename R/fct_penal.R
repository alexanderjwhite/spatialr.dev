#' Compute the penalty portion of the objective function
#'
#' @inheritParams fct_optimize
#' @param j2 matrix; matrix of 1s with dimension r by n
#'
#' @return penalization value
#' @export
#'
#' @examples
fct_penal <- function(u, w, j2, cores){
  return(sum(diag(t(w)%*%((u^2)%*%j2 + t(j2)%*%(t(u)^2) - u%*%t(u)))))
}