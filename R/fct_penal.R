#' Compute the penalty portion of the objective function
#'
#' @inheritParams fct_optimize
#' @param j matrix; matrix of 1s
#'
#' @return penalization value
#' @export
#'
#' @examples
fct_penal <- function(u, w, j, cores){
  return(2*sum(diag(t(w)%*%((u^2)%*%j - u%*%t(u)))))
}