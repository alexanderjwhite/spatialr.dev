#' Title
#'
#' @param x t
#' @param u t
#' @param v t
#' @param w t
#' @param j t
#' @param lambda t
#' @param cores t
#'
#' @return t
#'
#' @useDynLib spatialr.dev, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import RcppEigen
#' 
#' @export
#' 
#' @examples
fct_u_grad_c <- function(x, u, v, w, j, lambda, cores){
  mat_mul_c((t(x)-exp(mat_mul_c(u, t(v), cores))), v, cores)-2*lambda*(sum(diag(mat_mul_c(mat_mul_c(j, t(w), cores), u, cores)))- mat_mul_c(t(w), u, cores) - mat_mul_c(w, u, cores))
}
