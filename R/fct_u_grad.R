#' Internal U Gradient Function
#'
#' @inheritParams fct_optimize
#' @param uv_exp matrix; stored matrix multiplication
#' @param j matrix; matrix of 1s
#' @param one matrix; matrix of 1s with same shape as u
#'
#' @return matrix; u gradient
#' 
#' @export
#' 
#' @examples
fct_u_grad <- function(x, u, v, uv_exp, w, j, one, lambda, cores){
  (x-uv_exp)%*%v-2*lambda*(2*(w%*%t(j))*u-t(w)%*%u - w%*%u)
}
