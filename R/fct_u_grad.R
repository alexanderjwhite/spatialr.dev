#' Internal U Gradient Function
#'
#' @inheritParams fct_optimize
#' @param uv_exp matrix; stored matrix multiplication
#' @param j2 matrix; matrix of 1s with dimension r by n
#' @param one matrix; matrix of 1s with same shape as u
#'
#' @return matrix; u gradient
#' 
#' @export
#' 
#' @examples
fct_u_grad <- function(x, u, v, u_penal, uv_exp, w, j2, one, lambda, cores){
    (x-uv_exp)%*%v-lambda*(2*(w%*%t(j2))*u+(2*t(w)%*%t(j2))*u-t(w)%*%u - w%*%u) - u_penal
}
