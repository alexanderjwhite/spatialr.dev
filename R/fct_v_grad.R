#' Internal V Gradient Function
#'
#' @inheritParams fct_optimize
#' @param uv_exp matrix; stored matrix multiplication
#' 
#' @return matrix; v gradient
#' @export
#'
#' @examples
fct_v_grad <- function(x, u, v, v_penal, uv_exp, cores){
    t(x-uv_exp)%*%u - v_penal
}
