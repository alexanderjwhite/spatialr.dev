#' Internal V Gradient Function
#'
#' @inheritParams fct_optimize
#' @param uv_exp matrix; stored matrix multiplication
#' 
#' @return matrix; v gradient
#' @export
#'
#' @examples
fct_v_grad <- function(x, u, v, v_penal, eta, uv_exp, cores){
  
  if(v_penal){
    v_grad <- t(x-uv_exp)%*%u - eta*(4*v%*%t(v)%*%v - 4*v)
  } else {
    v_grad <- t(x-uv_exp)%*%u
  }
  return(v_grad)
}
