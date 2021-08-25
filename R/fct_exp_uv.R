#' Internal uv multiplication function
#'
#' @inheritParams fct_optimize 
#'
#' @return matrix;
#' @export
#'
#' @examples
fct_exp_uv <- function(u, v, cores){
  return(exp(u%*%t(v)))
}