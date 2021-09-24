#' Internal regularization function
#'
#' @inheritParams fct_optimize 
#' @param comp matrix; component to be regularized
#'
#' @return normalization
#' @export
#'
#' @examples
fct_reg <- function(comp, eta, core){
  eta*(4*comp%*%t(comp)%*%comp - 4*comp)
}