#' Internal normalization function
#'
#' @inheritParams fct_optimize
#' @param comp matrix; component to be regularized
#' @param ident matrix; identity matrix r by r
#'
#' @return normalization
#' @export
#'
#' @examples
fct_norm <- function(comp, eta, ident, cores){
  ctc <- t(comp)%*%comp - ident
  eta*sum(diag(t(ctc)%*%ctc))
}