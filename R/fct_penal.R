#' Compute the penalty portion of the objective function
#'
#' @inheritParams spatial_clust
#'
#' @return
#' @export
#'
#' @examples
fct_penal <- function(u, w, j, cores){
  return(2*sum(diag(t(w)%*%((u^2)%*%j - u%*%t(u)))))
}