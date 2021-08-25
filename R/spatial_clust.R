#' Main Spatial Clustering Method
#'
#' @param x matrix; data matrix
#' @param u_init matrix; u initialization matrix
#' @param v_init matrix; v initialization matrix
#' @param coords matrix or data frame; euclidean coordinates in separate columns
#' @param lambda numeric; penalization parameter
#' @param w matrix; distance matrix
#' @param optimizer internal function; only fct_amsgrad for now
#' @param epsilon numeric; convergence criterion
#' @param max_iter integer; maximum number of iterations
#' @param verbose TRUE or FALSE; print to screen?
#' @param fast TRUE or FALSE; use compiled c?
#' @param cores integer; number of cores to use
#'
#' @return list including u, v and convergence information
#' @export
#'
#' @examples 
spatial_clust <- function(x, u_init, v_init, coords, lambda, w = NULL, optimizer = "fct_amsgrad", epsilon = 1e-8, max_iter = 1e3, verbose = TRUE, fast = TRUE, cores = 1){
  
  # ToDo: Check X matrix for appropriate conditions
  if(is.null(w)){
    w <- fct_dist_matrix(coords, verbose = TRUE)
  }
  
  result <- fct_optimize(x, u_init, v_init, w, lambda, optimizer = optimizer, epsilon = epsilon, max_iter = max_iter, verbose = verbose, fast = fast, cores = cores)
  
  return(result)
}