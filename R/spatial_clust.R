#' Main Spatial Clustering Method
#'
#' @param x matrix; data matrix
#' @param u_init matrix; u initialization matrix
#' @param v_init matrix; v initialization matrix
#' @param coords matrix or data frame; euclidean coordinates in separate columns
#' @param lambda numeric; penalization parameter. set to NULL for lambda selection.
#' @param grid integer; lamabda search space size. used only if lambda = NULL.
#' @param w matrix; distance matrix. If null, computed on the fly. 
#' @param distance the distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given. 
#' @param method string; "dist" by default uses \eqn{w_{ij} = \exp(-\alpha \times distance_{ij})} in the w computation. Other options include "knn_1" which uses k nearest neighbors and uses neighbors as 1 and non-neighbors as 0. "knn_2" uses \eqn{w_{ij} = \exp(-\alpha \times distance_{ij})} for neighbors and 0 for non-neighbors.
#' @param k integer; if knn = TRUE, number of nearest neighbors to consider
#' @param alpha numeric; parameter to adjust distance calculation
#' @param optimizer internal optimization function; options include: fct_opt_amsgrad, fct_opt_adadelt, fct_opt_adam, fct_opt_adamax, fct_opt_grad_desc, fct_opt_nadam, fct_opt_rmsprop.
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
spatial_clust <- function(x, u_init, v_init, coords, lambda = NULL, grid = 5, w = NULL, distance = "euclidean", method = "dist", k = NULL, alpha = 1, optimizer = fct_opt_amsgrad, epsilon = 1e-8, max_iter = 1e3, verbose = TRUE, fast = TRUE, cores = 1){
  
  # ToDo: Check X matrix for appropriate conditions
  
  if(is.null(w)){
    w <- fct_dist_matrix(coords, distance = distance, method = method, k = k, alpha = alpha, verbose = verbose)
  }
  
  if(is.null(lambda)){
    if(verbose){print("No lambda specified, computing appropriate range...")}
    result <- NULL
    base_run <- fct_optimize(x, u_init, v_init, w, lambda = 0, optimizer = optimizer, epsilon = epsilon, max_iter = max_iter, verbose = verbose, fast = fast, cores = cores)
    r0 <- base_run$r
    lo <- 0.001*r0
    hi <- 100*r0
    if(verbose){print(paste("Done. r0 = ", round(r0, digits = 5)))}
    lambda_seq <- seq(lo, hi, length.out = grid)
    if(verbose){print(paste("Ranging lambda from",round(lo, 5),"to", round(hi, 5), "with a grid of", grid))}
    count <- 1
    for(lambda in abs(lambda_seq)){
      if(verbose){print(paste("Search",count))}
      count <- count + 1
      res <- fct_optimize(x, u_init, v_init, w, lambda, optimizer = optimizer, epsilon = epsilon, max_iter = max_iter, verbose = verbose, fast = fast, cores = cores)
      result <- append(list(res),result)
    }
  } else {
    result <- fct_optimize(x, u_init, v_init, w, lambda, optimizer = optimizer, epsilon = epsilon, max_iter = max_iter, verbose = verbose, fast = fast, cores = cores)
  }
  
  
  
  return(result)
}