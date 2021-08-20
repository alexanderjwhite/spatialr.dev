#' Internal Distance Matrix Computation
#'
#' @param coords matrix or data frame; euclidean coordinates in separate columns
#' @param verbose TRUE or FALSE; print to screen?
#'
#' @return distance matrix
#' @export
#' 
#' @import stats
#'
#' @examples fct_dist_matrix(matrix(rnorm(20), nrow = 10, ncol = 2), verbose = TRUE)
fct_dist_matrix <- function(coords, verbose = TRUE){
  # ToDo: conditions for coord dataframe
  if(verbose){print("Computing Distance Matrix...")}
  w <- as.matrix(1/stats::dist(coords, upper = TRUE, diag = TRUE))
  if(verbose){print("Done")}
  return(w)
}
