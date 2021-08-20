#' Internal Distance Matrix Computation
#'
#' @param coords matrix or data frame; euclidean coordinates in separate columns
#' @param verbose TRUE or FALSE; print to screen?
#'
#' @return distance matrix
#' @export
#'
#' @examples fct_dist_matrix(matrix(rnorm(20), nrow = 10, ncol = 2), verbose = TRUE)
fct_dist_matrix <- function(coords, verbose = TRUE){
  # ToDo: conditions for coord dataframe
  if(verbose){print("Computing Distance Matrix...")}
  size <- nrow(coords)
  w <- matrix(0, nrow = size, ncol = size)
  red_size <- (size^2-size)/2
  count <- 0
  for(j in 2:size){
    for(i in 1:(j-1)){
      count <- count+1
      if((count)%%(red_size/1e3)==0 & verbose){
        cat("\r","")
        cat("\r",paste0("Computing Distance Matrix: ", round((count/red_size),digits = 5)*100,"%"))
        }
      w[i,j] <- 1/norm(coords[i,]-coords[j,], type="2")
    }
  }
  w <- w + t(w)
  diag(w) <- 0
  if(verbose){cat("\n","Done \n")}
  return(w)
}
