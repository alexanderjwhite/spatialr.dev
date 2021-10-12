#' Internal Distance Matrix Computation
#'
#' @inheritParams spatial_clust
#'
#' @return distance matrix
#' @export
#' 
#' @import stats
#' @import FNN
#' @import Matrix
#'
#' @examples fct_dist_matrix(matrix(rnorm(20), nrow = 10, ncol = 2))
fct_dist_matrix <- function(coords, distance = "euclidean", method = "dist", k = NULL, alpha = 1, verbose = TRUE){
  # ToDo: conditions for coord dataframe
  if(verbose){print("Computing Distance Matrix...")}
  if(method == "dist"){
    dist_matrix <- as.matrix(stats::dist(coords, method = distance))
    w <- as.matrix(exp(-alpha*dist_matrix))
  } else {
    if(is.null(k)){k <- sqrt(nrow(coords))/2}
    if(verbose){print(paste("Using KNN method with k =",as.integer(k),"..."))}
    nn <- FNN::get.knn(coords, k = k)
    neighbours <- as.vector(t(nn$nn.index))
    #neighbours <- as.integer(apply(dist_matrix, 1, function(x) sort(x, index.return = TRUE)$ix[2:(k + 1)]))
    #adj_matrix <- (Matrix::sparseMatrix(i = rep(1:nrow(coords), each = k), j = neighbours, x = 1, dims = c(nrow(coords), nrow(coords))))
    if(method == "knn_1"){
      w <- Matrix::sparseMatrix(i = rep(1:nrow(coords), each = k), j = neighbours, x = 1, dims = c(nrow(coords), nrow(coords)))
      w <- as(w, "dgCMatrix")
    } else if(method == "knn_2"){
      w <- Matrix::sparseMatrix(i = rep(1:nrow(coords), each = k), j = neighbours, x = as.vector(t(nn$nn.dist)), dims = c(nrow(coords), nrow(coords)))
      w <- as(w, "dgCMatrix")
    }
  }
  if(verbose){print("Done")}
  return(w)
}
