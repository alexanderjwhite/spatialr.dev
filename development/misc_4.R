load("development/spatialLIBD.spe.RData")
library(dplyr)
library(ggplot2)
#library(spatialr.dev)
X <- expr_counts[["counts"]][1:1000,1:200]
X <- t(X)

dim(X)
N <- colSums(X)

dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}
dispersion_test(N)



# 
library(SingleCellExperiment)
# sample <- BayesSpace::getRDS(dataset="2018_thrane_melanoma", sample="ST_mel1_rep2")
sample <- BayesSpace::getRDS("2020_maynard_prefrontal-cortex", "151507")
x <- assay(sample)

sp_svd <- sparsesvd::sparsesvd(log(x+1),rank=20)
u <- sp_svd$u
v <- sp_svd$v %*% diag(sp_svd$d)

coords <- colData(sample)[, c("row", "col")]
nn <- FNN::get.knn(coords, k = 10)
index <- nn$nn.index
neighbours <- as.vector(t(index))
w <- Matrix::sparseMatrix(i = rep(1:nrow(coords), each = 10), j = neighbours, x = 1, dims = c(nrow(coords), nrow(coords)))
w <- as(w, "dgCMatrix")

time <- Sys.time()
model <- spatial_clust(input = sample, u_init = u, v_init = v, w = w, lambda = 1, epsilon = 5e-3, maxiter = 1e3)
Sys.time()-time
