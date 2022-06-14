library(dplyr)
# library(spatialr.dev)
library(Matrix)
input="151507"
lambda=1
k=20
normc=500
cv=100
seed=1
set.seed(seed)

sample <- BayesSpace::getRDS("2020_maynard_prefrontal-cortex", input)
if(is.null(normc)){norm_string <- "NULL"} else {norm_string <- normc}
in_string_mat <- paste0("/N/project/zhangclab/ShaCao/single_cell/spatial_matlab/dlpfc",input,"/dlpfc",input,".mat")
in_string_coords <- paste0("/N/project/zhangclab/ShaCao/single_cell/spatial_matlab/dlpfc",input,"/dlpfc",input,"coords.mat")
x <- as(R.matlab::readMat(in_string_mat)$x, "dgCMatrix")
coords <- as(R.matlab::readMat(in_string_coords)$coords, "dgCMatrix")
x <- as(R.matlab::readMat("G:/My Drive/Dissertation/Spatial/spatial_analysis/dlpfc151507.mat")$x, "dgCMatrix")
coords <- as(R.matlab::readMat("G:/My Drive/Dissertation/Spatial/spatial_analysis/dlpfc151507coords.mat")$coords, "dgCMatrix")

full_nn <- FNN::get.knn(coords, k = 6)
full_index <- full_nn$nn.index
full_dist <- full_nn$nn.dist

# Find non-neighboring 
test_cols <- NULL
neighbor_map <- NULL
n <- 0
iter <- 0
while ((n <= cv) & iter < nrow(coords)){
  sample_i <- sample(1:nrow(coords), 1, replace = FALSE)
  if (!(sample_i %in% neighbor_map)){
    test_cols <- c(test_cols,sample_i)
    neighbor_map <- c(neighbor_map,sample_i,full_index[sample_i,])
    n <- n + 1
  }
  iter <- iter + 1
}