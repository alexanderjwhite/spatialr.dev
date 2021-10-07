paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/BigRed3/r_4_0_4_library/",paths)
.libPaths(paths)
library(spatialr.dev)
X <- readRDS("/geode2/home/u100/whitealj/Carbonate/scripts/spatial_x.rds")
coords <- readRDS("/geode2/home/u100/whitealj/Carbonate/scripts/spatial_coords.rds")

sp_svd <- sparsesvd::sparsesvd(log(X+1),rank=20)

U0 <- sp_svd$u
V0 <- sp_svd$v %*% diag(sp_svd$d)

W <- fct_dist_matrix(coords, distance = "euclidean", method = "knn_1", k = NULL,alpha = 1, verbose = TRUE)

saveRDS(U0,"/geode2/home/u100/whitealj/Carbonate/scripts/spatial_u.rds")
saveRDS(V0,"/geode2/home/u100/whitealj/Carbonate/scripts/spatial_v.rds")
saveRDS(W,"/geode2/home/u100/whitealj/Carbonate/scripts/spatial_w.rds")

