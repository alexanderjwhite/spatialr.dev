load("spatialLIBD.spe.RData")
library(spatialr.dev)

X <- expr_counts[["counts"]]
X <- t(X)

sp_svd <- sparsesvd::sparsesvd(log(X+1),rank=20)

U0 <- sp_svd$u
V0 <- sp_svd$v %*% diag(sp_svd$d)

W <- fct_dist_matrix(coords, distance = "euclidean", method = "knn_1", k = NULL,alpha = 1, verbose = TRUE)

saveRDS(X, "spatial_x.rds")
saveRDS(coords, "spatial_coords.rds")
saveRDS(U0,"spatial_u.rds")
saveRDS(V0,"spatial_v.rds")
saveRDS(W,"spatial_w.rds")

