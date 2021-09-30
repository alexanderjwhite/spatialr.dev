paths <- .libPaths()
paths <- c("/geode2/home/u100/whitealj/BigRed3/r_4_0_4_library/",paths)
.libPaths(paths)
load("/geode2/home/u100/whitealj/BigRed3/scripts/spatialLIBD.spe.RData")
library(spatialr.dev)

X <- expr_counts[["counts"]]
X <- t(X)

sp_svd <-  sparsesvd::sparsesvd(log(X+1),rank=20)

U0 <- sp_svd$u
V0 <- sp_svd$v %*% diag(sp_svd$d)

W <- fct_dist_matrix(coords, distance = "euclidean", method = "knn_1", k = NULL,alpha = 1, verbose = TRUE)

test <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 15, norm_comp = "u", eta = 100, coords=NULL ,lambda = NULL, grid = 5, w=W, optimizer = fct_opt_amsgrad, fast = TRUE, cores = 16)
saveRDS(test, "/geode2/home/u100/whitealj/BigRed3/scripts/test.rds")
# no_norm <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 5e3, norm_comp = "none", eta = 100, coords=NULL ,lambda = NULL, grid = 10, w=W, optimizer = fct_opt_amsgrad, fast = TRUE, cores = 16)
# u_norm <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 5e3, norm_comp = "u", eta = 100, coords=NULL ,lambda = NULL, grid = 10, w=W, optimizer = fct_opt_amsgrad, fast = TRUE, cores = 16)
# v_norm <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 5e3, norm_comp = "v", eta = 100, coords=NULL ,lambda = NULL, grid = 10, w=W, optimizer = fct_opt_amsgrad, fast = TRUE, cores = 16)
#  
# saveRDS(no_norm, "/geode2/home/u100/whitealj/BigRed3/scripts/no_norm.rds")
# saveRDS(u_norm, "/geode2/home/u100/whitealj/BigRed3/scripts/u_norm.rds")
# saveRDS(v_norm, "/geode2/home/u100/whitealj/BigRed3/scripts/v_norm.rds")