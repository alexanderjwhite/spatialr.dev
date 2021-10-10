load("development/spatialLIBD.spe.RData")
library(dplyr)
library(ggplot2)
#library(spatialr.dev)
X <- expr_counts[["counts"]][1:100,1:200]
X <- t(X)

sp_svd <- sparsesvd::sparsesvd(log(X+1),rank=20)

U0 <- sp_svd$u
V0 <- sp_svd$v %*% diag(sp_svd$d)

test <- FNN::get.knn(coords[1:200,])
test_neighbours <- as.vector(t(test$nn.index))
Matrix::sparseMatrix(i = rep(1:nrow(coords[1:200,]), each = k), j = test_neighbours, x = 1, dims = c(200, 200))

dist_matrix <- as.matrix(stats::dist(coords[1:200,], method = "euclidean"))
neighbours <- as.integer(apply(dist_matrix, 1, function(x) sort(x, index.return = TRUE)$ix[2:(k + 1)]))
cur <- matrix(neighbours, ncol = 10, byrow = TRUE)
mat <- as.matrix(Matrix::sparseMatrix(i = rep(1:nrow(coords[1:200,]), each = k), j = neighbours, x = 1, dims = c(nrow(coords[1:200,]), nrow(coords[1:200,]))))

W=fct_dist_matrix(coords[1:200,], distance = "euclidean", method = "knn_1", k = NULL,alpha = 1, verbose = TRUE)

large_penal <- spatial_clust(x=(X), u_init=U0, v_init=V0, max_iter = 5, norm_comp = "u", eta = 100, coords=NULL ,lambda = 0, grid = 5, w=W, optimizer = fct_opt_amsgrad, fast = FALSE)


diag(apply(V0, 2, function(y){norm(y, type="2")}))

NMF::aheatmap(test2$u[which(W[1,]==1),], Colv = FALSE, Rowv = FALSE)
NMF::aheatmap(U0[which(W[1,]==1),], Colv = FALSE, Rowv = FALSE)


NMF::aheatmap(large_penal$u, Colv = FALSE, Rowv = FALSE)



test2 <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 100, norm_comp = "u", eta = 100, coords=NULL ,lambda = 1e3, w=W, optimizer = fct_opt_amsgrad, fast = FALSE)

purrr::map_dfr(large_penal, .f = function(.x){
  tibble(lik = .x$lik, penal = .x$penal, lam = .x$lambda)
}) %>% 
  mutate(lam_penal = lam*penal) %>% 
  ggplot() +
  geom_point(aes(x = log(lam), y = lik))


large_penal[[3]]

NMF::aheatmap(large_penal[[5]]$u[which(W[2,]==1),], Colv = FALSE, Rowv = FALSE)

test5 <- spatial_clust(x=(as.matrix(X)), u_init=U0, v_init=V0, max_iter = 1000, coords=coords1,lambda = 1e7, w=W, optimizer = fct_opt_amsgrad, fast = FALSE)


# Output the objective function

test4[[4]]$ which(W[1,]==1)




NMF::aheatmap(test5$u[which(W[4,]==1),], Colv = FALSE, Rowv = FALSE)
