library(dplyr)
library(spatialr.dev)
input=input
lambda=lambda
k=k

sample <- BayesSpace::getRDS("2020_maynard_prefrontal-cortex", input)
x <- assay(sample)

sp_svd <- sparsesvd::sparsesvd(log(x+1),rank=k)
u <- sp_svd$u %*% sqrt(diag(sp_svd$d))
v <- sp_svd$v %*% sqrt(diag(sp_svd$d))

coords <- SummarizedExperiment::colData(sample)[, c("row", "col")]
nn <- FNN::get.knn(coords, k = 6)
index <- nn$nn.index
neighbours <- as.vector(t(index))
w <- Matrix::sparseMatrix(i = rep(1:nrow(coords), each = 6), j = neighbours, x = 1/as.vector(t(nn$nn.dist)), dims = c(nrow(coords), nrow(coords)))
w <- as(w, "dgCMatrix")

start_time <- Sys.time()
model <- spatial_clust(input = sample, u_init = u, v_init = v, w = w, index = index, lambda = lambda, epsilon = 5e-3, maxiter = 1e3)
comp_time <- Sys.time()

time_diff <- difftime(comp_time, start_time, units = "secs")

out_string <- paste0("/N/project/zhangclab/alex/spatial/cv_2022_05_25/output/cv_",input,"_lambda_", lambda, "_k_",k,".rds")
result <- tibble(sample = input, time = as.numeric(time_diff), lambda = lambda, k = k, penal_lik = model$lik_penal)
saveRDS(result, out_string)