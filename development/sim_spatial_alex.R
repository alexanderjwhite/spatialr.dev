lam_param=lam_param
comp_param=comp_param
library(spatialr.dev)

X <- readRDS("spatial_x_sub.rds")
coords <- readRDS("spatial_coords_sub.rds")
W <- readRDS("spatial_w.rds")
U0 <- readRDS("spatial_u.rds")
V0 <- readRDS("spatial_v.rds")

result <- spatial_clust(x=as.matrix(X), u_init=U0, v_init=V0, max_iter = 5e3, norm_comp = comp_param, eta = 100, verbose = TRUE, coords=NULL ,lambda = lam_param, grid = 5, w=W, optimizer = fct_opt_amsgrad, fast = TRUE, cores = 4)
fn <- paste0("spatial_",lam_param,"_",comp_param,".rds")
saveRDS(result, fn)
