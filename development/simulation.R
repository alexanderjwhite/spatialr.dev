library(dials)
library(dplyr)
library(ggplot2)
library(gridExtra)

k <- 5
size <- c(2000,4000)
coef <- seq(5,8, length.out = k)
sigma <- 1
u_size <- c(size[1],k)
v_size <- c(size[2],k)
param_x <- new_quant_param(type = "double", range = c(0,1e3), inclusive = c(TRUE,TRUE),  label = c(x = "x"))
param_y <- new_quant_param(type = "double", range = c(0,1e3), inclusive = c(TRUE,TRUE),  label = c(y = "y"))

x <- dials::grid_latin_hypercube(param_x, size = size[1]) %>% 
  round(digits = 2) %>% 
  pull()
y <- dials::grid_latin_hypercube(param_y, size = size[1]) %>% 
  round(digits = 2) %>% 
  pull()



km <- kmeans(cbind(x,y), centers = k)

data <- tibble(x = x, y = y, id = km$cluster) %>% 
  arrange(id)

clust_assign_true <- data %>% pull(id)

coords <- data %>% 
  select(x,y) %>% 
  as.data.frame()

true_plot <- data %>% 
  ggplot() + 
  geom_point(aes(x = x, y = y, color = as.character(id))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())


u <- matrix(0,nrow = u_size[1], ncol = u_size[2])
v <- matrix(0,nrow = v_size[1], ncol = v_size[2])
r_start <- 1
c_start <- 1

for (j in 1:k) {
  n_rows <- data %>% filter(id == j) %>% nrow()
  r_end <- r_start + n_rows - 1
  c_end <- c_start + n_rows - 1
  cl_U <- qr.Q(qr(replicate(1, rnorm(n_rows, 0, 1))))
  cl_V <- t(qr.Q(qr(replicate(1, rnorm(n_rows, 0, 1)))))
  u[r_start:r_end,j] <- coef[j]*cl_U
  v[c_start:c_end,j] <- coef[j]*cl_V
  r_start <- r_end+1
  c_start <- c_end+1
}
signal <- u%*%t(v)
# noise <- matrix(rnorm(size[1]*size[2], mean = 0, sd = sigma), nrow = size[1], ncol = size[2])
uve <- signal #+ noise

param_mat <- exp(uve)
x <- rpois(length(param_mat), param_mat)
dim(x) <- dim(param_mat)
x <- as(x, "dgCMatrix")

NMF::aheatmap(u - min(u), Colv = NA, Rowv = NA)
NMF::aheatmap(v, Colv = NA, Rowv = NA)
NMF::aheatmap(signal, Colv = NA, Rowv = NA)
NMF::aheatmap(uve, Colv = NA, Rowv = NA)
NMF::aheatmap(x, Colv = NA, Rowv = NA)

W <- fct_dist_matrix(coords, distance = "euclidean", method = "knn_1", k = NULL,alpha = 1, verbose = TRUE)

sp_svd <- sparsesvd::sparsesvd(x,rank=20)

U0 <- sp_svd$u
V0 <- sp_svd$v %*% diag(sp_svd$d)

model <- spatial_clust(x=x, u_init=U0, v_init=V0, max_iter = 1e3, norm_comp = "u", eta = 100, coords=NULL ,lambda = 1e4, grid = 5, w=W, optimizer = fct_opt_amsgrad, fast = FALSE, 1)

clust <- kmeans(u - min(u), centers = 5)
clust_assign <- as.character(clust$cluster)
# viz <- tibble(x = coords[,1], y = coords[,2], label = clust_assign)
# viz %>% 
#   ggplot(aes(x = x, y = y, color = label)) +
#   geom_point() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())




shuffle <- clue::solve_LSAP(table(clust_assign_true, as.numeric(clust_assign)), maximum = TRUE)
tbl <- (table(clust_assign_true, as.numeric(clust_assign))[,shuffle])
tbl
sum(diag(tbl))/size[1]

iter <- 1
swap <- rep(0,size[1])
for (i in shuffle){
  val <- which(clust_assign == i)
  swap[val] <- iter
  iter <- iter + 1
}

viz <- tibble(x = coords[,1], y = coords[,2], label = swap, true = clust_assign_true) %>% 
  mutate(correct = ifelse(true==label,"True", "False"))
assess_plot <- viz %>% 
  ggplot(aes(x = x, y = y, color = as.character(label))) +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete()

correct_plot <- viz %>% 
  ggplot(aes(x = x, y = y, color = as.character(correct))) +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_discrete()


grid.arrange(true_plot, assess_plot)
