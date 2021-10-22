library(dplyr)
library(ggplot2)
result <- NULL
unrun <- NULL
lam_seq <- c(0, 0.01, 0.02, 0.04, 0.09, 0.18, 0.38, 0.78, 1.62, 3.36, 6.95, 14.38, 29.76, 61.58, 127.43, 263.67, 545.56, 1128.84, 2335.72, 4832.93, 10000)
comp_seq <- c("none", "u", "v")
safe_read <- purrr::safely(readRDS)
for(lam in lam_seq){
  for(comp in comp_seq){
    file <- safe_read(paste0("G:\\My Drive\\Dissertation\\Spatial\\spatialr.dev\\development\\exp_2021_10_10\\spatial_",lam,"_",comp,".rds"))
    if(!is.null(file$result)){
      file$result$comp <- comp
      result <- c(list(file$result), result)
    } else {
      unrun <- c(unrun, paste(comp,lam))
    }
  }
}
coords <- readRDS("development\\exp_2021_10_10\\spatial_coords_sub.rds")
saveRDS(result,"development\\exp_2021_10_10\\result.rds")

p1 <- 1:length(result) %>% 
  purrr::map_dfr(.f = function(.x){
    mem <- result %>% purrr::pluck(.x)
    tibble(id = .x, normalization = mem$comp, lam = mem[[1]]$lambda, lik = mem[[1]]$lik)
  }) %>% 
  arrange(normalization, lik) %>% #View()
  ggplot(aes(x = lam, y = lik, color = normalization)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = log_trans()) +
  scale_x_log10(labels = trans_format('log', function(x) exp(x))) +
  xlab("lambda") +
  ylab("Likelihood") +
  ggtitle("Likelihood with varying lambda")

p2 <- 1:length(result) %>% 
  purrr::map_dfr(.f = function(.x){
    mem <- result %>% purrr::pluck(.x)
    tibble(id = .x, normalization = mem$comp, lam = mem[[1]]$lambda, lik = mem[[1]]$lik)
  }) %>% 
  dplyr::filter(normalization != "v") %>% 
  arrange(normalization, lik) %>% #View()
  ggplot(aes(x = lam, y = lik, color = normalization)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = log_trans()) +
  scale_x_log10(labels = trans_format('log', function(x) exp(x))) +
  scale_y_continuous(labels = scales::scientific) +
  scale_color_manual(values=c("#F8766D", "#00BA38")) +
  xlab("lambda") +
  ylab("Likelihood") +
  ggtitle("Likelihood with varying lambda (zoomed)")

grid.arrange(p1, p2)

plotly::ggplotly(p)

W <- readRDS("development/exp_2021_10_09/spatial_w.rds")
NMF::aheatmap(result[[21]][[1]]$u[which(W[1,]==1),], Colv = FALSE, Rowv = FALSE)
NMF::aheatmap(result[[21]][[1]]$u[sample(1:2000, 25),], Colv = FALSE, Rowv = FALSE)
NMF::aheatmap(result[[9]][[1]]$u, Colv = FALSE, Rowv = FALSE)
NMF::aheatmap(cor(result[[21]][[1]]$u, method = "spearman"), Colv = FALSE, Rowv = FALSE)

id1 <- 1:length(result) %>% 
  purrr::map_dfr(.f = function(.x){
    mem <- result %>% purrr::pluck(.x)
    tibble(id = .x, comp = mem$comp, lam = mem[[1]]$lambda, lik = mem[[1]]$lik)
  }) %>% 
  dplyr::filter(comp == "v") %>% 
  arrange(lam) %>% 
  pull(id)

lamb <- 1:length(result) %>% 
  purrr::map_dfr(.f = function(.x){
    mem <- result %>% purrr::pluck(.x)
    tibble(id = .x, comp = mem$comp, lam = mem[[1]]$lambda, lik = mem[[1]]$lik)
  }) %>% 
  dplyr::filter(comp == "v") %>% 
  arrange(lam) %>% 
  pull(lam)

pdf(file = "development/none_bc.pdf", width = 30, height = 30) 
par(mfrow = c(7,3))
mem <- 1
for(i in id1){
  NMF::aheatmap(cor(result[[i]][[1]]$u, method = "spearman"), Colv = FALSE, Rowv = FALSE, legend = FALSE, labRow = NA, labCol = NA, fontsize = 35, main = paste("lambda =",lamb[mem]))
  mem <- mem + 1
}
par(mfrow = c(7,3))
mem <- 1
for(i in id1){
  NMF::aheatmap(cor(result[[i]][[1]]$v, method = "spearman"), Colv = FALSE, Rowv = FALSE, legend = FALSE, labRow = NA, labCol = NA, fontsize = 35, main = paste("lambda =",lamb[mem]))
  mem <- mem + 1
}
dev.off()
dev.off()

pdf(file = "development/v_comp.pdf", width = 30, height = 30) 

pdf(file = "development/none_bc.pdf", width = 30, height = 30)
par(mfrow = c(7,3))
for(i in id1){
  clust <- kmeans(result[[i]][[1]]$u, centers = 5)
  data <- tibble(x = coords[,1], y = coords[,2], label = as.character(clust$cluster))
  
  p <- data %>% 
    ggplot(aes(x = x, y = y, color = label)) +
    geom_point()
  print(p)
}
dev.off()


p_list <- NULL
for(i in id1){
  clust <- kmeans(result[[i]][[1]]$u, centers = 5)
  data <- tibble(x = coords[,1], y = coords[,2], label = as.character(clust$cluster))
  
  p <- data %>% 
    ggplot(aes(x = x, y = y, color = label)) +
    geom_point() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  p_list <- c(p_list,list(p))
  
}
pdf(file = "development/v_comp_bc.pdf", width = 30, height = 30)
do.call("grid.arrange", c(p_list, ncol=3))
dev.off()
x <- readRDS("development/exp_2021_10_10/spatial_x_sub.rds")
sp_svd <- sparsesvd::sparsesvd(log(x+1),rank=20)
sp_pc <- irlba::irlba(x, nu=20)
cluster_svd <- kmeans(sp_svd$u%*%diag(sp_svd$d)%*%t(sp_svd$v), centers=5)
cluster_pc <- kmeans(sp_pc$u%*%diag(sp_pc$d)%*%t(sp_pc$v), centers=5)
data <- tibble(x = coords[,1], y = coords[,2], label = as.character(cluster_svd$cluster))

data %>% 
  ggplot(aes(x = x, y = y, color = label)) +
  geom_point() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
