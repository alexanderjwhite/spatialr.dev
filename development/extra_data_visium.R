library(bluster)
library(mclust)
library(dplyr)
library(BayesSpace)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(Seurat)

results_raw <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_08_16_01/summary.csv")

colnames(results_raw) <- c("sample", "preproc","time","lambda","k", "penal_lik", "cv_like")
# lambda_filt <- c(5,10,15,20,25)
# k_filt <- c(10,20,30,40)
results <- results_raw %>% 
  tibble() %>% 
  # filter(lambda %in% lambda_filt & k %in% k_filt) %>% 
  # mutate(sample = as.numeric(sample)) %>% 
  # group_by(sample, lambda, k) %>% 
  # dplyr::slice(n = 1) %>% 
  # arrange(sample, lambda, k) %>% 
  group_by(sample, preproc, time, lambda, k, penal_lik) %>% 
  summarize(cv_like_sum = sum(cv_like),
            cv_like_mean = mean(cv_like),
            cv_like_median = median(cv_like)) %>% 
  group_by(sample, preproc, lambda, k) %>% 
  dplyr::slice(n = 1) %>% 
  arrange(sample, preproc, k, lambda)

results %>% 
  group_by(sample, preproc, k) %>% 
  slice_max(order_by = cv_like_mean, n = 1) %>% 
  group_by(sample, preproc) %>% 
  slice_max(order_by = penal_lik, n = 1) %>% 
  View()

input <- "data_1"
pre <- "before"
load(paste0("G:/My Drive/Dissertation/Spatial/data/LuLabSpatialData/LuLabSpatial_seuObj_",input,".RData"))
sce <- as.SingleCellExperiment(visiumData)
# cdata <- visiumData@images$B1@coordinates
# cdata[,"spatial.cluster"] <- 
colData(sce)<- cbind(colData(sce),visiumData@images$A1@coordinates)
# sce <- visiumData
sce <- spatialPreprocess(sce, platform="Visium", skip.PCA=TRUE)
# bs_colnames <- c("col","row","spatial.cluster")

n_clust <- 2:8
p_grid <- results %>%
  filter(sample == input, preproc == pre) %>%
  ggplot(aes(x = lambda, y = cv_like_mean)) +
  geom_point() +
  geom_line() +
  facet_wrap(~k, scales = "free") +
  ylab("Mean of CV Likelihood")

p_k <- results %>%
  filter(sample == input, preproc == pre) %>%
  group_by(k) %>%
  slice_max(order_by = cv_like_mean, n = 1) %>%
  ungroup() %>%
  mutate(selected = ifelse(penal_lik == max(penal_lik), "Yes", "No")) %>%
  ggplot(aes(x = k, y = penal_lik, label = lambda)) +
  geom_point() +
  geom_line() +
  geom_label(aes(color = selected)) +
  ylab("Penalized Likelihood")

lambda <- results %>%
  filter(sample == input) %>%
  group_by(k) %>%
  slice_max(order_by = cv_like_mean, n = 1) %>% 
  ungroup() %>% 
  slice_max(order_by = penal_lik, n = 1) %>% 
  pull(lambda)

k <- results %>%
  filter(sample == input) %>%
  group_by(k) %>%
  slice_max(order_by = cv_like_mean, n = 1) %>% 
  ungroup() %>% 
  slice_max(order_by = penal_lik, n = 1) %>% 
  pull(k)
norm_string <- "NULL"
seed <- 1
set.seed(seed)
# in_string <- paste0("/N/project/zhangclab/alex/spatial/cv_2022_06_29_01/output/clust_",input, "_",k,"_", scales::number(lambda, accuracy = 0.1),"_",norm_string,"_",seed,".rds")
# 
# in_string

model_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_07_13_03/output/clust_",input, "_",k,"_", scales::number(lambda, accuracy = 0.1),"_",norm_string,"_",seed,".rds")
out_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/lulab_selection/")
result <- readRDS(model_path)
model <- result$model
data <- model$v


# Compile cluster plots

plots <<- NULL
# SpatialDimPlot(sce)
n_clust %>% 
  purrr::walk(.f = function(.x){
    clust <- Mclust(data, G = .x, modelNames = "EEE")$classification
    saveRDS(clust, paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/lulab_selection/cluster/",input,"/",input,"_",.x,".rds"))
    # p <- clusterPlot(sce, label=clust, palette=NULL, size=0.05, color = NA) +
    #   scale_fill_viridis_d(option = "B") +
    #   labs(title=paste("k =",.x)) +
    #   theme(axis.title.x = element_blank(),
    #         axis.title.y = element_blank())
    # plots <<- c(plots,list(p))
  })



pdf(file = paste0(out_path,input,".pdf"), width = 10, height = 10)
do.call("grid.arrange", c(list(p_grid), ncol=1))
do.call("grid.arrange", c(list(p_k), ncol=1))
do.call("grid.arrange", c(plots, ncol=2))
dev.off()

# 
# pdf(file = paste0(out_path,input,"_",k,"_",lambda,".pdf"), width = 10, height = 10)
# # do.call("grid.arrange", c(list(p1), ncol=1))
# # do.call("grid.arrange", c(list(p2), ncol=1))
# do.call("grid.arrange", c(plots, ncol=3))
# dev.off()
# 
# 
# set.seed(1)
# M <- 30
# clust_range <- 1:M %>% 
#   purrr::map_dfr(.f = function(.x){
#     print(.x)
#     clust <- mclust::Mclust(data, G = centers)$classification
#     # clust <- kmeans(data,centers, nstart = 50, iter.max = 100)$cluster
#     ari <- safe_ari(clust,clue::solve_LSAP(table(true, clust), maximum = TRUE)[true])
#     tibble(id = .x, ari = ari$result)
# })
# 
# clust_range %>% 
#   summarize(med_ari = median(ari),
#             mean_ari = mean(ari),
#             sd_ari = sd(ari),
#             max_ari = max(ari),
#             min_ari = min(ari))
# 
# 
