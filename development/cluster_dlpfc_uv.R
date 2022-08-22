library(bluster)
library(mclust)
library(dplyr)
library(BayesSpace)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(gridExtra)

safe_ari <- purrr::safely(mclust::adjustedRandIndex)
results_raw_1 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_27_01/summary.csv")
results_raw_2 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_07_02_01/summary.csv")
results_raw_3 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_07_06_01/summary.csv")
results_raw_4 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_07_06_02/summary.csv")
results_raw <- rbind(rbind(rbind(results_raw_1, results_raw_2),  results_raw_3), results_raw_4)
colnames(results_raw) <- c("sample","time","lambda","k", "penal_lik", "cv_like")
lambda_filt <- c(5,10,15,20,25)
k_filt <- c(10,20,30,40)
results <- results_raw %>% 
  tibble() %>% 
  filter(lambda %in% lambda_filt & k %in% k_filt) %>% 
  mutate(sample = as.numeric(sample)) %>% 
  # group_by(sample, lambda, k) %>% 
  # dplyr::slice(n = 1) %>% 
  # arrange(sample, lambda, k) %>% 
  group_by(sample, time, lambda, k, penal_lik) %>% 
  summarize(cv_like_sum = sum(cv_like),
            cv_like_mean = mean(cv_like),
            cv_like_median = median(cv_like)) %>% 
  group_by(sample, lambda, k) %>% 
  dplyr::slice(n = 1) %>% 
  arrange(sample, k, lambda)


# results %>% filter(k==5,sample=="151673")


# Compile selection plots

input <- "151676"
if(input %in% c("151669", "151670", "151671", "151672")){
  centers <- 5
} else {
  centers <- 7
}
set.seed(1)

p_grid <- results %>%
  filter(sample == as.numeric(input)) %>%
  ggplot(aes(x = lambda, y = cv_like_mean)) +
  geom_point() +
  geom_line() +
  facet_wrap(~k, scales = "free") +
  ylab("Mean of CV Likelihood")

p_k <- results %>%
  filter(sample == input) %>%
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
in_string <- paste0("/N/project/zhangclab/alex/spatial/cv_2022_06_29_01/output/clust_",input, "_",k,"_", scales::number(lambda, accuracy = 0.1),"_",norm_string,"_",seed,".rds")

in_string

model_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_29_01/output/clust_",input, "_",k,"_", scales::number(lambda, accuracy = 0.1),"_",norm_string,"_",seed,".rds")
out_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/2022_07_08_selection/")
result <- readRDS(model_path)
model <- result$model
saveRDS(model,paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/smoothlrc/models/",input,".rds"))


data <- model$v

true_data <- getRDS("2020_maynard_prefrontal-cortex", input)
true <- as.numeric(true_data$layer_guess_reordered)


# Compile cluster plots

plots <- NULL
# clust <- kmeans(data,centers, nstart = 50, iter.max = 100)$cluster
# shuffled <- clue::solve_LSAP(table(clust, true), maximum = TRUE)[clust]
# ari <- safe_ari(true,shuffled)
# p0 <- clusterPlot(true_data, label=shuffled, palette=NULL, size=0.05) +
#   scale_fill_viridis_d(option = "A") +
#   labs(title=paste("kmeans ari =",round(ari$result, 4))) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")
# plots <- c(plots,list(p0))

clust <- Mclust(data, G = centers, modelNames = "EEE")$classification
shuffled <- clue::solve_LSAP(table(clust, true), maximum = TRUE)[clust]
saveRDS(shuffled, paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/smoothlrc/",input,".rds"))
ari <- safe_ari(true,shuffled)
ari$result


p1 <- clusterPlot(true_data, label=shuffled, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A") +
  labs(title=paste("mclust ari =",round(ari$result, 4))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
plots <- c(plots,list(p1))

clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_t.rds"))
shuffled <- clue::solve_LSAP(table(clust, true), maximum = TRUE)[clust]
ari <- safe_ari(true,shuffled)
p2 <- clusterPlot(true_data, label=shuffled, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A") +
  labs(title=paste("bayes t ari =",round(ari$result, 4))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
plots <- c(plots,list(p2))

clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_norm.rds"))
shuffled <- clue::solve_LSAP(table(clust, true), maximum = TRUE)[clust]
ari <- safe_ari(true,shuffled)
p3 <- clusterPlot(true_data, label=shuffled, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A") +
  labs(title=paste("bayes n ari =",round(ari$result, 4))) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
plots <- c(plots,list(p3))


p4 <- clusterPlot(true_data, label=true, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A") +
  labs(title=paste("true ari = 1")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")
plots <- c(plots,list(p4))


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
