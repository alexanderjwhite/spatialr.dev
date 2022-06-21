library(bluster)
library(mclust)
library(dplyr)
library(BayesSpace)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(gridExtra)

input <- "151676"
folder <- "cv_2022_06_18_11"
folder_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/",folder,"/output/")
file <- list.files(folder_path)
file <- file[stringr::str_detect(file, ".rds")]
file_path <- paste0(folder_path,file)
out_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/2022_06_21_median_selection/")


# ehub <- ExperimentHub::ExperimentHub()

safe_ari <- purrr::safely(mclust::adjustedRandIndex)
clust_fun <- c("walktrap", "louvain", "infomap", "fast_greedy", "leading_eigen", "kmeans", "bayesspace_t", "bayesspace_n", "true")

result <- readRDS(file_path)
model <- result$model
data <- model$v


true_data <- getRDS("2020_maynard_prefrontal-cortex", input)
# true_data <- spatialCluster(true_data, q=7, d=2, platform='Visium', nrep=5000, gamma=3, save.chain=TRUE)
true <- as.numeric(true_data$layer_guess_reordered)

# true <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data/truth/",sample,".rds"))
# true_data <- BayesSpace::getRDS("2020_maynard_prefrontal-cortex", sample)

if(input %in% c("151669", "151670", "151671", "151672")){
  centers <- 5
} else {
  centers <- 7
}

results_raw <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_14_02/summary.csv")
colnames(results_raw) <- c("sample","time","lambda","k", "penal_lik", "cv_like")
results <- results_raw %>% 
  as_tibble() %>% 
  mutate(sample = as.numeric(sample)) %>% 
  arrange(sample, lambda, k) %>% 
  group_by(sample, lambda, k, penal_lik) %>% 
  summarize(cv_like_sum = sum(cv_like),
            cv_like_mean = mean(cv_like),
            cv_like_median = median(cv_like)) %>% 
  arrange(sample, k, lambda)

# Compile selection plots

p1 <- results %>% 
  filter(sample == as.numeric(input)) %>%
  ggplot(aes(x = lambda, y = cv_like_median)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = scales::log_trans()) +
  facet_wrap(~k) +
  ylab("Median of CV Likelihood")

p2 <- results %>% 
  filter(sample == input) %>% 
  group_by(k) %>% 
  slice_max(order_by = cv_like_median, n = 1) %>% 
  ungroup() %>% 
  mutate(selected = ifelse(penal_lik == max(penal_lik), "Yes", "No")) %>% 
  ggplot(aes(x = k, y = penal_lik, label = lambda)) +
  geom_point() +
  geom_line() +
  geom_label(aes(color = selected)) +
  ylab("Penalized Likelihood")


# Compile cluster plots

plots <<- NULL
clust_results <- clust_fun %>% 
  purrr::map_dfr(.f = function(.x){
    print(.x)
    
    
    
    if(.x == "kmeans"){
      clust <- kmeans(data,centers)$cluster
    } else if(.x == "true"){
      clust <- true
    } else if(.x == "bayesspace_t"){
      clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_t.rds"))
    } else if(.x == "bayesspace_n"){
      clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_norm.rds"))
    } else {
      clust <- as.numeric(clusterRows(data, BLUSPARAM=NNGraphParam(shared = FALSE, cluster.fun = .x)))
      
    }
    ari <- safe_ari(clust,clue::solve_LSAP(table(true, clust), maximum = TRUE)[true])
    
    if(!is.null(ari$result)){
      p <- clusterPlot(true_data, label=clust, palette=NULL, size=0.05) +
        scale_fill_viridis_d(option = "A") +
        labs(title=paste(.x,"ari =",round(ari$result, 4))) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
      plots <<- c(plots,list(p))
    } 
    
    return(tibble(cluster = .x, ari = ari$result))
  })


pdf(file = paste0(out_path,input,".pdf"), width = 10, height = 10)
do.call("grid.arrange", c(list(p1), ncol=1))
do.call("grid.arrange", c(list(p2), ncol=1))
do.call("grid.arrange", c(plots, ncol=3))
dev.off()

