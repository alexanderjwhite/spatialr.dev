# results_raw1 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_23_04/summary.csv", fill = TRUE)
# results_raw2 <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_24_01/summary.csv")
# results_raw <- rbind(results_raw1, results_raw2)
# colnames(results_raw) <- c("sample","time","lambda","k", "ari")
# 
# results_best <- results_raw %>% 
#   tibble() %>% 
#   group_by(sample) %>% 
#   slice_max(order_by = ari, n = 1, with_ties = FALSE)


# mclust
library(dplyr)
library(mclust)
library(BayesSpace)
library(ggplot2)
library(gridExtra)
safe_ari <- purrr::safely(mclust::adjustedRandIndex)
clust_fun <- c("mclust", "bayesspace_t", "bayesspace_n", "true")



results_mclust <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_26_01/summary.csv", fill = TRUE)
colnames(results_mclust) <- c("sample","lambda","k", "ari_mean", "ari_sd")
results_mclust_best <- results_mclust %>% 
  tibble() %>% 
  group_by(sample) %>% 
  slice_max(order_by = ari_mean, n = 1, with_ties = FALSE)

results %>% 
  tibble() %>% 
  distinct() %>% 
  group_by(sample, k) %>% 
  summarize(n = n()) %>% 
  View()

results_mclust %>% 
  tibble() %>% 
  select(sample, lambda,k) %>% 
  distinct() %>% 
  filter(sample == 151509) %>% 
  arrange(k,lambda)
  

set.seed(1)
out_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/2022_06_26_selection/")
results_mclust_best %>% 
  pull(sample) %>% 
  # "151510" %>% 
  purrr::walk(.f = function(.x){
    print(.x)
    input <- .x
    if(input %in% c("151669", "151670", "151671", "151672")){
      centers <- 5
    } else {
      centers <- 7
    }
    true_data <- getRDS("2020_maynard_prefrontal-cortex", input)
    true <- as.numeric(true_data$layer_guess_reordered)
    lambda <- results_mclust_best %>% filter(sample == input) %>% pull(lambda)
    k <- results_mclust_best %>% filter(sample == input) %>% pull(k)
    norm_string <- "NULL"
    seed <- 1
    
    model_path <- paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_23_04/output/clust_",input, "_",k,"_", lambda,"_",norm_string,"_",seed,".rds")
    result <- readRDS(model_path)
    model <- result$model
    data <- model$v
    
    plots <- NULL
    clust <- Mclust(data, G = centers)$classification
    shuffled <- clue::solve_LSAP(table(clust, true), maximum = TRUE)[clust]
    ari <- safe_ari(true,shuffled)
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
    do.call("grid.arrange", c(plots, ncol=2))
    dev.off()
    
  })
