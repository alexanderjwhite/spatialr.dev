library(dplyr)
library(BayesSpace)
path <- "/N/project/zhangclab/alex/spatial/cv_2022_06_01/output/"
files <- list.files(path)
files <- files[stringr::str_detect(files, "_ari.rds$")]
files <- sample(files, length(files), replace = FALSE)
files %>% 
  purrr::walk(.f = function(.x){
    res <- readRDS(paste0(path,.x)) %>% 
      select(-model)
    
    write.table(res,
              file = "/N/project/zhangclab/alex/spatial/cv_2022_06_01/summary.csv", 
              append = TRUE,
              row.names = FALSE,
              col.names = FALSE)
    
  })
# saveRDS(results,paste0("/N/project/zhangclab/alex/spatial/cv_2022_06_01/summary.rds"))

path_pre <- "/N/project/zhangclab/alex/spatial/cv_2022_06_03/output/clust_"
results <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_03/summary.csv") %>% unique()
colnames(results) <- c("sample", "centers", "K", "norm", "lambda", "clust", "seed", "ari_shared", "ari_k", "ari_km")
clean <- results %>% 
  tibble() %>% 
  mutate(sample = as.numeric(sample),
         centers = as.numeric(centers),
         K = as.numeric(K),
         lambda = as.numeric(lambda),
         seed = as.numeric(seed),
         ari_k = as.numeric(ari_k),
         ari_shared = as.numeric(ari_shared),
         ari_km = as.numeric(ari_km)) %>% 
  mutate(ari_max := pmax(ari_k, ari_shared, ari_km, na.rm = TRUE)) %>% 
  mutate(path = paste0(path_pre,sample,"_",K,"_",lambda,"_",ifelse(is.null(norm),"NULL",norm),"_",seed,"_",clust,"_k")) %>% 
  filter(clust != "label_prop") %>%
  group_by(sample) %>% 
  slice_max(order_by = ari_max, n = 2, with_ties = FALSE) 

# summarized <- clean %>% 
#    
#   mutate(path = paste0(path_pre,sample,"_",K,"_",lambda,"_",ifelse(is.null(norm),"NULL",norm),"_",seed,"_",clust,"_k")) %>% 
#   group_by(sample, centers, K, norm, lambda, clust) %>% 
#   mutate(ari_shared = max(ari_shared, na.rm = TRUE),
#          ari_k = max(ari_k, na.rm = TRUE),
#          ari_km = max(ari_km, na.rm = TRUE),
#          ari_max = max(ari_max, na.rm = TRUE)) %>% 
#   ungroup() %>% 
#   group_by(sample, clust) %>% 
#   slice_max(order_by = ari_max, n = 3, with_ties = FALSE) 

sample <- BayesSpace::getRDS("2020_maynard_prefrontal-cortex", "151507")
# labels <- dplyr::recode(dlpfc$spatial.cluster, 3, 4, 5, 6, 2, 7, 1)

## View results
clusterPlot(sample, label=clust_res, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "A", labels = 1:7) +
  labs(title="BayesSpace")
