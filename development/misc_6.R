results_raw <- read.table("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_09/summary.csv") 
colnames(results_raw) <- c("sample", "centers", "k","norm_string","lambda","clust_func","rand_seed","ari1","ari2","ari3")

results_raw %>% 
  tibble() %>% 
  group_by(sample) %>% 
  slice_max(n = 1, order_by = ari2, with_ties = FALSE)
  # mutate(path = out_string <- paste0("/N/project/zhangclab/alex/spatial/cv_2022_06_08/output/clust_",sample, "_",k,"_", lambda,"_",norm_string,"_",rand_seed,"_",clust_func)) %>% 
  # pull(path)

res <- results_raw[,c(1,5)] %>% unique() %>% tibble() %>% rename("sample" = V1, "lambda" = V5)
res %>% 
  arrange(sample, lambda) %>% 
  View()
