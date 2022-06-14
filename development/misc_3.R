library(ggplot2)
files <- list.files("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_05_25/output")
results <- files %>% 
  purrr::map_dfr(.f = function(.x){
    readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_05_25/output/",.x))
  })
results %>% 
  ggplot(aes(x = k, y = penal_lik)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample)


# lambda selection

files <- list.files("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_05_30/output")
results <- files %>% 
  purrr::map_dfr(.f = function(.x){
    readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_05_30/output/",.x))
  })

p2 <- results %>% 
  group_by(sample, lambda, k) %>% 
  summarize(like = mean(cv_like_val), time = mean(time)) %>% 
  arrange(sample, k, lambda) %>% 
  ggplot(aes(x = lambda, y = like)) +
  geom_point() +
  geom_line() +
  facet_wrap(~sample)
plotly::ggplotly(p2)
