set.seed(1)
input <- "151507"
path <- "G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_21_02/output/"
files <- list.files(path)
true_data <- getRDS("2020_maynard_prefrontal-cortex", input)
cv <- 10
coords <- SummarizedExperiment::colData(true_data)[, c("row", "col")]

full_nn <- FNN::get.knn(coords, k = 6)
full_index <- full_nn$nn.index
full_dist <- full_nn$nn.dist

# Find non-neighboring 
test_cols <- NULL
neighbor_map <- NULL
n <- 0
iter <- 0
while ((n <= cv) & iter < nrow(coords)){
  sample_i <- sample(1:nrow(coords), 1, replace = FALSE)
  if ((!(sample_i %in% neighbor_map)) & (sum(full_index[sample_i,] %in% test_cols) == 0)){
    test_cols <- c(test_cols,sample_i)
    neighbor_map <- c(neighbor_map,sample_i,full_index[sample_i,])
    n <- n + 1
  }
  iter <- iter + 1
}

compare_neighbors <- files %>% 
  purrr::map_dfr(.f = function(.x){
    print(.x)
    result <- readRDS(paste0(path,.x))
    model <- result$model
    u <- model$u
    v <- model$v
    x <- u%*%t(v)
    pwd <- rep(0,cv)
    for (i in 1:cv){
      hood <- c(test_cols[i], full_index[test_cols[i],])
      hood_x <- t(x[,hood])
      pwd[i] <- sum(dist(hood_x))
    }
    tibble(sample = unique(result$sample), lambda = unique(result$lambda)) %>% 
      bind_cols(hood = 1:10) %>% 
      bind_cols(dist = pwd)
  })

compare_neighbors %>% 
  ggplot(aes(x = lambda, y = dist)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = scales::log_trans()) +
  facet_wrap(~hood)
