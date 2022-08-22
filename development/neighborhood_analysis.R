set.seed(1)
input <- "151507"
path <- "G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_21_04/output/"
files <- list.files(path)
files <- files[stringr::str_detect(files,".rds")]
true_data <- getRDS("2020_maynard_prefrontal-cortex", input)
cv <- 100
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
      bind_cols(hood = 1:cv) %>% 
      bind_cols(dist = pwd)
  })

p1 <- compare_neighbors %>% 
  mutate(lambda = as.factor(lambda)) %>% 
  ggplot(aes(x = lambda, y = dist)) +
  geom_boxplot() +
  scale_y_continuous(trans = scales::log_trans()) +
  ggtitle(paste("Pairwise Distance of",cv,"neighborhoods"))


lambdas <- compare_neighbors %>% select(lambda) %>% distinct() %>% pull() %>% sort()
lambdas <- lambdas[1:6]
clust_fun <- c("kmeans","walktrap", "louvain", "leading_eigen", "bayesspace_t", "bayesspace_n", "true")

combos <- expand.grid(lambdas = lambdas, func = clust_fun)
combos <- combos[c(1:24, 25,31,37 ),]

plots <<- NULL
clust_results <- list(combos$lambdas, as.character(combos$func)) %>% 
  purrr::pmap_dfr(.f = function(.x, .y){
    print(.y)
    res <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_21_04/output/clust_151507_20_",.x,"_NULL_1.rds"))
    model <- res$model
    data <- model$v
    
    
    if(.y == "kmeans"){
      clust <- kmeans(data,centers)$cluster
    } else if(.y == "true"){
      clust <- true
    } else if(.y == "bayesspace_t"){
      clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_t.rds"))
    } else if(.y == "bayesspace_n"){
      clust <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_17_01/output/clust_",input,"_1_norm.rds"))
    } else {
      clust <- as.numeric(clusterRows(data, BLUSPARAM=NNGraphParam(shared = FALSE, cluster.fun = .y)))
      
    }
    ari <- safe_ari(clust,clue::solve_LSAP(table(true, clust), maximum = TRUE)[true])
    
    if(!is.null(ari$result)){
      p <- clusterPlot(true_data, label=clust, palette=NULL, size=0.05) +
        scale_fill_viridis_d(option = "A") +
        labs(title=paste("lambda=",.x,"| cluster=",.y,"| ari =",round(ari$result, 4))) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
      if(length(plots) < which(clust_fun == .y)){
        plots[[which(clust_fun == .y)]] <<- list(p)
      } else {
        plots[[which(clust_fun == .y)]] <<- c(plots[[which(clust_fun == .y)]],list(p))
      }
      
    } 
    
    return(tibble(lambda = .x, cluster = .y, ari = ari$result))
  })

pdf(file = paste0(path,input,".pdf"), width = 12, height = 12)
do.call("grid.arrange", c(list(p1), ncol=1))
for(i in 1:length(clust_fun)){
  do.call("grid.arrange", c(plots[[i]], ncol=3))
}
dev.off()
