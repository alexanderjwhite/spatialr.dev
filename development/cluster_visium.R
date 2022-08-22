library(bluster)
library(mclust)
library(dplyr)
library(BayesSpace)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(Seurat)
library(smoothLRC)

names <- paste(1:11)
pre <- c("before", "after")
n_clust <- 2:8
grid <- expand.grid(names = names, pre = pre)[-18,]

files <- list.files("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_08_17_01/output/", full.names = TRUE)
main_out <- "G:/My Drive/Dissertation/Spatial/spatial_analysis/visium/"

for (step in 1:nrow(grid)){
  
  print(grid[step,])
  coords <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/data_08_16_2022/data_",grid[step,"names"],"_",grid[step,"pre"],"_coords.rds"))
  coords$row <- as.numeric(coords$row)
  coords$col <- as.numeric(coords$col)

  
  file <- files[stringr::str_detect(files, paste0("_",grid[step,"names"],"_",grid[step,"pre"]))]
  res_path <- paste0("data_",grid[step,"names"])
  out_path <- paste0(main_out,res_path)
    if (!file.exists(out_path)){
      dir.create(file.path(out_path))
    } 
  model <- readRDS(file)$model
  saveRDS(model, paste0(out_path,"/",grid[step,"pre"],"_model.rds"))
  data <- model$v
  plots <<- NULL
  n_clust %>% 
    purrr::walk(.f = function(.x){
      set.seed(1)
      clust <- Mclust(data, G = 2:3, modelNames = "EEE")$classification
      saveRDS(clust, paste0(out_path,"/",grid[step,"pre"],"_clust.rds"))
      plot_data <- coords
      plot_data <- cbind(plot_data, data.frame(cluster = as.factor(clust)))
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = col, y = -row, color = cluster)) +
        ggplot2::geom_point(size = 0.25) +
        ggplot2::scale_color_viridis_d() +
        ggplot2::labs(title=paste("k =",.x)) +
        ggplot2::theme(axis.title.x = element_blank(),
                       axis.title.y = element_blank())
      plots <<- c(plots,list(p))
      
      
    })
  pdf(file = paste0(out_path,"/",grid[step,"pre"],"_clust.pdf"), width = 10, height = 10)
  do.call("grid.arrange", c(plots, ncol=2))
  dev.off()
}
