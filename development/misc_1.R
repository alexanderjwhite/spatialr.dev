library(ggplot2)
library(dplyr)
library(ggnewscale)
library(randomcoloR)
library(scales)
library(umap)

x <- readRDS("development/exp_bc_2021_10_15/spatial_x_sub.rds")
coords <- readRDS("development/exp_bc_2021_10_15/spatial_coords_sub.rds")

sp_svd <- sparsesvd::sparsesvd(log(x+1),rank=10)

plot_spatialClusterMap <- function(coordinates,labels){
  #coordinates is n*2
  #labels is a vector, length n
  #start plotting
  maxLabel=NULL
  maxLabel=max(labels)
  distinctColors=NULL
  distinctColors=distinctColorPalette(maxLabel)
  #    show_col(distinctColors)
  
  spatialMap_allPoints=NULL
  spatialMap_allPoints=data.frame("x"=coordinates[,1],"y"=coordinates[,2])
  
  p_spatialMap=NULL
  p_spatialMap <- ggplot(data = spatialMap_allPoints, mapping = aes(x = x, y = y))
  p_spatialMap = p_spatialMap + theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black"))+
    geom_jitter(alpha = 0.3)
  for(label_i in 1:maxLabel){
    print(paste("Label:",label_i," color:",distinctColors[label_i],sep = ''))
    cluster_coordinates=NULL
    cluster_coordinates=coordinates[(rownames(coordinates)[which(labels==label_i)]),]
    cluster_coordinates=data.frame('x'=cluster_coordinates[,1],'y'=cluster_coordinates[,2])
    
    p_spatialMap = p_spatialMap+
      geom_point(data = cluster_coordinates,color=distinctColors[label_i])
  }
  print(p_spatialMap)
}

cluster_svd <- kmeans(coords, centers=10)
# tmp <- coords
# tmp_labels <- cluster_svd$cluster
# plot_spatialClusterMap(tmp,tmp_labels)

data <- tibble(x = coords[,1], y = coords[,2], label = as.character(cluster_svd$cluster)) %>% 
  arrange(label)

data %>% 
  ggplot(aes(x = x, y = y, color = label)) +
  geom_point()

data %>% arrange(label)

u <- matrix(0,nrow = nrow(data), ncol = 10)
r_start <- 1
coef <- abs(rnorm(10, 20, sd  = 5))
for (j in 1:10) {
  n_rows <- data %>% filter(label == j) %>% nrow()
  r_end <- r_start + n_rows - 1
  cl_U <- qr.Q(qr(replicate(1, rnorm(n_rows, 0, 1))))
  u[r_start:r_end,j] <- coef[j]*cl_U
  r_start <- r_end+1
}


NMF::aheatmap(u, Rowv = NA, Colv = NA)
norm(matrix(0,nrow=3,ncol=3),type="F")
