library(dplyr)
library(devtools)
load_all()
obj <- readRDS("G:/My Drive/Dissertation/Spatial/spatial_analysis/cv_2022_06_14_02/debug.rds")

cv_like(obj$test_x, obj$u, obj$v, obj$test_index, obj$index_map)
