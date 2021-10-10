library(dplyr)
library(ggplot2)
result <- NULL
unrun <- NULL
lam_seq <- c(0, 0.01, 0.03, 0.08, 0.22, 0.60, 1.67, 4.64, 12.92, 35.94, 100)
comp_seq <- c("none", "u", "v")
safe_read <- purrr::safely(readRDS)
for(lam in lam_seq){
  for(comp in comp_seq){
    file <- safe_read(paste0("G:\\My Drive\\Dissertation\\Spatial\\spatialr.dev\\development\\exp_2021_10_09\\spatial_",lam,"_",comp,".rds"))
    if(!is.null(file$result)){
      file$result$comp <- comp
      result <- c(list(file$result), result)
    } else {
      unrun <- c(unrun, paste(comp,lam))
    }
  }
}

1:length(result) %>% 
  purrr::map_dfr(.f = function(.x){
    mem <- result %>% purrr::pluck(.x)
    tibble(id = .x, comp = mem$comp, lam = mem[[1]]$lambda, lik = mem[[1]]$lik)
  }) %>% 
  arrange(comp, lik) 
  ggplot(aes(x = log(lam), y = lik, color = comp)) +
  geom_point() +
  geom_line()

W <- readRDS("development/exp_2021_10_09/spatial_w.rds")
NMF::aheatmap(result[[22]][[1]]$u[which(W[2,]==1),], Colv = FALSE, Rowv = FALSE)


for(i in 1:length(result)){
  
}
