library(spatialLIBD)
library(BayesSpace)
library(mclust)
library(ggplot2)
library(dplyr)
set.seed(1)
safe_ari <- purrr::safely(mclust::adjustedRandIndex)
safe_lsap <- purrr::safely(clue::solve_LSAP)

samples <- c("151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676")
sce <- spatialLIBD::fetch_data(type="sce")
compile_results <- samples %>% 
  purrr::map_dfr(.f = function(.x){
    print(.x)
    sample <- .x
    dlpfc <- getRDS("2020_maynard_prefrontal-cortex", sample)
    true <- as.numeric(dlpfc$layer_guess_reordered)
    cell_names <- colnames(dlpfc)
    dec <- scran::modelGeneVar(dlpfc)
    top <- scran::getTopHVGs(dec, n = 2000)
    
    dlpfc <- scater::runPCA(dlpfc, subset_row=top)
    dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)
    
    if(sample %in% c("151669", "151670", "151671", "151672")){
      centers <- 5
    } else {
      centers <- 7
    }
    
    # q <- 7  # Number of clusters
    d <- 15  # Number of PCs
    clusters <- list()
    
    Y <- reducedDim(dlpfc, "PCA")[, seq_len(d)]
    
    km <- kmeans(Y, centers = centers)$cluster
    shuffled <- clue::solve_LSAP(table(km, true), maximum = TRUE)[km]
    clusters[["k-means"]] <- shuffled
    
    
    g.jaccard = scran::buildSNNGraph(dlpfc, use.dimred="PCA", type="jaccard")
    lv <- igraph::cluster_louvain(g.jaccard)$membership
    n_mem <- unique(lv)
    if(length(n_mem) > centers){
      tbl <- table(lv,true)
      for(j in (centers+1):length(n_mem)){
        tbl <- cbind(tbl,rep(0,length(n_mem)))
      }
      shuffled <- clue::solve_LSAP(tbl, maximum = TRUE)[lv]
    } else {
      shuffled <- clue::solve_LSAP(table(lv,true), maximum = TRUE)[lv]
    }
    
    clusters[["Louvain"]] <- shuffled
    
    
    mc <- Mclust(Y, G = centers, modelNames = "EEE")$classification
    shuffled <- clue::solve_LSAP(table(mc, true), maximum = TRUE)[mc]
    clusters[["mclust"]] <- shuffled
    
    # clusters[["stLearn (HVGs)"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".stLearn_HVGs.csv") , package = "BayesSpace"))$pca_kmeans + 1
    # clusters[["stLearn (markers)"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".stLearn_markers.csv"), package = "BayesSpace"))$pca_kmeans + 1
    # clusters[["Giotto"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".Giotto_HMRF.csv"), package = "BayesSpace"))$HMRF_km_Delaunay_k7_b.9
    
    clusters[["stLearn"]] <- read.table(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/stlearn/",sample,".tsv"), fill = TRUE, header = TRUE)$X_pca_kmeans+1
    giotto_results <- read.table(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/Giotto/",sample,".tsv"), sep='\t', header=TRUE)
    giotto_clust <- giotto_results %>% 
      tibble() %>% 
      select(cell_ID, HMRF_cluster)
    clusters[["Giotto"]] <- tibble(cell_ID = cell_names) %>% 
      left_join(giotto_clust, by = "cell_ID") %>% 
      pull(HMRF_cluster)
    
    ## Fetch original SCE with cluster labels
    # sce <- spatialLIBD::fetch_data(type="sce")
    #> Loading objects:
    #>   sce
    # sce <- sce[, sce$sample_name == sample]

    ## Walktrap and ground truth (manual annotations)
    
    wt_hvg <- sce[, sce$sample_name == sample]$HVG_PCA_spatial
    n_mem <- unique(wt_hvg)
    if(length(n_mem) > centers){
      tbl <- table(wt_hvg,true)
      for(j in (centers+1):length(n_mem)){
        tbl <- cbind(tbl,rep(0,length(n_mem)))
      }
      shuffled <- clue::solve_LSAP(tbl, maximum = TRUE)[wt_hvg]
    } else {
      shuffled <- clue::solve_LSAP(table(wt_hvg,true), maximum = TRUE)[wt_hvg]
    }
    
    
    clusters[["Walktrap (HVGs)"]] <- shuffled
    
    wt_mrk <- sce[, sce$sample_name == sample]$markers_PCA_spatial
    n_mem <- unique(wt_mrk)
    if(length(n_mem) > centers){
      tbl <- table(wt_mrk,true)
      for(j in (centers+1):length(n_mem)){
        tbl <- cbind(tbl,rep(0,length(n_mem)))
      }
      shuffled <- clue::solve_LSAP(tbl, maximum = TRUE)[wt_mrk]
    } else {
      shuffled <- clue::solve_LSAP(table(wt_mrk,true), maximum = TRUE)[wt_mrk]
    }
    
    
    clusters[["Walktrap (markers)"]] <- shuffled
    
    
    clusters[["BayesSpace (n)"]] <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/bayesspace/clust_",sample,"_1_norm.rds"))
    clusters[["BayesSpace (t)"]] <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/bayesspace/clust_",sample,"_1_t.rds"))
    
    clusters[["smoothLRC"]] <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/smoothlrc/",sample,".rds"))
    
    clusters[["Manual annotation"]] <- sce$layer_guess_reordered
    
    clusters %>% 
      purrr::map_dfr(.id = "method", .f = function(.x){
        clust <- as.numeric(.x)
        shuffled_assign <- safe_lsap(table(clust, true), maximum = TRUE)$result
        if(is.null(shuffled_assign)){
          shuffled_assign <- safe_lsap(table(true, clust), maximum = TRUE)$result
          shuffled <- shuffled_assign[true]
        }
        shuffled <- shuffled_assign[clust]
        
        ari <- safe_ari(true,shuffled)
        tibble(input = sample, part = list(shuffled), ari = ari$result)
      })
    
  })

sp <- c("smoothLRC", "BayesSpace (n)", "BayesSpace (t)", "Giotto", "stLearn")

p <- compile_results %>% 
  filter(method != "Manual annotation") %>% 
  filter(method != "Walktrap (markers)") %>% 
  mutate(method = ifelse(method == "Walktrap (HVGs)", "Walktrap", method)) %>% 
  mutate(sp = ifelse(method %in% sp, "Spatial", "Non-Spatial")) %>% 
  mutate(method = factor(method, levels = c("k-means",
                                            "mclust",
                                            "Walktrap (markers)",
                                            "Walktrap",
                                            "Louvain",
                                            "Giotto",
                                            "stLearn",
                                            "BayesSpace (n)",
                                            "BayesSpace (t)",
                                            "smoothLRC"))) %>% 
  ggplot(aes(x = method, y = ari, fill = sp)) +
  geom_violin() +
  geom_point(size = 2.25) +
  theme_bw() +
  # theme(text = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 15, hjust=1), legend.position = "none",text = element_text(size = 12)) +
  scale_fill_viridis_d() +
  # ggtitle("Twelve-sample clustering accuracy acomparison") +
  xlab("") +
  ylab("ARI")

compile_results <- readRDS("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/compiled.rds")
# saveRDS(compile_results, "G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/compiled.rds")




pdf(file = paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/figures/ari_compare.pdf"), width = (6*mult), height = (3.75*mult))
p
dev.off()

# dlpfc <- getRDS("2020_maynard_prefrontal-cortex", sample)
# true <- as.numeric(dlpfc$layer_guess_reordered)
# cell_names <- colnames(dlpfc)
# dec <- scran::modelGeneVar(dlpfc)
# top <- scran::getTopHVGs(dec, n = 2000)
# 
# dlpfc <- scater::runPCA(dlpfc, subset_row=top)
# dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)
# 
# q <- 7  # Number of clusters
# d <- 15  # Number of PCs
# clusters <- list()
# 
# Y <- reducedDim(dlpfc, "PCA")[, seq_len(d)]
# 
# clusters[["k-means"]] <- kmeans(Y, centers = q)$cluster
# g.jaccard = scran::buildSNNGraph(dlpfc, use.dimred="PCA", type="jaccard")
# clusters[["Louvain"]] <- igraph::cluster_louvain(g.jaccard)$membership
# clusters[["mclust"]] <- Mclust(Y, G = 7, modelNames = "EEE")$classification
# 
# # clusters[["stLearn (HVGs)"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".stLearn_HVGs.csv") , package = "BayesSpace"))$pca_kmeans + 1
# # clusters[["stLearn (markers)"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".stLearn_markers.csv"), package = "BayesSpace"))$pca_kmeans + 1
# # clusters[["Giotto"]] <- read.csv(system.file("extdata", "2020_maynard_prefrontal-cortex", paste0(sample,".Giotto_HMRF.csv"), package = "BayesSpace"))$HMRF_km_Delaunay_k7_b.9
# 
# clusters[["stLearn"]] <- read.table(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/stlearn/",sample,".tsv"), fill = TRUE, header = TRUE)$X_pca_kmeans+1
# giotto_results <- read.table(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/Giotto/",sample,".tsv"), sep='\t', header=TRUE)
# giotto_clust <- giotto_results %>% 
#   tibble() %>% 
#   select(cell_ID, HMRF_cluster)
# clusters[["Giotto"]] <- tibble(cell_ID = cell_names) %>% 
#   left_join(giotto_clust, by = "cell_ID") %>% 
#   pull(HMRF_cluster)
# 
# ## Fetch original SCE with cluster labels
# sce <- spatialLIBD::fetch_data(type="sce")
# #> Loading objects:
# #>   sce
# sce <- sce[, sce$sample_name == sample]
# 
# ## Walktrap and ground truth (manual annotations)
# clusters[["Walktrap (HVGs)"]] <- sce$HVG_PCA_spatial
# clusters[["Walktrap (markers)"]] <- sce$markers_PCA_spatial
# 
# 
# clusters[["BayesSpace (n)"]] <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/bayesspace/clust_",sample,"_1_norm.rds"))
# clusters[["BayesSpace (t)"]] <- readRDS(paste0("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/bayesspace/clust_",sample,"_1_t.rds"))
# 
# clusters[["Manual annotation"]] <- sce$layer_guess_reordered
# 
# clusters %>% 
#   purrr::map_dfr(.id = "method", .f = function(.x){
#     clust <- as.numeric(.x)
#     shuffled_assign <- safe_lsap(table(clust, true), maximum = TRUE)$result
#     if(is.null(shuffled_assign)){
#       shuffled_assign <- safe_lsap(table(true, clust), maximum = TRUE)$result
#       shuffled <- shuffled_assign[true]
#       }
#     shuffled <- shuffled_assign[clust]
#     
#     ari <- safe_ari(true,shuffled)
#     tibble(part = list(shuffled), ari = ari$result)
#   })
# clusters[["stLearn"]] %>% unique() %>% length()
