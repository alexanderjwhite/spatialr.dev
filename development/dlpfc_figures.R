library(BayesSpace)
library(ggplot2)
library(dplyr)
library(gridExtra)
sample <- "151509"
mult <- 1.75
true_data <- getRDS("2020_maynard_prefrontal-cortex", sample)
true_label <- true_data$layer_guess_reordered
# true_levels <- levels(true_label)
# true_label <- as.character(true_label)
# true_label[is.na(true_label)] <- ""
# true_label <- factor(true_label, levels = c(true_levels,""))
true <- as.numeric(true_label)
out_path <- "G:/My Drive/Dissertation/Spatial/spatial_analysis/figures/"



# p_true <- clusterPlot(true_data, label=true_label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title="Sample 151671 manual annotation") +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())



compile_results <- readRDS("G:/My Drive/Dissertation/Spatial/spatial_analysis/final_results/compiled.rds")
sample_results <- compile_results %>% filter(input == sample) %>% 
  filter(method != "Manual annotation") %>% 
  filter(method != "Walktrap (markers)") %>% 
  mutate(method = ifelse(method == "Walktrap (HVGs)", "Walktrap", method)) %>% 
  mutate(method = factor(method, levels = c("k-means",
                                            "mclust",
                                            "Walktrap (markers)",
                                            "Walktrap",
                                            "Louvain",
                                            "Giotto",
                                            "stLearn",
                                            "BayesSpace (n)",
                                            "BayesSpace (t)",
                                            "smoothLRC")))

clust_plot <- function(sce, label){
  col_data <- colData(sce) %>% 
    as_tibble() %>% 
    select(row, col) %>% 
    bind_cols(tibble(label = label))
  
  col_data %>% 
    ggplot(aes(x = col, y = -row, color = label)) +
    geom_point(size = 0.5) +
    scale_color_viridis_d(na.translate = FALSE) +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title=element_blank(),
          text = element_text(size = 10),
          legend.position = "none",
          panel.background = element_blank())
}



# plots <- NULL
# plots <- ggplotGrob(p_true)
# plots <- c(plots,list(p_true))
p_true <- clust_plot(true_data, true_label) + labs(title=paste0(sample," annotation"))

alg <- "k-means"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p1 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p1 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")
# g2 <- ggplotGrob(p)
# plots <- c(plots,list(p))

alg <- "mclust"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p2 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p2 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "Walktrap"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1)
label[is.na(label)] <- 6
label <- as.factor(label)
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p3 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p3 <- clusterPlot(true_data, label=shuffled, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "Louvain"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1)
label[is.na(label)] <- 6
label <- as.factor(label)
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p4 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p4 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "Giotto"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p5 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p5 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "stLearn"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p6 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p6 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "BayesSpace (n)"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p7 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p7 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "BayesSpace (t)"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p8 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p8 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")

alg <- "smoothLRC"
label <- sample_results %>% filter(method==alg) %>% pull(part) %>% purrr::pluck(1) %>% as.factor()
ari <- sample_results %>% filter(method==alg) %>% pull(ari) %>% scales::number(accuracy = 0.001)
p9 <- clust_plot(true_data, label = label) + labs(title=paste0(alg," (", ari, ")"))
# p9 <- clusterPlot(true_data, label=label, palette=NULL, size=0.05, color = NA) +
#   scale_fill_viridis_d(option = "A", na.translate = FALSE) +
#   labs(title=alg) +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         legend.position = "none")


# do.call("grid.arrange", c(plots, layout_matrix = rbind(c(NA,1,NA),
#                                                        c(2,3,4),
#                                                        c(5,6,7),
#                                                        c(8,9,10))))
# do.call("grid.arrange",c(plots, layout_matrix = rbind(c(NA,1,NA),
#                                                       c(2,3,4),
#                                                       c(5,6,7),
#                                                       c(8,9))))





pdf(file = paste0(out_path,sample,".pdf"), width = (6*mult), height = (3.75*mult))
grid.arrange(p_true, p1, p2, p3, p4, p5, p6, p7, p8, p9, layout_matrix = cbind(c(1,6),
                                                                               c(2,7),
                                                                               c(3,8),
                                                                               c(4,9),
                                                                               c(5,10)))
dev.off()
