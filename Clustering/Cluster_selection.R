# Code to pick the clustering resolution 
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############
library(Seurat)
library(tidyverse)
library(data.table)
library(caret)


################# Load data ###############
# Data generated with Clustering.R

results_path <- "clustering_output_path/"
resolutions <- list.files(results_path,pattern = "silhouette_grouped_") %>%
  gsub("silhouette_grouped_|.rds","",.) %>%
  as.numeric() %>%
  round(digit=3) %>%
  as.numeric()
resolutions


################# Silhouette scores per cluster resolution ###############

scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()

scores$res <- scores$res %>% # fix for weird issue
  as.numeric() %>%
  round(digit=3) %>%
  as.numeric()


################# PPV, sensitivity and F1 per cluster resolution ###############
# Data generated with Clustering.R

cluster_data <- readRDS(paste0(results_path, "clustered_data.rds"))

list_res <- list()

for (i in resolutions) {
  
  list_res[[as.character(i)]] <- readRDS(paste0(results_path, "silhouette_", i, ".rds")) %>%
    as.data.frame.matrix() %>%
    dplyr::select(cluster) %>%
    cbind(cluster_data@meta.data,.) %>%
    group_by(cluster) %>%
    mutate(PPV = posPredValue(factor(pred_class), factor(obs), positive="1"))  %>%
    mutate(sensitivity = sensitivity(factor(pred_class), factor(obs), positive="1")) %>%
    mutate(F1 = (2 * PPV * sensitivity) / (PPV + sensitivity)) %>%
    mutate(res = i) %>%
    dplyr::select(res,cluster,PPV,sensitivity,F1) 
  
}

res_all <- do.call(rbind,list_res)

res_all <- res_all %>%
  group_by(res,cluster) %>%
  mutate(N_clust = n()) %>%
  distinct() %>%
  as.data.frame()

res_all$PPV[is.nan(res_all$PPV)] <- 0
res_all$F1[is.nan(res_all$F1)] <- 0

################# Average metrics per cluster resolution ###############

wt_mean_scores <- res_all %>%
  pivot_longer(-c(res,cluster,N_clust), names_to = "var", values_to = "values") %>%
  group_by(var,res) %>%
  summarise(wt_mean_score = matrixStats::weightedMedian(values, N_clust, na.rm=TRUE))

avg_silhouettes <- scores %>%
  group_by(res) %>%
  left_join(res_all) %>%
  summarise(wt_mean_sil = Hmisc::wtd.mean(avg_sil, N_clust, na.rm=TRUE),
            n_clusters = n()) %>%
  arrange(-wt_mean_sil)

################# Average metrics per cluster resolution ###############

# The optimal resolution parameter was selected manually (with the table below + supplementary figure 5)
# to maximise cluster robustness, number of clusters and per-cluster F1 values 
# (resolution = 0.5, n = 14 clusters)
wt_mean_scores %>%
  left_join(avg_silhouettes) %>%
  mutate(cluster_metric =  wt_mean_sil / (1-wt_mean_score)) %>%
  mutate(var = str_to_title(var)) %>% 
  filter(var == "F1") %>%
  arrange(-cluster_metric,-n_clusters)

optimal_cluster <- wt_mean_scores %>%
  left_join(avg_silhouettes) %>%
  mutate(cluster_metric =  wt_mean_sil / (1-wt_mean_score)) %>%
  mutate(var = str_to_title(var)) %>%
  filter(res == 0.5)

optimal_cluster


