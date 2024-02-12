# Code to generate clusters (multiple resolutions tested)
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############
library(Seurat)
library(tidyverse)
library(data.table)
library(ggplot2)
library(pbapply)
library(parallel)
`%>%` <- magrittr::`%>%`
pdf(file = NULL)
set.seed(42)

# from https://github.com/rbpatt2019/chooseR
source("cluster.stability/R/pipeline.R")


################# Load ML output data to generate a seurat object ###############
################# The object contains SHAP value as counts, and SHAP and RAW values as metadata for future uses ###############




# Load XGBoost model output files: SHAP values
# Rows = individuals (rownames = UKBB eids)
# Columns = features (colnames = feature names)
# Values = XGBoost output SHAP values
SHAP_values <- fread("xgboost_output_path/SHAP_values_CV5_testset.txt", header = T, sep = "\t")
SHAP_values <- as.data.table(SHAP_values)

# Load XGBoost model output files: RAW values
# Rows = individuals (rownames = UKBB eids)
# Columns = features (colnames = feature names)
# Values = XGBoost input RAW values
RAW_values <- fread("xgboost_output_path/SHAP_raw_data_values_CV5_testset.txt",header = T, sep = "\t")
RAW_values <- as.data.frame(lapply(RAW_values, function(x) as.numeric(as.character(x))))

# Appending column names to differentiate SHAP (output) and RAW (input) values
colnames(SHAP_values) <- paste0("SHAP_",colnames(SHAP_values))
colnames(RAW_values) <- paste0("RAW_",colnames(RAW_values))

# Load XGBoost model output files: prediction values
# Rows = individuals
# Columns = eid, obs, pred_class, pred_pob 
# Example: Fig_1b/Source_Data_Fig1B_Individual_level_pred.txt
Individual_level_pred <- fread("xgboost_output_path/Individual_level_pred.txt", header = T, sep = "\t")
head(Individual_level_pred)


# Seurat object metadata (SHAP + RAW)

ALL_DATA <- cbind(SHAP_values,RAW_values, Individual_level_pred) %>%
  as.data.frame()
rownames(ALL_DATA) <- paste0("ind_",c(1:dim(ALL_DATA)[1]))

# Seurat object counts (SHAP)
ALL_DATA_mat <- t(ALL_DATA[,colnames(ALL_DATA) %in% colnames(SHAP_values)]) %>%
  as.data.frame()

# Create Seurat object 
SHAP_seurat <- CreateSeuratObject(ALL_DATA_mat, meta.data = ALL_DATA)
head(SHAP_seurat@meta.data)

# Run PCA and UMAP
SHAP_seurat <- SHAP_seurat%>%
  #NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(do.scale = F,
            do.center = F) %>% 
  RunPCA(npcs = 40, seed.use = 42) 

# DimPlot(SHAP_seurat,reduction="pca")
# ElbowPlot(SHAP_seurat, ndims = 39)
npcs <- 10
SHAP_seurat <- SHAP_seurat %>% 
  RunUMAP(1:npcs, seed.use = 42)  
SHAP_seurat <- SHAP_seurat %>% 
  FindNeighbors(dims = 1:npcs)



################# Load ML output data to generate a seurat object ###############

# Adapted from https://github.com/rbpatt2019/chooseR/blob/master/examples/1_seurat_pipeline.R
# Set all resolutions to be tested, with 100 repeat, in parallel

obj <- SHAP_seurat
DefaultAssay(obj) <- "RNA"
n_repeat <- 100
assay <- "RNA"
reduction <- "pca"
results_path <- "clustering_output_path/"
dir.create(results_path)

resolutions <- c(0.01,0.05,seq(0.1,2,0.1))
resolutions <- resolutions %>%
  as.numeric() %>%
  round(digit=2) %>%
  as.numeric()
length(resolutions)


n_par <- 4
cl <- parallel::makeCluster(n_par,type="FORK")

function_test = function(x){
  
  obj <- find_clusters(
    obj,
    reduction = reduction,
    npcs = npcs,
    assay = assay,
    resolution = x
  )
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{x}")]]
  
  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
    obj,
    n = n_repeat,
    size = 0.8,
    npcs = npcs,
    res = x,
    reduction = reduction,
    assay = assay
  )
  
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }
  
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )
  
  # Now calculate silhouette scores
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", x, ".rds"))
  
  # Finally, calculate grouped metrics
  message(paste0("Grouping ", x, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", x, ".rds"))
  sil <- group_sil(sil, x)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", x, ".rds"))
  
  
}

res_list <- pblapply(resolutions, function_test, cl = cl)

parallel::stopCluster(cl)

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))
