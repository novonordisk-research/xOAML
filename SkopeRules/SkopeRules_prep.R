# Code to prepare the data for SkopeRules
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############
library(Seurat)
library(tidyverse)
library(data.table)


################# Load cluster data ###############
# Data generated with Clustering.R

results_path <- "clustering_output_path/"
cluster_data <- readRDS(paste0(results_path, "clustered_data.rds"))

################# Generate the data for SkopeRule ###############

# Save cluster IDs
cluster_data_temp <- cluster_data
n_clu = 0.5
cluster_data_temp$cluster_ID <- readRDS(paste0(results_path, "silhouette_", n_clu, ".rds")) %>%
  as.data.frame.matrix() %>%
  dplyr::pull(cluster)
table(cluster_data_temp$cluster_ID)

cluster_data_temp@meta.data %>%
  dplyr::select(eid,cluster_ID) %>%
  write.table(.,
              "skoperule_input/Clin_model_cluster_ID.tsv",
              row.names = F,
              col.names = T,
              quote = F,
              sep = "\t"
  )


# Save raw values
RAW_values  <- cluster_data@meta.data %>%
  dplyr::select(contains("RAW_"))

write.table(RAW_values,
            "skoperule_input/skrules_raw_values.tsv",
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")


# Save phenotypes
write.table(cluster_data@meta.data[,"obs",drop=F],
            "skoperule_input/skrules_obs.tsv",
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")

