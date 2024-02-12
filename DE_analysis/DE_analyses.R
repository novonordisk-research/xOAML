# Code to run the differential expression (DE) analysis
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############
library(Seurat)
library(tidyverse)
library(data.table)

################# Load proteomics data ###############

# Load the UKBB Olink data (same processing as used in the ML model / described in methods)
# Rows = genes (rownames = gene names)
# Columns = individuals (colnames = UKBB eids)
# Values = Olink NPX values
UKBB_prot_data = fread("proteomics_input/Data_Olink.txt") %>%
  as.data.frame()

################# Load cluster data ###############
# Data generated with Clustering.R

results_path <- "clustering_output_path/"
cluster_data <- readRDS(paste0(results_path, "clustered_data.rds"))

################# Create Seurat Object with Olink data and necessary metadata ###############

# Only keep samples that have both eids in proteomics and the clustered data
cluster_data_temp <- cluster_data
cluster_data_temp@meta.data <- cluster_data_temp@meta.data %>%
  mutate(eid = as.character(eid))
eid_to_keep <- intersect(cluster_data_temp@meta.data$eid,colnames(UKBB_prot_data))
length(eid_to_keep)

# Subset datasets
UKBB_prot_data_subset <- UKBB_prot_data[match(eid_to_keep,colnames(UKBB_prot_data))]
cluster_data_subset <- subset(cluster_data_temp,subset = eid %in% eid_to_keep)
rownames(cluster_data_subset@meta.data) <- cluster_data_subset@meta.data$eid

# Safety checks
dim(UKBB_prot_data_subset)
cluster_data_subset
mean(cluster_data_subset@meta.data$eid %in% colnames(UKBB_prot_data))

# Create Seurat object with the proteomics data for DE
seurat_obj <- CreateSeuratObject(UKBB_prot_data_subset, meta.data = cluster_data_subset@meta.data)
head(seurat_obj@meta.data)
table(seurat_obj@meta.data$cluster_ID)


################# DE analyses ###############

Idents(seurat_obj) <- "cluster_ID"

# Cases only
# Adjusted for the Sex variable
seurat_obj_cases <- subset(seurat_obj, obs == 1)
table(seurat_obj_cases$cluster_ID)
seurat_DE_cases <- FindAllMarkers(seurat_obj_cases,   
                                  min.pct = 0,
                                  logfc.threshold = 0,
                                  only.pos = FALSE,
                                  random.seed = 42,
                                  return.thresh = 1,  
                                  test.use = "LR",
                                  latent.vars = c("RAW_Sex"))


################# Save output ###############

saveRDS(seurat_DE_cases,
        "DE_output/seurat_DE_cases.rds")
