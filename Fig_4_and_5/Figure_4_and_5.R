# Code to generate Figure 4 and 5
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############

library(Seurat)
library(tidyverse)
library(data.table)
library(caret)
library(ComplexHeatmap)
library(ggpmisc)
library(cowplot)
library(ggcorrplot)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

sessionInfo()

################# Load cluster data ###############
# Data generated with Clustering.R

results_path <- "clustering_output_path/"
cluster_data <- readRDS(paste0(results_path, "clustered_data.rds"))

# Selected cluster resolution for analyses
n_clu = 0.5
cluster_data$cluster_ID <- readRDS(paste0(results_path, "silhouette_", n_clu, ".rds")) %>%
  as.data.frame.matrix() %>%
  dplyr::pull(cluster)

################# Define average prediction metrics per clusters ###############

cluster_data@meta.data$PPV_pred <-  cluster_data@meta.data %>%
  group_by(cluster_ID) %>%
  mutate(PPV = posPredValue(factor(pred_class), factor(obs), positive="1")) %>%
  pull(PPV)

cluster_data@meta.data$sensitivity_pred <-  cluster_data@meta.data %>%
  group_by(cluster_ID) %>%
  mutate(sensitivity = sensitivity(factor(pred_class), factor(obs), positive="1")) %>%
  pull(sensitivity)

cluster_data@meta.data$F1_pred <-  cluster_data@meta.data%>%
  group_by(cluster_ID) %>%
  mutate(F1 = (2 * PPV_pred * sensitivity_pred) / (PPV_pred + sensitivity_pred)) %>%
  pull(F1)

cluster_data$PPV_pred[is.nan(cluster_data$PPV_pred)] <- 0
cluster_data$F1_pred[is.nan(cluster_data$F1_pred)] <- 0


# Data wrangling
data_to_plot <- cbind(cluster_data@reductions$umap@cell.embeddings,
                      cluster_data@meta.data)  %>%
  mutate(eid = as.character(eid)) %>%
  dplyr::select(-orig.ident,-nCount_RNA,-nFeature_RNA,-contains("RNA_snn_res")) %>%
  set_names(~ str_replace_all("RAW_", "RAW.",.) %>%
              str_replace_all("\\.", "_") %>%
              str_replace_all("RAW_", "RAW.") %>%
              str_replace_all("SHAP_", "SHAP.")) %>% 
  pivot_longer(-c(eid,cluster_ID,obs,pred_class,pred_prob,UMAP_1,UMAP_2,
                  PPV_pred,sensitivity_pred,F1_pred), 
               names_to = c("type", ".value"), 
               names_sep = '[.]')

################# UMAP plots for Figure 4 ###############

# Prepare the data for the UMAP plots
umap_plot_data <- data_to_plot %>%
  filter(type == "SHAP")

umap_plot_data_summary <- umap_plot_data %>%
  group_by(cluster_ID) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))

# UMAP plot for the clusters
G1 <- umap_plot_data %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, colour = factor(cluster_ID))) +
  geom_point(size=0.2) +
  theme_bw() +
  # theme(legend.key.size = unit(1, 'npc')) +
  labs(colour = "Cluster") +
  guides(colour=guide_legend(ncol=2,
                             override.aes = list(size = 3))) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  geom_label(data = umap_plot_data_summary,
             aes(x = UMAP_1, 
                 y = UMAP_2, 
                 fill = factor(cluster_ID),
                 label = cluster_ID),
             colour = "black",
             size = 3,
             alpha = 0.5,
             show.legend = F)

# UMAP plot for the prediction probabilities
G2 <- umap_plot_data %>%
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             z = pred_prob)) +
  stat_summary_hex(fun = mean, 
                   bins = 250, 
                   col = "black",
                   linetype = "blank",
                   linewidth = 0.00) +
  theme_bw() +
  scale_fill_gradient(low = "yellow", high = "red", na.value = NA) +
  labs(fill = "Prediction\nProbability") +
  xlab("UMAP 1") +
  ylab("UMAP 2")  +
  theme(legend.title = element_text(size=10))


# UMAP plot for the SHAP values
G3 <- umap_plot_data %>%
  dplyr::select(UMAP_1,UMAP_2,Age_OA_diagnosis_controlsmatcheddate,BMI,NSAIDs_value_pre_1yrs,walk_pace,health_rating,Vitamin_D) %>%
  dplyr::rename(Age = Age_OA_diagnosis_controlsmatcheddate,
                `NSAIDs pre 1yrs` = NSAIDs_value_pre_1yrs,
                `Walking pace` = walk_pace,
                `Health rating` = health_rating,
                `Vitamin D` = Vitamin_D) %>%
  pivot_longer(-c(UMAP_1,UMAP_2), names_to = "var", values_to = "values") %>%
  mutate(var = factor(var, levels = c("Age",
                                      "NSAIDs pre 1yrs",
                                      "BMI",
                                      "Health rating",
                                      "Vitamin D",
                                      "Walking pace"))) %>%
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             z = values)) +
  stat_summary_hex(fun = mean, 
                   bins = 250,
                   col = "black",
                   linetype = "blank",
                   linewidth = 0.00) +
  theme_bw() +
  
  scale_fill_gradientn(colours = c("blue4","blue2","grey","red1","red4"),
                       values = scales::rescale(c(-1.5,-0.75,0,0.75,1.5)),
                       limits=c(-1.5,1.5),
                       breaks = c(-1.5,-0.75,0,0.75,1.5)) +
  facet_wrap(~var, ncol =2)  +
  labs(fill = "SHAP values\n") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(legend.title = element_text(size=10)) 


################# Data prep for the various heatmap ###############

# Data wrangling
data_to_plot2 <- data_to_plot %>%
  filter(type=="RAW") %>% 
  dplyr::select(eid,cluster_ID,UMAP_1,UMAP_2,type,obs,pred_class,pred_prob,
                PPV_pred,sensitivity_pred,F1_pred,
                Age_OA_diagnosis_controlsmatcheddate,NSAIDs_value_pre_1yrs,
                BMI,Vitamin_D,health_rating,walk_pace) %>%
  pivot_longer(-c(eid,cluster_ID,UMAP_1,UMAP_2,type,obs,pred_class,pred_prob,
                  PPV_pred,sensitivity_pred,F1_pred,
                  NSAIDs_value_pre_1yrs,walk_pace,health_rating),
               names_to = "var_name",
               values_to = "value") %>%
  group_by(cluster_ID) %>%
  mutate(mean_obs = mean(obs),
         mean_pred_class = mean(pred_class),
         mean_pred_prob = mean(pred_prob),
         mean_NSAIDs_value_pre_1yrs = mean(NSAIDs_value_pre_1yrs),
         mean_walk_pace = mean(walk_pace),
         mean_health_rating = mean(health_rating),
         mean_obs = mean(obs),
         mean_pred_class = mean(pred_class)) %>%
  dplyr::select(-obs,-pred_class,-pred_prob,-NSAIDs_value_pre_1yrs,-walk_pace,-health_rating) %>%
  group_by(var_name) %>%
  mutate(value = scale(value)) %>%
  group_by(cluster_ID,var_name) %>%
  mutate(mean_value = mean(value)) %>%
  dplyr::select(-eid,-UMAP_1,-UMAP_2,-value,-type) %>%
  distinct()  %>%
  pivot_wider(names_from = "var_name",
              values_from = "mean_value") %>%
  dplyr::select(-mean_pred_class,-mean_pred_class) %>%
  column_to_rownames("cluster_ID") %>%
  arrange(-mean_pred_prob) %>%
  dplyr::select(F1_pred,PPV_pred,sensitivity_pred,mean_pred_prob,mean_obs,everything()) %>%
  # rename_all(funs(stringr::str_replace_all(., 'mean_', ''))) %>%
  dplyr::arrange(-F1_pred)


################# Load cluster rules (generated with SkopeRules) ###############
# Data generated with code from SkopeRules.py

# Data wrangling
rules_1 <- read.table("skoperules_data/OA_cluster_rules_all.txt",sep="~") %>%
  mutate(ind = rep(c(1, 2),length.out = n())) %>%
  group_by(ind) %>%
  mutate(id = row_number()) %>% 
  spread(ind, V1) %>%
  select(-id) %>%
  magrittr::set_colnames(c("cluster","temp_col"))  %>% 
  mutate(temp_col = gsub("\\)\\),.*","))]",temp_col)) %>%
  mutate(temp_col = gsub("\\[|\\]|\\(|\\)","",temp_col)) %>%
  separate(col = temp_col, sep = ",",  into = c("rules", "PPV_clusters", "sensitivity_clusters"), extra = "drop", fill = "warn") %>%
  mutate(PPV_clusters = round(as.numeric(PPV_clusters), 2),
         sensitivity_clusters = round(as.numeric(sensitivity_clusters), 2)) %>%
  as.data.frame()

################# Data wrangling ###############

# Save a copy of the dataframe before formatting, for later use
rules_1_temp <- rules_1

# Formatting
list_var_rules_1 <- rules_1$rules %>%
  str_split(pattern=" and ") %>%
  unlist() %>%
  gsub(" .*|RAW_","",.) %>%
  unique()

# Define F1 values
rules_1$F1_clusters <-  rules_1 %>%
  group_by(cluster) %>%
  mutate(F1 = (2 * PPV_clusters * sensitivity_clusters) / (PPV_clusters + sensitivity_clusters)) %>%
  pull(F1)

# Variable formatting
rules_1 <- rules_1 %>%
  select(everything(),F1_clusters,PPV_clusters,sensitivity_clusters) %>%
  mutate(rules = gsub(" and ", ", ", rules))%>%
  # mutate(rules = gsub(" and ", "\n", rules)) %>%
  mutate(rules = gsub("(.*[0-9]+\\.[0-9][0-9])[0-9]*(.*)", "\\1\\2", rules)) %>%
  mutate(rules = gsub("walk_pace > 1.5","Not-slow walking pace",rules)) %>%
  mutate(rules = gsub("walk_pace <= 1.5","Slow walking pace",rules)) %>%
  mutate(rules = gsub("NSAIDs_value_pre_1yrs > 0.5","Take NSAIDs (pre-1yr)",rules)) %>%
  mutate(rules = gsub("NSAIDs_value_pre_1yrs <= 0.5","Do not take NSAIDs (pre-1yr)",rules)) %>%
  mutate(rules = gsub("NSAIDs_value_pre_2yrs > 0.5","Take NSAIDs  (pre-2yrs)",rules)) %>%
  mutate(rules = gsub("NSAIDs_value_pre_2yrs <= 0.5","Do not take NSAIDs  (pre-2yrs)",rules)) %>%
  mutate(rules = gsub("Age_OA_diagnosis_controlsmatcheddate","Age",rules)) %>%
  mutate(rules = gsub("health_rating > 1.5","Health rating: Good or poorer",rules)) %>%
  mutate(rules = gsub("health_rating <= 1.5","Health rating: Excellent",rules)) %>%
  mutate(rules = gsub("health_rating > 2.5","Health rating: Fair or poorer",rules)) %>%
  mutate(rules = gsub("health_rating <= 2.5","Health rating: Good or better",rules)) %>%
  mutate(rules = gsub("health_rating > 3.5","Health rating: Poor",rules)) %>%
  mutate(rules = gsub("health_rating <= 3.5","Health rating: Fair or better",rules)) %>%
  mutate(rules = ifelse(grepl("Health rating: Good or poorer",rules) &
                          grepl("Health rating: Good or better",rules),
                        gsub(" or poorer| or better","",rules),rules)) %>%
  mutate(rules = ifelse(grepl("Health rating: Fair or poorer",rules) &
                          grepl("Health rating: Fair or better",rules),
                        gsub(" or poorer| or better","",rules),rules)) %>%
  mutate(rules = gsub("RAW_Health rating: Good, RAW_Health rating: Good","RAW_Health rating: Good",rules)) %>%
  mutate(rules = gsub("RAW_Health rating: Fair, RAW_Health rating: Fair","RAW_Health rating: Fair",rules)) %>%
  mutate(rules = gsub("Townsend_deprivation_index_at_recruitment.0.0","Townsend deprivation index",rules)) %>%
  mutate(rules = gsub("College_University_degree <= 0.5","Does not have College or University degree",rules)) %>%
  mutate(rules = gsub("College_University_degree > 0.5","Has College or University degree",rules)) %>%
  mutate(rules = gsub("Hand_grip_strength_.right.","Hand grip strength (Right)",rules)) %>%
  mutate(rules = gsub("IGF_1","IGF-1",rules)) %>%
  mutate(rules = gsub("LDL_direct","Direct LDL cholesterol",rules))  


# Define cluster size
cluster_sizes <- cluster_data@meta.data %>%
  group_by(cluster_ID) %>%
  tally() %>%
  mutate(n = n/sum(n)) %>%
  magrittr::set_colnames(c("cluster","cluster_proportion")) %>%
  mutate(cluster = as.character(cluster))


################# Merge data needed for heatmaps ###############

# Data wrangling
colnames(rules_1) <- c(colnames(rules_1)[1],paste0(colnames(rules_1)[-1], "_R1"))

data_to_plot3 <- data_to_plot2  %>%
  rownames_to_column("cluster") %>%
  left_join(rules_1) %>%
  left_join(cluster_sizes)%>%
  arrange(-mean_pred_prob)


################# Heatmaps for Figure 4 ###############

# Heatmap: prediction probabilities
Pred_heatmap_min <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains("mean_pred_prob")) %>%
  rename(`Average\nPrediction\nProbability` = mean_pred_prob) %>%
  t() %>%
  Heatmap(name = "Probability",
          column_title = " ",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 8),
          cluster_columns = F,
          cluster_rows = F,
          row_names_side = "left",
          col=colorRampPalette(c("grey", "cornflowerblue","purple"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 8))},
          rect_gp = gpar(col = "black", lwd = 1), 
          height = 0.2,
          show_column_names = F,
          column_names_rot = 0,
          heatmap_legend_param = list(direction = "horizontal",
                                      title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)))

# Heatmap: categorical values
Binary_heatmap <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains("mean"),-contains("obs"),-contains("prob")) %>%
  mutate(mean_walk_pace = 3 - mean_walk_pace) %>%
  mutate(mean_NSAIDs_value_pre_1yrs = range01(mean_NSAIDs_value_pre_1yrs)) %>%
  mutate(mean_walk_pace = range01(mean_walk_pace)) %>%
  mutate(mean_health_rating = range01(mean_health_rating)) %>%
  rename(`NSAIDs pre-1yr` = mean_NSAIDs_value_pre_1yrs,
         `Walking pace` = mean_walk_pace,
         `Health rating` = mean_health_rating) %>%
  t() %>%
  Heatmap(name = "Categorical values",
          column_title = "Average categorical\nvalues (rescaled)",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 8),
          cluster_columns = F,
          cluster_rows = F,
          # col = structure(c("TRUE","FALSE"), names = c("TRUE","FALSE")),
          row_names_side = "left",
          col=colorRampPalette(c("white", "yellow","red"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", .[i, j]), x, y, gp = gpar(fontsize = 8))},
          rect_gp = gpar(col = "black", lwd = 1), 
          height = 0.3,
          show_column_names = F,
          column_names_rot = 0,
          heatmap_legend_param = list(direction = "horizontal", at = seq(0, 1, 0.25),
                                      title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)))

# Heatmap: continuous values
Continuous_heatmap <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(-contains("mean"),-contains("pred"),-contains("obs"),-contains("cluster"),-contains("rules")) %>%
  # mutate_all(function(x) ifelse(x == 1,TRUE,FALSE)) %>%
  rename(Age = Age_OA_diagnosis_controlsmatcheddate,
         `Vitamin D` = Vitamin_D) %>%
  t() %>%
  Heatmap(name = "Continuous values",
          column_title = "Average continuous\nvalues (Z-scores)",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 8),
          cluster_columns = F,
          cluster_rows = F,
          show_column_names = F,
          # col = structure(c("TRUE","FALSE"), names = c("TRUE","FALSE")),
          row_names_side = "left",
          col=colorRampPalette(c("blue", "grey","red"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 8))},
          rect_gp = gpar(col = "black", lwd = 1), 
          height = 0.3,
          column_names_rot = 0,
          heatmap_legend_param = list(direction = "horizontal",
                                      title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)))

# Heatmap: cluster numbers
cluster_number_flipped <- data_to_plot3 %>%
  dplyr::select(cluster) %>%
  dplyr::rename(Cluster = cluster) %>%
  t() %>%
  Heatmap(name = "Cluster",
          column_title = " ",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F, 
          show_column_names = F,
          rect_gp = gpar(col = "black", fill = "white"),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(.[i, j], x, y,gp = gpar(fontsize = 8, fontface = "bold"))},
          height = 0.1,
          row_names_side = "left",
          show_heatmap_legend = F)



################# Extra heatmaps for Figure 4: proteomics ###############
# Data generated with DE_analyses.R

# Load proteomics DEresults
seurat_DE_cases <- readRDS("DE_output/seurat_DE_cases.rds")

# Extract top genes per clusters
gene_list <- seurat_DE_cases %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj <= 0.05) %>%
  slice_min(order_by = p_val, n=2) %>% 
  pull(gene) %>%
  unique()

# Average NPX values per clusters
summarised_markers <- seurat_obj_cases@assays$RNA@counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene,names_to = "eid",values_to = "NPX") %>%
  filter(gene %in% gene_list) %>%
  left_join(seurat_obj_cases@meta.data[,c("eid","cluster_ID")]) %>%
  group_by(gene) %>%
  mutate(NPX = scale(NPX)) %>%
  group_by(cluster_ID,gene) %>%
  summarise(NPX = mean(NPX,na.rm=T)) %>%
  pivot_wider(names_from = "gene", values_from = "NPX") %>%
  dplyr::rename(cluster = cluster_ID) %>%
  mutate(cluster = as.character(cluster))

# Heatmap annotation (add an asterisk * if significant)
Olink_annotation <- data_to_plot3 %>%
  left_join(summarised_markers)  %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains(gene_list)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cluster", values_to = "value") %>%
  left_join(seurat_DE[,c("gene","cluster","p_val_adj")]) %>%
  mutate(signif = ifelse(!is.na(p_val_adj) & p_val_adj <= 0.05, " *","")) %>%
  mutate(value = paste0(round(value,2),signif)) %>%
  dplyr::select(gene,cluster,value) %>%
  pivot_wider(names_from = "cluster", values_from = "value") %>%
  column_to_rownames("gene")

# Heatmap for the Olink data, rescale per gene
Olink_markers <- data_to_plot3 %>%
  left_join(summarised_markers)  %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains(gene_list)) %>%
  t() %>%
  Heatmap(name = "Protein levels",
          column_title = "Proteomics biomarkers (Olink)",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          column_names_gp = grid::gpar(fontsize = 10),
          row_names_gp = grid::gpar(fontsize = 8),
          cluster_columns = F,
          cluster_rows = T,
          show_column_names = F,
          show_row_dend = FALSE,
          # col = structure(c("TRUE","FALSE"), names = c("TRUE","FALSE")),
          row_names_side = "left",
          col=colorRampPalette(c("deeppink2", "grey","chartreuse3"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(Olink_annotation[i, j], x, y, gp = gpar(fontsize = 8))},
          rect_gp = gpar(col = "black", lwd = 1), 
          height = 2,
          column_names_rot = 0,
          heatmap_legend_param = list(direction = "horizontal",
                                      title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)))


# Heatmap: number of cases per clusters
case_numbers <- seurat_obj_cases@meta.data %>%
  dplyr::select(cluster_ID,obs) %>%
  group_by(cluster_ID) %>%
  tally() %>%
  dplyr::rename(`Cases w/ Olink (n)` = n) %>%
  dplyr::rename(cluster = cluster_ID) %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(select(data_to_plot3,cluster),.) %>%
  column_to_rownames("cluster") %>%
  t() %>%
  Heatmap(name = "case_N",
          column_title = " ",
          column_title_gp = gpar(fontsize = 8, fontface = "bold"),
          row_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F, 
          show_column_names = F,
          rect_gp = gpar(col = "black", fill = "white"),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(.[i, j], x, y,gp = gpar(fontsize = 8, fontface = "bold"))},
          height = 0.08,
          row_names_side = "left",
          show_heatmap_legend = F)




################# Combine all heatmaps together for Figure 4 ###############

heatmap_grob <- grid.grabExpr(draw(cluster_number_flipped %v% Pred_heatmap_min %v% Binary_heatmap %v% Continuous_heatmap %v% case_numbers %v% Olink_markers, 
                                   merge_legend = TRUE, heatmap_legend_side = "bottom", 
                                   annotation_legend_side = "bottom"))



################# Generate and save Figure 4 ###############

left_panel <- plot_grid(G1,G2, ncol = 1,labels =c("A","B"), align = "v")
right_panel <- plot_grid(G3, ncol = 1 ,labels =c("C"))
bottom_panel <- plot_grid(NULL, heatmap_grob, NULL, ncol = 3 ,labels =c("D","",""), rel_widths = c(0.05,1,0.15))

Figure4_temp <- plot_grid(left_panel,
                             right_panel,
                             ncol = 2,
                             rel_widths = c(1,1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

Figure4 <- plot_grid(Figure4_temp,
                        bottom_panel,
                        ncol = 1,
                        rel_heights = c(3,3)) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(filename="Figure4_cluster_plots_07022024.pdf", 
       plot=UMAP_plots,
       device="pdf",
       path="./figures/",
       width = 180,
       height = 180,
       units = "mm",
       scale = 1.8,
       dpi=1000)



################# Apply SkopeRules rules to the validation set ###############

# Load validation dataset 
# (UKBB data containing the features used in SkopeRules, for all individuals in the validation cohort)
# Rows = individuals
# Columns = features (colnames = feature names)
# Values = RAW values
validation_data <- fread("validation_data_path/OA_validation.txt")
validation_data_filt <- validation_data %>%
  drop_na()
colnames(validation_data_filt) <- paste0("RAW_",colnames(validation_data_filt))

# change "and" to "&" to make as a R-interpretable rule 
rules_1_temp$rules <- gsub(" and "," & ",rules_1_temp$rules) 

# Apply the rules to the validation cohort
list_res <- list()
for (i in 1:dim(rules_1_temp)[1]){
  
  list_res[[i]] <- validation_data_filt %>%
    filter_(rules_1_temp$rules[i]) %>%
    mutate(cluster = rules_1_temp$cluster[i]) %>%
    dplyr::select(RAW_eid,cluster) 
  
}

res_all <- do.call(rbind,list_res)
dim(res_all)
mean(res_all$RAW_eid %in% validation_data_filt$RAW_eid)
mean(!unique(validation_data_filt$RAW_eid) %in% unique(res_all$RAW_eid))

################# Validation statistics, for reporting ###############

# Some are not mapped
RAW_eid_none <- unique(validation_data_filt$RAW_eid)[!unique(validation_data_filt$RAW_eid) %in% unique(res_all$RAW_eid)]

# Some are mapped to multiple clusters (2 at most)
res_all %>%
  group_by(RAW_eid) %>%
  tally() %>%
  pull(n) %>%
  max()

RAW_eid_wrong <- res_all %>%
  group_by(RAW_eid) %>%
  tally() %>%
  filter(n > 1) %>%
  pull(RAW_eid)

# How many are not mapped to a cluster (%)
length(RAW_eid_none)/dim(validation_data_filt)[1]

# How many are mapped to multiple (%)
length(RAW_eid_wrong)/dim(validation_data_filt)[1]

# How many are uniquely assigned to a cluster (%)
(dim(validation_data_filt)[1] - length(RAW_eid_none) - length(RAW_eid_wrong))/dim(validation_data_filt)[1]



################# Validation results wrangling ###############

# Order cluster by case % to minimise false negatives
cluster_ordered <- data_to_plot3 %>%
  arrange(-mean_obs) %>%
  pull(cluster)

res_all %>%
  filter(RAW_eid %in% RAW_eid_wrong) %>%
  mutate(cluster = factor(cluster, levels = cluster_ordered, ordered = T)) %>%
  arrange(RAW_eid) %>%
  group_by(RAW_eid) %>%
  slice_max(order_by = desc(cluster), n=1) %>%
  group_by(cluster) %>%
  tally()

# Uniquely assign to highest risk cluster when there are 2
res_all_temp <- res_all %>%
  filter(RAW_eid %in% RAW_eid_wrong) %>%
  mutate(cluster = factor(cluster, levels = cluster_ordered, ordered = T)) %>%
  arrange(RAW_eid) %>%
  group_by(RAW_eid) %>%
  slice_max(order_by = desc(cluster), n=1) 

res_all <- res_all %>%
  filter(!RAW_eid %in% RAW_eid_wrong) %>%
  rbind(res_all_temp)

# Sanity check
dim(res_all)
dim(validation_data_filt)
dim(res_all)[1] - dim(validation_data_filt)[1]
length(RAW_eid_none)
left_join(validation_data_filt,res_all) %>%
  filter(is.na(cluster))


################# Heatmaps for Figure 5 ###############

# Discovery cohort: prediction metrics heatmap
Pred_heatmap <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(mean_pred_prob,contains("pred"),contains("obs")) %>%
  rename(F1 = F1_pred,
         PPV = PPV_pred,
         Sensitivity = sensitivity_pred,
         `Avg Pred Prob` = mean_pred_prob,
         `Case %` = mean_obs) %>%
  Heatmap(name = "OA study:\nvalues",
          column_title = "           OA study population",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F,
          row_names_side = "left",
          col=colorRampPalette(c("grey", "cornflowerblue","purple"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 10))},
          rect_gp = gpar(col = "black", lwd = 1), 
          width = 0.5)


# Discovery cohort: cluster size heatmap
Size_heatmap <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains("proportion")) %>%
  mutate(cluster_proportion = round(cluster_proportion * 100,2)) %>%
  dplyr::rename(`Cluster size (%)` = cluster_proportion) %>%
  Heatmap(name = "OA study:\ncluster\nsize",
          column_title = " ",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F,
          row_names_side = "left",
          col=colorRampPalette(c("grey", "pink","purple"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 10))},
          rect_gp = gpar(col = "black", lwd = 1), 
          width = 0.1)

# Validation cohort: prediction metrics heatmap
Pred_VAL_heatmap <- left_join(validation_data_filt,res_all) %>%
  group_by(cluster) %>%
  drop_na() %>%
  summarise(mean_obs_val = mean(RAW_PHENOTYPE,na.rm=T)) %>%
  left_join(data_to_plot3,.) %>%
  column_to_rownames("cluster") %>%
  dplyr::select(mean_obs_val) %>%
  rename(`Case %` = mean_obs_val) %>%
  Heatmap(name = "Validation population",
          column_title = " ",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F,
          row_names_side = "left",
          col=colorRampPalette(c("grey", "cornflowerblue","purple"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 10))},
          rect_gp = gpar(col = "black", lwd = 1), 
          width = 0.1,
          heatmap_legend_param = list(title = "Validation:\nvalues"))


# Validation cohort: cluster size heatmap
Size_VAL_heatmap <- left_join(validation_data_filt,res_all) %>%
  group_by(cluster) %>%
  drop_na() %>%
  summarise(N_val = n()) %>%
  mutate(Per_val =  100 *N_val/sum(N_val)) %>% 
  left_join(data_to_plot3,.) %>%
  column_to_rownames("cluster") %>%
  dplyr::rename(`Cluster size (%)` = Per_val) %>%
  dplyr::select(`Cluster size (%)`) %>%
  Heatmap(name = "Validation:\ncluster\nsize",
          column_title = "Validation population              ",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F,
          row_names_side = "left",
          col=colorRampPalette(c("grey", "pink","purple"))(100),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", .[i, j]), x, y, gp = gpar(fontsize = 10))},
          rect_gp = gpar(col = "black", lwd = 1), 
          width = 0.1)


# Cluster numbers
cluster_number <- data_to_plot3 %>%
  dplyr::select(cluster) %>%
  Heatmap(name = "Cluster",
          column_title = "Cluster",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          cluster_rows = F, 
          show_row_names = F,
          show_column_names = F,
          rect_gp = gpar(col = "black", fill = "white"),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(.[i, j], x, y,gp = gpar(fontsize = 10, fontface = "bold"))},
          width = 0.1,
          row_names_side = "left",
          show_heatmap_legend = F)


# Rules to plot alongside heatmaps
Rules_heatmap_R1 <- data_to_plot3 %>%
  column_to_rownames("cluster") %>%
  dplyr::select(contains("rules_R1")) %>%
  rename_all(funs(stringr::str_replace_all(., '_R1', ''))) %>%
  mutate(rules = gsub("RAW_","",rules)) %>%
  mutate(rules = gsub("([^,]+,[^,]+),", "\\1,\n", rules)) %>%
  mutate(rules = gsub("IGF-1 <= 34.85.*4,", "IGF-1 <= 34.85,", rules)) %>% #Manual fix
  Heatmap(name = "    ",
          column_title = " ",
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          cluster_columns = F,
          show_column_names = F,
          cluster_rows = F, 
          rect_gp = gpar(col = "black", fill = "white"),
          cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(.[i, j], x, y,gp = gpar(fontsize = 10))},
          width = 1,
          row_names_side = "left",
          show_heatmap_legend = F)



################# Combine all heatmaps together and save Figure 5 ###############

metrics_grob <- grid.grabExpr(draw(cluster_number + Pred_heatmap + Size_heatmap + Rules_heatmap_R1  + Pred_VAL_heatmap + Size_VAL_heatmap))
Figure5 <- plot_grid(metrics_grob) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))

ggsave(filename="Figure5_metrics_plots_07022024.pdf", 
       plot=metrics_plot,
       device="pdf",
       path="./figures/",
       width = 180,
       height = 150,
       units = "mm",
       scale = 1.8,
       dpi=1000)



