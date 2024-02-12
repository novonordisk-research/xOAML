# Script that utilises derived SHAP values and raw values in the ML models to 
# i) generate SHAP beeswarm plot 
# ii) generate ranked list of most important features for prediction
# SHAP beeswarm plots for ClinSNP, ClinGRS, ClinPath, ClinMet and ClinPro (top 40) were maually inspected and documented for Figure 7A
# See examples for required data format in 
# SourceData7C_SHAP_raw_data_values_CV5_testset.txt -> load as SHAP_raw_data_values_CV5_testset.txt
# SourceData7C_SHAP_values_CV5_testset.txt -> load as SHAP_values_CV5_testset.txt



# set dir to the specifc model of interest e.g. ClinPro
set.seed(432)
dir <- ()
setwd(dir) 

###### load library ###
library(haven)
library(datasets)
library(caret)
library(data.table)
library(pROC)
library(ggplot2)
library(plyr)
library(lubridate)
library(Rmisc)
library(SHAPforxgboost)
library(xgboost)
library(here)
library(tidyr)
library(tidyverse)
library(parallel)
library(shapviz)
source("SHAP_functions_cleaned.R")


## load data
combined_shap_values <- read.table(paste0(dir,"SHAP_values_CV5_testset.txt"), header = T, sep = "\t")
combined_raw_values <- read.table(paste0(dir,"SHAP_raw_data_values_CV5_testset.txt"),header = T, sep = "\t")
combined_raw_values <- as.data.frame(lapply(combined_raw_values, function(x) as.numeric(as.character(x))))
combined_shap_values <- as.data.table(combined_shap_values)
mean_shap_score <- colMeans(abs(combined_shap_values))[order(colMeans(abs(combined_shap_values)), decreasing = T)]
shap_results <- list(shap_score = combined_shap_values, mean_shap_score = (mean_shap_score))


#print SHAP feature importance 
png(filename="SHAP_mean_top40_test_.png") 
var_importance(shap_results, top_n=40) 
dev.off()


# plot SHAP beeswarm
shap_long = shap.prep(shap = shap_results,
                      X_train = combined_raw_values,
                      top_n = 40) 

png(filename = "SHAP_beeswarm_top40.png",width = 960, height = 750) 
plot.shap.summary(data_long = shap_long)
dev.off()

shap_long = shap.prep(shap = shap_results,
                      X_train = combined_raw_values,
                      top_n =length(which(shap_results$mean_shap_score > 0)))

write.table(as.data.frame(table(shap_long$variable)),"Ranking_features.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")
