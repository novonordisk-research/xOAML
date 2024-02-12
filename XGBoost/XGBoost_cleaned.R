set.seed(432)
setwd("")

## Load libraries
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
library(parallelly)
source("SHAP_functions_cleaned.R")

############### Data load and processing #####################
# load data and rename id
# IMPORTANT!
# ID should be called 'eid'
# target outcomes need to be named 'PHENOTYPE'
# PHENOTYPE should be encoded as 0 for controls, 1 for cases.
datadir <- c("")
data_model1 <- fread(datadir)
print(table(data_model1$PHENOTYPE))

# remove features with > 0.50 data missing to mimic
missing_ratio <- sapply(data_model1, function(x) (sum(is.na(x)))/nrow(data_model1))
missing_50 <- which(missing_ratio > 0.50)
print(paste("REMOVED", colnames(data_model1)[missing_50],sep = " "))
data_model1 <- data_model1[,-missing_50]
print(table(data_model1$PHENOTYPE))

# remove features where there near0variance
feat_near0var <- colnames(data_model1)[nearZeroVar(data_model1)]
print(paste0("REMOVED DUE TO NEARZEROVAR ",feat_near0var))
data_model1 <- data_model1[,-nearZeroVar(data_model1)]
########################################################



############### Modelling parameters ###################
# Cross-validation settings
cv.fold.outer <- 5
cv.fold.inner <- 5

# hyperparameters
booster_xg <- c("gbtree")
eta_xg <- c(0.05, 0.10, 0.15, 0.20, 0.25,0.30) 
max_depth_xg <- c(10) 
n_rounds_xg <- c(50, 100, 200, 300, 500, 700, 1000) 
min_child_weight <- c(50) 

# define cross-validation for outer CV fold
index_outer <- createFolds(data_model1$PHENOTYPE, k=cv.fold.outer, returnTrain=T)

# store performance metrics per outer CV fold
outer_models_a <- matrix("NA", ncol = 9, nrow = cv.fold.outer)
colnames(outer_models_a) <- c("opt_eta_xg","opt_n_rounds_xg","ROC_AUC", "Sensitivity", "Specificity","PPV","NPV","MCC","F1")

# store inner performance of model
AUC_store <- list()

# Normalisation
cvfold.inner.median.4norm <- list() 
cvfold.outer.median.4norm <- list() 
########################################################



##################### XGboost #########################
print(paste0("Start XGBoost for controls(0)/cases(1) in XGboost"))

for (cv.out in 1:cv.fold.outer) {
  fold = paste0("Fold", cv.out)
  print(paste("*** outer CV", fold, " ***"))
  
  #### Split data outer fold #####
  training    <- data_model1[ index_outer[[fold]], ]
  training <- subset(training, select=-c(eid))
  validation  <- data_model1[-index_outer[[fold]], ]
  index_inner <- createMultiFolds(training$PHENOTYPE, k=cv.fold.inner, times = 1) # times = cv.rep.inner
  
  ## Store prediction output with eid information ##
  Individual_pred <- matrix("NA", ncol = 4, nrow = nrow(validation))
  colnames(Individual_pred) <- c("eid","obs","pred_class","pred_prob")
  Individual_pred <- as.data.frame(Individual_pred)
  Individual_pred[,1] <- validation[,c("eid")]
  Individual_pred[,2] <- validation[,c("PHENOTYPE")]
  validation <- subset(validation, select=-c(eid))
  
  # reset performance and hyperparamters 
  AUC <- c(0)
  opt_eta_xg <- NA
  opt_n_rounds_xg <-NA
  
  # Inner cross-validation
  for (cv.in in 1:cv.fold.inner) {
    index.tmp = paste("Fold", cv.in,".Rep1", sep="")
    print(paste("outer fold: ", fold, index.tmp))
    train_inner <- training[ index_inner[[index.tmp]],]
    test_inner <- training[-index_inner[[index.tmp]],]
    
    # Imputation by median value
    numdata_train_inner <- as.data.frame(sapply(train_inner, as.numeric) )
    train_inner$PHENOTYPE<- as.character(train_inner$PHENOTYPE)
    test_inner$PHENOTYPE <- as.character(test_inner$PHENOTYPE)
    median_numdata_train_inner <- sapply(numdata_train_inner, median, na.rm = T)
    cvfold.inner.median.4norm[[fold]][[index.tmp]] <- sapply(numdata_train_inner, median, na.rm = T)
    for (i in 1:ncol(train_inner)) {
      name <- colnames(train_inner)[i]
      train_inner[which(is.na(train_inner[,name]) == T),name] <- cvfold.inner.median.4norm[[fold]][[index.tmp]][[name]]
    }
    
    for (i in 1:ncol(test_inner)) {
      name <- colnames(test_inner)[i]
      test_inner[which(is.na(test_inner[,name]) == T),name] <- cvfold.inner.median.4norm[[fold]][[index.tmp]][[name]]
    }
    
    
    # Grid search on ROC_AUC
    AUC_store[[index.tmp]]<-matrix(NA,nrow = length(eta_xg), ncol = length(n_rounds_xg))
    
    for (n in 1:length(eta_xg)) {
      for (m in 1:length(n_rounds_xg)) {
        model <- xgboost(data = as.matrix(train_inner[,-which(names(train_inner) %in% c("PHENOTYPE"))]),
                         label = as.vector(train_inner$PHENOTYPE), booster = booster_xg,
                         max.depth = max_depth_xg, eta = eta_xg[n], nthread = parallelly::availableCores() - 2, nrounds = n_rounds_xg[m], eval_metric=c("logloss"),
                         objective = "binary:logistic", subsample = 0.8, colsample_bytree = 1, verbose=0, sampling_method=c("uniform"),min_child_weight =min_child_weight)
        
        Predictions <- predict(model,as.matrix(subset(test_inner,select=-c(PHENOTYPE))), type = "prob")
        observations <- test_inner[,c("PHENOTYPE")]
        roc_obj <- roc(observations, as.vector(Predictions))
        tmp.auc <-as.numeric(roc_obj$auc)
        AUC_store[[index.tmp]][n,m] <- tmp.auc
      }
    }
  }
  
  # Best parameters given ROC-AUC
  AUC_df <- matrix(NA, nrow = length(eta_xg) * length(n_rounds_xg) , ncol = (5 + cv.fold.inner))
  AUC_df <- as.data.frame(AUC_df)
  colnames(AUC_df)[1:5] <- c("eta_xg", "n_rounds_xg", "mean_AUC", "lower95CI_AUC","upper95CI_AUC")
  set = 0
  for (i in 1:length(eta_xg)) {
    for (j in 1:length(n_rounds_xg)) {
      aucs <- sapply(AUC_store,function(x){x[i,j]})
      set = set + 1
      AUC_df[set,1] <- eta_xg[i]
      AUC_df[set,2] <- n_rounds_xg[j]
      AUC_df[set,3] <- mean(aucs)
      tmp <- CI(aucs)
      AUC_df[set,4] <- as.numeric(tmp[1])
      AUC_df[set,5] <- as.numeric(tmp[3])
      AUC_df[set,6] <- aucs[1]
      AUC_df[set,7] <- aucs[2]
      AUC_df[set,8] <- aucs[3]
      AUC_df[set,9] <- aucs[4]
      AUC_df[set,10] <- aucs[5]
    }
  }
  auc_hp_index <- which.max(AUC_df[,3])
  opt_eta_xg <- AUC_df[auc_hp_index,1]
  opt_n_rounds_xg <- AUC_df[auc_hp_index,2]
  
  ## Imputation by median value
  training$PHENOTYPE<- as.character(training$PHENOTYPE)
  validation$PHENOTYPE <- as.character(validation$PHENOTYPE)
  numdata_training <- as.data.frame(sapply(training, as.numeric) )
  median_numdata_training <- sapply(numdata_training, median, na.rm = T)
  cvfold.outer.median.4norm[[fold]][[index.tmp]] <- sapply(numdata_training, median, na.rm = T)
  for (i in 1:ncol(training)) {
    name <- colnames(training)[i]
    training[which(is.na(training[,name]) == T),name] <- cvfold.outer.median.4norm[[fold]][[index.tmp]][[name]]
  }
  
  for (i in 1:ncol(validation)) {
    name <- colnames(validation)[i]
    validation[which(is.na(validation[,name]) == T),name] <- cvfold.outer.median.4norm[[fold]][[index.tmp]][[name]]
  }
  
  best_model_xg <- xgboost(data = as.matrix(training[,-which(names(training) %in% c("PHENOTYPE"))]),
                           label = as.vector(training$PHENOTYPE), booster = booster_xg,
                           max.depth = max_depth_xg, eta = opt_eta_xg, nthread = parallelly::availableCores() - 2, nrounds = opt_n_rounds_xg, eval_metric=c("logloss"),
                           objective = "binary:logistic", subsample = 0.8, colsample_bytree = 1, verbose=0, sampling_method=c("uniform"), min_child_weight =min_child_weight)
  
  
  # save model
  xgb.save(best_model_xg, fname =paste0("best_model_xg",fold))

  
  ## SHAP on test data ##
  shap_values_testset <- predict(best_model_xg, as.matrix(validation[,-which(names(validation) %in% c("PHENOTYPE"))]), predcontrib = T, approxcontrib = F)
  shap_values_testset <- as.data.table(shap_values_testset)
  shap_values_testset[,BIAS:=NULL]
  mean_shap_score_test <- colMeans(abs(shap_values_testset))[order(colMeans(abs(shap_values_testset)), decreasing = T)]
  shap_results_test <- list(shap_score = shap_values_testset, mean_shap_score_test = (mean_shap_score_test))
  
  assign(paste0("testSHAP_calc_matrix_test", cv.out),shap_values_testset)
  assign(paste0("testSHAP_calc_data_test",cv.out), as.data.table(validation[,-which(names(validation) %in% c("PHENOTYPE"))]))
  assign(paste0("testSHAP_results_test",cv.out), shap_results_test)
  rm(shap_values_testset)
  rm(shap_results_test)
  
  
  ## Store predictions ##
  Individual_pred[,3] <- predict(object = best_model_xg, newdata = as.matrix(subset(validation,select=-c(PHENOTYPE))), type = 'class')
  Individual_pred$pred_class[ Individual_pred$pred_class> 0.5 ]  <- 1
  Individual_pred$pred_class[ Individual_pred$pred_class<= 0.5 ]  <- 0
  Individual_pred[,4] <- predict(object = best_model_xg, newdata = as.matrix(subset(validation,select=-c(PHENOTYPE))), type = 'prob')
  assign(paste0("Individual_pred", cv.out),Individual_pred)
  
  TP <- c(0)
  FP <- c(0)
  FN <- c(0)
  TN <- c(0)
  for (i in 1:nrow(Individual_pred)) {
    if (Individual_pred[i,3] == "1" && Individual_pred[i,2] == "1") {
      TP = TP + 1
    } else if (Individual_pred[i,3] == "1" && Individual_pred[i,2] == "0") {
      FP = FP + 1
    } else if (Individual_pred[i,3] == "0" && Individual_pred[i,2] == "0") {
      TN = TN + 1
    } else if (Individual_pred[i,3] == "0" && Individual_pred[i,2] == "1") {
      FN = FN + 1
    }
  }
  
  # get performance measurements
  sens_a <- TP/(TP+FN)
  spec_a <- TN/(TN+FP)
  ppv_a <- TP/(TP + FP)
  npv_a <- TN/(FN+TN)
  MCC_a <- ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  F1_a <- TP /(TP+(0.5*(FP+FN)))
  AUC_a <-auc(roc(as.vector(validation[,c("PHENOTYPE")]), as.vector(predict(object = best_model_xg, newdata = as.matrix(subset(validation,select=-c(PHENOTYPE))), type = 'prob'))))
  
  # print validation results out for each fold
  outer_models_a[cv.out,1] <- opt_eta_xg
  outer_models_a[cv.out,2] <- opt_n_rounds_xg
  outer_models_a[cv.out,3] <- AUC_a
  outer_models_a[cv.out,4] <- sens_a
  outer_models_a[cv.out,5] <- spec_a
  outer_models_a[cv.out,6] <- ppv_a
  outer_models_a[cv.out,7] <- npv_a
  outer_models_a[cv.out,8] <- MCC_a
  outer_models_a[cv.out,9] <- F1_a
  rm(Individual_pred)
  
}

############### print model insights ##################

# print performance 
write.table(outer_models_a,paste0("parameters_performance_fiveouterCV.txt"), row.names = FALSE, quote = F, col.names = T, sep = "\t")

# print individual performance
MCC_list <-do.call("list",mget(grep("Individual_pred",names(.GlobalEnv),value=TRUE)))
combined_MCC <- bind_rows(MCC_list)
write.table(combined_MCC,"Individual_level_pred.txt", col.names = T, row.names = F, quote = F, sep = "\t")

# print SHAP
shap_list_testset <- do.call("list",mget(grep("testSHAP_calc_matrix_test",names(.GlobalEnv),value=TRUE)))
combined_shap_values <- bind_rows(shap_list_testset)
raw_list_testset <- do.call("list",mget(grep("testSHAP_calc_data_test",names(.GlobalEnv),value=TRUE)))
combined_raw_values <- bind_rows(raw_list_testset)
combined_shap_values <- as.data.table(combined_shap_values)
combined_raw_values <- as.data.frame(lapply(combined_raw_values, function(x) as.numeric(as.character(x))))

# Store raw data and derived SHAP values for downstream analyses
write.table(combined_raw_values,"SHAP_values_CV5_testset.txt", row.names = F, quote = F, sep = "\t")
write.table(combined_shap_values,"SHAP_raw_data_values_CV5_testset.txt", row.names = F, quote = F, sep = "\t")


# calculate and order by absolute colMeans across all cross-validations
combined_shap_values <- as.data.table(combined_shap_values)
# result from individual SHAPs averaged across all outer CVs for the five different models
mean_shap_score <- colMeans(abs(combined_shap_values))[order(colMeans(abs(combined_shap_values)), decreasing = T)]
sum(mean_shap_score)#  [1] 2.437648
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
#######################################################
