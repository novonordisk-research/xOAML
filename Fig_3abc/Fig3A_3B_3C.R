setwd("")
getwd()


################## Load libraries ###############
library(pROC)
library(caret)
library(yardstick)
library(Rmisc)


########### Load data ###############
MCC_subset <- fread("/Figure3A3B.txt", header = T)
MCC_subset$obs <- as.factor(MCC_subset$obs)

######### fig 3A ########################
MCC_arm <- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Arm"))
MCC_arm <- MCC_arm %>% 
  unique()


MCC_foot  <- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Foot"))
MCC_foot <- MCC_foot %>% 
  unique()

MCC_hip <- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Hip"))
MCC_hip <- MCC_hip %>% 
  unique()

MCC_knee <- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Knee"))
MCC_knee <- MCC_knee %>% 
  unique()

MCC_spine<- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Spine"))
MCC_spine <- MCC_spine %>% 
  unique()


MCC_unOA <- subset(MCC_subset, is.na(MCC_subset$OAdef) | MCC_subset$OAdef == c("Unspecified OA"))
MCC_unOA <- MCC_unOA %>% 
  unique()


combined <- rbind(MCC_arm,MCC_foot,MCC_hip,MCC_knee,MCC_spine,MCC_unOA)

MCC_pr <- roc_curve(MCC, obs, pred_prob,event_level = "second") 
mccroc <- roc(MCC$obs, MCC$pred_prob)
MCC_pr$Model <- rep(paste0("OA all (",round(mccroc$auc[1],digits = 2),")"), nrow(MCC_pr))


MCC_pr_arm <- roc_curve(MCC_arm, obs, pred_prob,event_level = "second") 
mccroc_arm <- roc(MCC_arm$obs, MCC_arm$pred_prob)
MCC_pr_arm$Model <- rep(paste0("Arm (",round(mccroc_arm$auc[1],digits = 2),")"), nrow(MCC_pr_arm))

MCC_pr_foot <- roc_curve(MCC_foot, obs, pred_prob,event_level = "second") 
mccroc_foot <- roc(MCC_foot$obs, MCC_foot$pred_prob)
MCC_pr_foot$Model <- rep(paste0("Foot (",round(mccroc_foot$auc[1],digits = 2),")"), nrow(MCC_pr_foot))


MCC_pr_hip <- roc_curve(MCC_hip, obs, pred_prob,event_level = "second") 
mccroc_hip <- roc(MCC_hip$obs, MCC_hip$pred_prob)
MCC_pr_hip$Model <- rep(paste0("Hip (",round(mccroc_hip$auc[1],digits = 2),")"), nrow(MCC_pr_hip))


MCC_pr_knee <- roc_curve(MCC_knee, obs, pred_prob,event_level = "second") 
mccroc_knee <- roc(MCC_knee$obs, MCC_knee$pred_prob)
MCC_pr_knee$Model <- rep(paste0("Knee (",round(mccroc_knee$auc[1],digits = 2),")"), nrow(MCC_pr_knee))

MCC_pr_spine <- roc_curve(MCC_spine, obs, pred_prob,event_level = "second") 
mccroc_spine <- roc(MCC_spine$obs, MCC_spine$pred_prob)
MCC_pr_spine$Model <- rep(paste0("Spine (",round(mccroc_spine$auc[1],digits = 2),")"), nrow(MCC_pr_spine))

mccroc <- roc(MCC$obs, MCC$pred_prob)
mcc2roc <- roc(MCC2$obs, MCC2$pred_prob)


MCC_combined <- rbind(MCC_pr,MCC_pr_arm,
                      MCC_pr_foot,MCC_pr_hip,
                      MCC_pr_knee,MCC_pr_spine) 
table(MCC_combined$Model)



pp <- ggplot(MCC_combined, aes(x=(1-specificity), y=sensitivity,colour=Model)) + 
  geom_line(size=1.5,alpha = 0.6) +
  theme_minimal() +
  xlab('1 - Specificity') +
  ylab('Sensitivity') + geom_abline(intercept=0, slope=1, color="grey", size = 1, alpha=0.6,linetype="dashed") +
  theme(axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20,face="bold"),
        axis.text.x=element_text(size=20, angle = 45, hjust=1),
        axis.title.x=element_text(size=20,face="bold"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))
pp
png(filename = "ROC_curve_oneModel.png", width = 650, height = 650)
print(pp)
dev.off()




########### Fig 3B ######
MCC <- MCC_subset
MCC$obs <- as.factor(MCC$obs )
MCC_pr <- pr_curve(MCC, obs, pred_prob,event_level = "second") 
mccroc <- roc(MCC$obs, MCC$pred_prob)
MCC_pr$Model <- rep(paste0("OA all"), nrow(MCC_pr))


MCC_pr_arm <- pr_curve(MCC_arm, obs, pred_prob,event_level = "second")
mccroc_arm <- roc(MCC_arm$obs, MCC_arm$pred_prob)
MCC_pr_arm$Model <- rep(paste0("Arm"), nrow(MCC_pr_arm))

MCC_pr_foot <- pr_curve(MCC_foot, obs, pred_prob,event_level = "second") 
mccroc_foot <- roc(MCC_foot$obs, MCC_foot$pred_prob)
MCC_pr_foot$Model <- rep(paste0("Foot"), nrow(MCC_pr_foot))


MCC_pr_hip <- pr_curve(MCC_hip, obs, pred_prob,event_level = "second") 
mccroc_hip <- roc(MCC_hip$obs, MCC_hip$pred_prob)
MCC_pr_hip$Model <- rep(paste0("Hip"), nrow(MCC_pr_hip))


MCC_pr_knee <- pr_curve(MCC_knee, obs, pred_prob,event_level = "second") 
mccroc_knee <- roc(MCC_knee$obs, MCC_knee$pred_prob)
MCC_pr_knee$Model <- rep(paste0("Knee"), nrow(MCC_pr_knee))

MCC_pr_spine <- pr_curve(MCC_spine, obs, pred_prob,event_level = "second") 
mccroc_spine <- roc(MCC_spine$obs, MCC_spine$pred_prob)
MCC_pr_spine$Model <- rep(paste0("Spine"), nrow(MCC_pr_spine))

MCC_pr_unOA <- pr_curve(MCC_unOA, obs, pred_prob,event_level = "second") 
mccroc_unOA <- roc(MCC_unOA$obs, MCC_unOA$pred_prob)
MCC_pr_unOA$Model <- rep(paste0("Unspecified OA"), nrow(MCC_pr_unOA))



MCC_combined <- rbind(MCC_pr,MCC_pr_arm,
                      MCC_pr_foot,MCC_pr_hip,
                      MCC_pr_knee,MCC_pr_spine) 




ppp <- ggplot(MCC_combined, aes(x=recall, y=precision,colour=Model)) +
  geom_line(size=1.5,alpha = 0.6) +
  theme_minimal() +
  xlab('Sensitivity') +
  ylab('PPV') +
  theme(axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20,face="bold"),
        axis.text.x=element_text(size=20, angle = 45, hjust=1),
        axis.title.x=element_text(size=20,face="bold"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20))
ppp


png(filename = "PR_curve.png")
print(ppp)
dev.off()




#### Figure 3C #######
CV_perf <- read.table("/SourceDataFile_parameters_performance_fiveouterCV.txt", header = T)
as.data.frame(apply(CV_perf, 2, CI))







