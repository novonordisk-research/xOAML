# Code to reproduce Figure 1B, UpSet plot of OA diagnoses
#Nielsen et al, Nat Comms 2024

# PLease note due to individual level data sensitive data, real source data has not been provided, just examples to reproduce the plot. 
# some eid in the real data is duplicated due to mulitple diagnosis. Source data eid have been named by the row number and hence not containing duplicated eids representing people with multiple diagnoses

dir <- c("")
setwd(dir)


################# Load libraries ###############
library(ComplexHeatmap)

########### Load data ########################

OA_diagnoses <- read.table("/Source_Data_Fig1B_OAdiagnoses.txt", header = T, sep = "\t")
Perf <- read.table("Source_Data_Fig1B_Individual_level_pred.txt", header = T)
Perf1 <- merge(Perf,OA_diagnoses, by = c("eid"), all.x=T)

# subset to OA-diagnosed patients 
Perf1_cases <- subset(Perf1,Perf1$obs == 1)

# indicate which diagnoses a patient had from primary+secondary EHR 
Perf1_cases_enc <-
  Perf1_cases %>% group_by(eid) %>%
  mutate(Arm=if_else(OAdef == "Arm", true = 1, false = 0)) %>%
  mutate(Foot=if_else(OAdef == "Foot", true = 1, false = 0)) %>%
  mutate(Hip=if_else(OAdef == "Hip", true = 1, false = 0)) %>%
  mutate(Knee=if_else(OAdef == "Knee", true = 1, false = 0)) %>%
  mutate(Spine=if_else(OAdef == "Spine", true = 1, false = 0)) %>%
  mutate(Unspecified_OA=if_else(OAdef == "Unspecified OA", true = 1, false = 0))

# Upset plot
hh<-split(Perf1_cases_enc$eid, Perf1_cases_enc$OAdef)
n_sets <-length(split(Perf1_cases_enc$eid, Perf1_cases_enc$OAdef))
m3 <- make_comb_mat(hh, mode = "intersect") 
UpsetPlot <- UpSet((m3),
      comb_order=order(comb_size(m3),decreasing = T),
      top_annotation = upset_top_annotation((m3), add_numbers = TRUE), 
      right_annotation = upset_right_annotation((m3), add_numbers = TRUE) )


print(UpsetPlot)

