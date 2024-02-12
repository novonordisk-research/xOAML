# Script to create tables in Fig 1C, Nielsen et al, Nat Comms 2024

set.seed(432)
setwd("")
# load libraries
library(haven)
library(datasets)
library(data.table)
library(ggplot2)
library(plyr)
library(lubridate)
library(here)
library(tidyr)
library(tidyverse)
library(caret)
library(table1)

pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}


############# load data ####################
data <- fread("")

# Data most be formatted in five columns
# OA Age      BMI Sex status
# see example SourceDataFile_Fig1C_example.txt

########### Create table in Fig 1C ########
data$Sex <- as.factor(data$Sex)
data$`OA` <- as.factor(data$`OA`)
data$status <- rep("OA",nrow(data))
table1(~ Sex + Age + BMI | status * OA, data=data, overall=F, extra.col=list(`P-value`=pvalue))



