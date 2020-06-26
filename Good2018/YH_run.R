#!/usr/bin/R
# Authors: Eric Urbansky, Felix Lohrke and Yen Hoang
# DRFZ 2019 - 2020
# Goood2018

rm(list = ls())
options(max.print = 100)
### - - - - - - - - - - input
# working.station   FL, YH; different path settings
# conditional     TRUE/FALSE, if apply Ria's conditions
# dataset.id      1:5, for Basal, BCR, IL7, Pervanadate, TSLP
# subset.id       1:4 for full, func, func3, func6; 
# subject.in.id   1:13 for "df", "quadrant", "ROC", "AUC", "iAUC",  "CompTotal", "Variables", "Heat", "CompAll", "KM", "Summary"
# feature.id      1:4 for "absRange", "variance", "freq.green", "mean"
# model.id        1:4 for "cox", "RF", "boosting"
# lambda.id       1:2, 1 for lambda.min (lmin), 2 for lambda.1se (l1se)
# set.alpha       0:1, set alpha which was found to have lowest error rate, see find_alpha.R
# sampling.size   1:100, number of CV iterations
# cluster.size    1:12, find out maximum cluster size with detectCores()
working.station <- "YH"
# initdate <- "YH200410"
initdate <- "YH200526"
today <- "YH200609"
# today <- "YH200610"
# today <- "YH200623"
dataset.id <- 1
subset.id <- 4
subject.in.id <- 2
feature.id <- 1
model.id <- 2
# model.id <- 2
# comment.in <- "manSec.cof0.2"
comment.in <- "autoSec.cof0.2"
comment.out <- comment.in
if (dataset.id == 6) {
  comment.in <- comment.out <- ""
}
lambda.id <- 1
set.alpha <- 1
sampling.size <- 1
top.vars <- 20
cofactor <- 0.2
cluster.size <- 4
conditional <- FALSE
# either randomize or spikeins
randomize.label <- FALSE
# spikeIns introduces spikeIns 0 for no, 1 for 0/1 and 2 for 0/1+4 the last two columns
spikeIns <- 0

### Parameters 
# use optimal p-value via log rank test
pVal.opt = FALSE
# if above is FALSE, use manual SENS and FP values
# FP <= value
FP.value <- 0.2
# Sensitivity  >=  value
SENS.value <- 1.0
# threshold to select (should multiple apply to parameters)
selected.thresh <- "last"
# selected.thresh <- 1
### Output activated?
## generate output that doesnt change in the model (collection of ROC values, ROC curve, variables, heat map)
## generally only apply once
output.model <- FALSE
## generate output for different selected thresholds (threshold selected + info, AUC)
output.thresholds <- FALSE
lambda <- c("lmin", "l1se")


if (FALSE) {
  initdate <- "YH200526"
  comment.in <- "autoSec.cof0.2"
  model.id <- 2
  dataset.id <- 1
  subset.id <- 4
  randomize.label <- FALSE
  sampling.size <- 1
}


if (model.id == 1 | model.id == 3) {
  if (pVal.opt) {
    comment.out <- sprintf("%s.pValOpt", comment.in, FP.value, SENS.value, selected.thresh)
  } else {
    comment.out <- sprintf("%s.FP%s.SENS%s.Thresh%s", comment.in, FP.value, SENS.value, selected.thresh)
  }
  if (model.id == 1) {
    comment.out <- paste(comment.out, lambda[lambda.id], sep=".")
  }
}

if (randomize.label) {
  comment.out <- paste0(comment.out, ".rand")
}
if (spikeIns > 0) {
  comment.out <- paste0(comment.out, ".spikes", spikeIns)
}
if (randomize.label && spikeIns > 0) stop("Use either randomized labels or SpikeIns in a run!")

#############################################################






### - - - - - - - - - - load R packages
load.libraries <- c("survival", "readxl", "dplyr", "tidyverse", "xlsx", "survminer", "pROC", "survAUC")
#install.packages(load.libraries, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(load.libraries, require, character.only = TRUE)
#a useful print function
printf <- function(...) invisible(print(sprintf(...)))
### load custom R script to generate standardized file names
source(file.path("D:", "drfz", "Good2018", "YH_func_naming.R"))

### set path automatically
if (working.station == "FL") {
    Rdata.path <- file.path("", "home", "felix", "AG_Baumgrass", "Scripts", "Pri", "Pri_good_established", "Rds")
    Text.path <- file.path("", "home", "felix", "AG_Baumgrass", "Data", "Good", "Marker_combinations")
    Project.path <- file.path("", "home", "felix", "AG_Baumgrass", "Data", "Good", "Basal_for_Felix", "Basal_for_Felix")
    Output.path <- file.path("", "home", "felix", "AG_Baumgrass", "Results", "Good_BCR", "NoCond")
} else if (working.station == "YH") {
  Project.path <- file.path("D:", "drfz", "Good2018")
  Rdata.path <- file.path("D:", "drfz", "Good2018", "Rdata", dataset[dataset.id])
  Text.path <- Cohort.path <- file.path("D:", "drfz", "Good2018", "tables")
  if (dataset.id == 6) {
    Output.path <- file.path(Project.path, model[model.id], dataset[dataset.id])
  } else {
    Output.path <- file.path(Project.path, model[model.id], dataset[dataset.id], subset[subset.id])
  }
    
}

######## INPUT FILES #####################
# subgroup <- c("Training", "Validation", "Total")
training.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(1))
validation.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(2))
# total data as rds (if it does not exist a new will one be created at path location)
total.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(3))
##########################


######## INITIATE OUTPUT FILES #####################
### create folder first
dir.create(Output.path, recursive = TRUE, showWarnings = FALSE)
### subject = c("df", "quadrant", "ROC", "AUC", "iAUC", "Comp.Val", "Comp.Train", "Comp.Total", "Variables", "Heat", "Comp.Error", "KM", "Summary")
# ROC / AUC Curve pdf
roc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(3))
auc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(4))
iauc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(5))
# comparisons of validation, training and total xlsx
comp.total.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(6))
#### variables xlsx
variables.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(7))
# heatmap pdf
heat.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(8))
# comparisons of Sensitivity, specificity and FP-rate to associated thresholds
comp.SESPFP.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(9))
# kaplan-meier curve
kaplan.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(10))
# summary txt
summary.path <- sprintf("%s/%s.log", Output.path, output_filenaming(11))
# tree pdf
tree.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(12))
# VarImportance PDF
varimp.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(13))
################################


### set seed for reproduction 
seed.vec <- sample(100:(100 + sampling.size))
### ------------------- load PRI features and patient meta data and preprocess
# patient metadata 
patient.data.path <- file.path(Text.path, "patient_cohort.xlsx")
### load custom R script to load patient metadata, which covers in detail
# loading RDS files
# create df.total and save
# load patient metadata and convert into numerical data
# add column "Survival Time"
# convert data frames as matrices
# convert NaNs and NAs
source(file.path("D:", "drfz", "Good2018", "YH_PRIfeature_preprocession.R"))
source(file.path("D:", "drfz", "Good2018", "YH_mem.alloc.R"))



if (FALSE) {

### ---------------  run model
if (model.id == 1) {
  source(file.path("D:", "drfz", "Good2018", "coxhazard", "YH_model.coxhazard.R"))
} else if (model.id == 2) {
  source(file.path("D:", "drfz", "Good2018", "randomforest", "YH_model.randomforest.R"))
} else if (model.id == 3) {
  source(file.path("D:", "drfz", "Good2018", "boostingtree", "YH_model.boostingtree.R"))
} else if (model.id == 4) {
  source(file.path("D:", "drfz", "Good2018", "YH_rf-cox.model.R"))
} else {
  print("Chose model.id between 1:4 for cox, random forest, boosting tree or random forest-cox.")
}
}