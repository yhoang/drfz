#!/usr/bin/R
### - - - - - - - - - - load R packages
add.libraries <- c("stringr")
#install.packages(load.libraries, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(add.libraries, require, character.only = TRUE)

dataset <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP", "DDPR")
feature <- c("absRange", "variance", "freq.green", "mean")
subgroup <- c("Training", "Validation", "Total")
subset <- c("full", "func", "func3", "func6")
subject <- c("df", "quadrant", "ROC", "AUC", "iAUC", 
"Error.Total", "Variables", "Heat", "SESP", "KM", 
"Summary", "Model", "VarImp", "ModelComp", "PredTotal",
"PredOOB", "ErrorList", "RMSEList")
model <- c("coxhazard", "randomforest", "boostingtree")
if (!exists("today")) today <- paste0(working.station, substring(str_replace_all(Sys.Date(), "-", ""), 3))


input_filenaming <- function(
  subgroup.var = 1,
  data = dataset[dataset.id],
  subg = subgroup,
  subs = subset[subset.id],
  subj = subject[subject.in.id],
  feat = feature[feature.id],
  com = comment.in,
  init = initdate) {
    if (data == "DDPR") {
      paste(data, sep="_")
    } else {
      paste(data, subg[subgroup.var], subs, subj, feat, com, init, sep="_")
    }
}

output_filenaming <- function(
  subject.var = 1,
  data = dataset[dataset.id],
  subs = subset[subset.id],
  subj = subject,
  feat = feature[feature.id],
  mod = model[model.id],
  com = comment.out,
  init = today) {
    if (data == "DDPR") {
      paste(data, subj[subject.var], com, mod, init, sep="_")
    }
    else {
      paste(data, subs, subj[subject.var], feat, com, mod, init, sep = "_")
    }
}