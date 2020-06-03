#!/usr/bin/R

dataset <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")
feature <- c("absRange", "variance", "freq.green", "mean")
subgroup <- c("Training", "Validation", "Total")
subset <- c("full", "func", "func3", "func6")
subject <- c("df", "quadrant", "ROC", "AUC", "iAUC", "Comp.Total", "Variables", "Heat", "Comp.Error", "KM", "Summary", "Tree", "VarImp")
model <- c("coxhazard", "randomforest", "boostingtree")
today <- paste0(working.station, substring(str_replace_all(Sys.Date(), "-", ""), 3))


input_filenaming <- function(
  subgroup.var = 1,
  data = dataset[dataset.id],
  subg = subgroup,
  subs = subset[subset.id],
  subj = subject[subject.in.id],
  feat = feature[feature.id],
  com = comment.in,
  init = initdate) {
    paste(data, subg[subgroup.var], subs, subj, feat,  com, init, sep="_")
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

    paste(data, subs, subj[subject.var], feat, com, mod, init, sep = "_")
}