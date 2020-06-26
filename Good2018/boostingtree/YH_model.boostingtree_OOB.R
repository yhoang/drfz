#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

# in command
#cd C:Program\ Files/R/R-3.6.1/bin
#./R.exe --max-ppsize 500000
# or short:
# "C:Program Files\R\R-3.6.1\bin\R.exe" --max-ppsize 500000
# source(file.path("D:", "drfz", "Good2018", "YH_run.R"))
options(expressions = 5e5)
# memory.limit(size = 5000000)

# install.packages("gbm")
add.libraries <- c("gbm", "doParallel", "ggplot2")
lapply(add.libraries, require, character.only = TRUE)
### no need at OOB
if (exists("df.training")) rm(df.training)
if (exists("df.validation")) rm(df.validation)

print("### Running bossting tree model with OOB on df.total..")
### ------------------ boosting tree model 
it.total <- 0
it.total <- it.total + 1 
set.seed(seed.vec[it.total])

timeStart <- Sys.time()
printf("### Start %s with %s on dataset=%s comment.in=%s subset=%s random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], comment.in, subset[subset.id], randomize.label, spikeIns)
  
### fits generalized boosted regression models
# n.trees             number of gardient boosting interation, default = 100
# shrinkage           learning rate, default = 0.001
# bag.fraction        subsampling frac to select the next tree in the expansion, default = 0.5
# interaction.depth   number of splits it has to perform on a tree, default = 1
# n.minobsinnode      the minimum number of observations in trees' terminal nodes, default = 10
# As each split increases the total number of nodes by 3 and number of terminal nodes by 2, the total number of nodes in the tree will be 3∗N+1 and the number of terminal nodes 2∗N+1
ntree.val <- 100
shrink.val <- 0.001

cluster.size <- 2
cl <- makeCluster(cluster.size)
registerDoParallel(cl)
ptm <- proc.time()
for (oob in 1:nrow(df.total)) {
  printf("Run model=%s with n.trees=%s, shrinkage=%s and oob=%s..", model[model.id], ntree.val, shrink.val, oob)

  comment.tmp <- paste(comment.in, paste0("ntree", ntree.val), paste0("shrink", shrink.val), paste0("oob", oob), sep=".")
  model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, com=comment.tmp))

  if (!file.exists(model.path)) {
    gbm.cox <- vector()
    gbm.cox <- gbm(formula = Surv(Survivaltime, RelapseStatus) ~ . - DDPRStatus - PatientID,
                  data = as.data.frame(df.total[-oob, ]),
                  distribution = "coxph",
                  n.trees = ntree.val,
                  shrinkage = shrink.val, 
                  verbose = FALSE,
                  n.cores = cluster.size
                  )
    saveRDS(gbm.cox, model.path)
    # free up unused RAM
    rm(gbm.cox)
    gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)  #This will 
  }
}

printf("Done with creating model=%s with n.trees=%s, shrinkage=%s.", model[model.id], ntree.val, shrink.val)
print(proc.time() - ptm)
stopCluster(cl)
