#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018
# install.packages("gbm")
add.libraries <- c("gbm", "doParallel")
lapply(add.libraries, require, character.only = TRUE)

print("Running bossting tree model..")

### ------------------ boosting tree model 
it.total <- 0
all.activ.Index <- all.coef.value <- vector()
variabl.count <- vector()
min.cvm <- 100

timeStart <- Sys.time()
ptm <- proc.time()
printf("### Start %s with %s on dataset=%s comment.in=%s subset=%s  random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], comment.in, subset[subset.id], randomize.label, spikeIns)

while (it.total < sampling.size) {
  it.total <- it.total + 1 
  set.seed(seed.vec[it.total])
  
  # in command
  #cd C:Program\ Files/R/R-3.6.1/bin
  #./R.exe --max-ppsize 500000
  ### fits generalized boosted regression models
  # n.trees             number of gardient boosting interation, default = 100
  # shrinkage           learning rate, default = 0.001
  # bag.fraction        subsampling frac to select the next tree in the expansion, default = 0.5
  # interaction.depth   number of splits it has to perform on a tree, default = 1
  # n.minobsinnode      the minimum number of observations in trees' terminal nodes, default = 10

  # As each split increases the total number of nodes by 3 and number of terminal nodes by 2, the total number of nodes in the tree will be 3∗N+1 and the number of terminal nodes 2∗N+1
cluster.size <- 3
cl <- makeCluster(cluster.size)
registerDoParallel(cl)
ptm <- proc.time()
  gbm.cox <- gbm(formula = Surv(Survivaltime, RelapseStatus) ~ . - DDPRStatus,
    data = as.data.frame(df.training),
    distribution = "coxph",
    n.trees = 500,
    cv.folds = 4,
    shrinkage = 0.01, 
    verbose = FALSE,
    n.cores = cluster.size
  )
print(proc.time() - ptm)
stopCluster(cl)

# estimates the optimal number of boosting iterations 
best.iter <- gbm.perf(gbm.cox, method="OOB", plot.it = TRUE, oobag.curve = TRUE)  # returns out-of-bag estimated best number of trees
print(best.iter)
# best.iter <- gbm.perf(gbm.cox, method="cv") # returns test set estimate of best number of trees
# print(best.iter)
# best.iter <- gbm.perf(gbm.cox, method="test") # returns test set estimate of best number of trees
# print(best.iter)

all.coef.value <- summary(gbm.cox, n.trees = best.iter, plotit = FALSE) # using estimated best number of trees
top.coef.value <- all.coef.value[1:top.vars, ]

# Plot relative influence of each variable
if (output.model) {
  pdf(varimp.path, width = 12)
  # par(mfrow = c(1, 2), oma = 4, 1, 1, 1)
  par(mfrow = c(1, 2), oma = 4, 1, 1, 1, mar=c(3, 3, 4, 2))
  # summary(gbm.cox, n.trees = 1, cBars = 20) # using first tree
  summary(gbm.cox, n.trees = best.iter, cBars = 20) # using estimated best number of trees
  dev.off()
}

  # Compactly print the first and last trees for curiosity
  print(pretty.gbm.tree(gbm.cox, i.tree = 1))
  print(pretty.gbm.tree(gbm.cox, i.tree = gbm.cox$n.trees))

pred.training <- predict(gbm.cox, newdat=as.data.frame(df.training[, - c(1:3)]), n.trees=100)
pred.validation <- predict(gbm.cox, newdat=as.data.frame(df.validation[, - c(1:3)]), n.trees=100)


### --------------- calculat threshold for cutoff #####################
# "low Risk" < treshold > "high risk"
### scale RR between 0 -> 1
pred.training.shift = pred.training - min(pred.training)
scale.pred.training <- pred.training.shift / (max(pred.training) - min(pred.training))

#performance test for set of thresholds
threshold <- seq(0, 1, 0.01)
predictions.roc <- data.frame()
#safe all prediction cutoffs as df  
for (it in 1:length(threshold)){
  newline <- findInterval(scale.pred.training, threshold[it])
  predictions.roc <- bind_cols(as.data.frame(newline), predictions.roc)
  names(predictions.roc)[1]= threshold[it]
}
#mirror df to set index right 
predictions.roc= predictions.roc[, ncol(predictions.roc):1]

###calculate error
  #true prositive = prediction 1 & relapse status 1
  #calculate ssensitivity and false positiv rate 
  all.spec <- all.sens <- all.FP.rate <- AUC <- vector()
  for (i in 1:length(threshold)){
    TP <- length(which(predictions.roc[, i] == 1 & df.training[, 2] == 1))
    TN <- length(which(predictions.roc[, i] == 0 & df.training[, 2] == 0))
    FP <- length(which(predictions.roc[, i] == 1 & df.training[, 2] == 0))
    ALLP <- length(which(df.training[, 2] == 1))
    ALLN <- length(which(df.training[, 2] == 0))
    SENS <- TP / ALLP
    SPEC <- TN / ALLN
    all.spec <- c(all.spec, SPEC)
    all.sens <- c(all.sens, SENS)
    FP.rate <- FP / ALLN
    all.FP.rate <- c(all.FP.rate, FP.rate)
    AUC <- c(AUC, auc(roc(df.training[, 2], predictions.roc[, i])))
  }

# save  ROC plot as pdf
  if (output.model){
    pdf(roc.pdf.path)
    par(mgp = c(2, 0.5, 0), pty = "s", mar = c(4, 4, 3, 1), cex.lab=1.2)
    plot(all.FP.rate, all.sens, type = "l", ylab = "TP Rate (Sensitivity)", 
    xlab = "FP Rate", ylim = c(0, 1), xlim = c(0, 1), main = "ROC Curve")
    dev.off()
  }
# saving all FP - rates, SENS, SPEC and associated thresholds
  if (output.model){
    comp.all <- data.frame(data.frame(threshold * max(pred.training)), data.frame(all.sens), data.frame(all.spec), data.frame(all.FP.rate))
    names(comp.all)[1] <- "threshold"
    names(comp.all)[2] <- "SENS"
    names(comp.all)[3] <- "SPEC"
    names(comp.all)[4] <- "FP-rate"
    write.xlsx(comp.all, comp.SESPFP.path)
  }

if (pVal.opt) {
    ###########################################################################
    # #treshold by log - rang test. find "optimal" p - value for log rank
    p.value <- c()
    for (i in 1:ncol(predictions.roc)) {
      df.new <- data.frame(df.training[, 1], df.training[, 2], predictions.roc[, i])
      
      names(df.new)[1] <- "Time"
      names(df.new)[2] <- "Status"
      names(df.new)[3] <- "Prediction"
      # calculate estimate of survival curve for censored data using Kaplan-Meier 
      km.type = survfit(
        Surv(df.new$Time, df.new$Status) ~ df.new$Prediction,
        data = df.new,
        type = "kaplan-meier")
      # or try coxph( Surv(df.new$Time, df.new$Status) ~ df.new$Prediction, data = df.new)
      # calculate log rank test and extract p-value
      tmp <- surv_pvalue(km.type, method = "1")$pval
      p.value <- c(p.value, tmp)
    }
    p.value[which(is.na(p.value) == TRUE)] <- 1
    op.thresh <- threshold[which(p.value == min(p.value))[1]]
    op.thresholds <- op.thresh
    ##########################################################################
  } else {
    ### Eric's interpretation BUT DO NOT USE FURTHER ON
    #finde best threshold
    # temp <- which(all.sens[which(all.FP.rate == FP.value)] <= SENS.value)[1]
    # op.thresh <- which(all.sens == all.sens[which(all.FP.rate == FP.value)][temp])
    # op.thresholds <- op.thresh

    ### SENS = 1 and FP lowest
    # does not work here
    sens.idx <- which(all.sens >= SENS.value)
    ## get lowest FP
    FP.min.idx <- which(all.FP.rate == min(all.FP.rate[sens.idx]))
    op.thresholds <- FP.min.idx
    printf("Positions that have been selected as possible thresholds: %s", op.thresholds)

    ## selected threshold via position selected.thresh
    op.thresh <- threshold[op.thresholds][selected.thresh]
  }
  printf("Varied threshold Values: %s", op.thresh)
  printf("RR over %s are interpret as Status 1 (High Risk) and under as Status 0(Low Risk)", op.thresh * max(pred.training))


  ## saving threshold data in txt
  if (output.thresholds){
    sink(file = summary.path, append = TRUE)
    print(timeStart)
    printf("Sampling size: %s", sampling.size)
    print("Variable distribution in models:")
    print(table(variabl.count))
    printf("Best model >0: %s features", op.variabl.count)
    printf("Random sampling of condition: %s", randomize.label)
    printf("Introduce Spike Ins: %s", spikeIns)
    printf("optimal p-Value search via log-rank: %s", pVal.opt)
    printf("Selected Parameters for FP rate <= %s and for SENS  >=  %s", FP.value, SENS.value)
    printf("Possible thresholds for this configuration: %s", op.thresholds)
    printf("Specifically selected Threshold for this configuration: %s at position: %s", op.thresholds[selected.thresh], selected.thresh)
    printf("RR for this threshold: %s", op.thresh * max(pred.training))
    sink()
  }
  
  #pediction RR Values aprox Relapsstaus with calculates Threshold
  predicted.validation=findInterval(pred.validation / max(pred.training), op.thresh)
  #print(predicted.validation)
  predicted.training=findInterval(pred.training / max(pred.training), op.thresh)
  #print(predicted.training)

MSE.training <- mean((scale.pred.training - df.training[, 2]) ^ 2)
MSE.validation <- mean((scale.pred.validation - df.validation[, 2]) ^ 2)
MSE.total <- mean((c(scale.pred.training, scale.pred.validation) - c(df.training[, 2], df.validation[, 2])) ^ 2)
}