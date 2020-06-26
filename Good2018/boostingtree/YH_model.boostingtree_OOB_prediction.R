#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

# in command
# cd C:Program\ Files/R/R-3.6.1/bin
# ./R.exe --max-ppsize 500000
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

### ------------------ boosting tree model 
it.total <- 0
it.total <- it.total + 1 
set.seed(seed.vec[it.total])

timeStart <- Sys.time()
printf("### Start %s with %s on dataset=%s comment.in=%s subset=%s selected.thresh=%s random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], comment.in, subset[subset.id], selected.thresh, randomize.label, spikeIns)
  
### fits generalized boosted regression models
# n.trees             number of gardient boosting interation, default = 100
# shrinkage           learning rate, default = 0.001
# bag.fraction        subsampling frac to select the next tree in the expansion, default = 0.5
# interaction.depth   number of splits it has to perform on a tree, default = 1
# n.minobsinnode      the minimum number of observations in trees' terminal nodes, default = 10
# As each split increases the total number of nodes by 3 and number of terminal nodes by 2, the total number of nodes in the tree will be 3∗N+1 and the number of terminal nodes 2∗N+1
ntree.val <- 100
shrink.val <- 0.001
status.idx <- which(colnames(df.total) == "RelapseStatus")
ddpr.idx <- which(colnames(df.total) == "DDPRStatus")
time.idx <- which(colnames(df.total) == "Survivaltime")
sample.idx <- which(colnames(df.total) == "PatientID")

cluster.size <- 2
cl <- makeCluster(cluster.size)
registerDoParallel(cl)

### ---------------- run for prediction -------------------------------- ###
global.error.pred <- global.error.oob <- global.RMSE.training <- global.RMSE.OOB <- global.RMSE.total <- vector()

 # error and RMSE list
comment.model.pred <- paste(comment.out, paste0("ntree", ntree.val), paste0("shrink", shrink.val), sep=".")
errorlist.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(17, 
    com=comment.model.pred))
RMSElist.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(18, 
    com=comment.model.pred))

for (oob in 1:nrow(df.total)) {
# for (oob in 1:1) {
  local.thresh <- selected.thresh
  ###---------------- initiate output file names --------------------------###
  comment.tmp <- paste(comment.in, paste0("ntree", ntree.val), paste0("shrink", shrink.val), paste0("oob", oob), sep=".")
  comment.model.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("shrink", shrink.val), paste0("oob", oob), sep=".")
  roc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(3, com =  
    comment.model.tmp))
  #### variables xlsx
  variables.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(7, com 
    = comment.model.tmp))
  # comparisons of Sensitivity, specificity and FP-rate to associated thresholds
  comp.SESPFP.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(9, 
    com = comment.model.tmp))
  kaplan.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(10, 
    com=comment.model.tmp))
  summary.path <- sprintf("%s/%s.log", Output.path, output_filenaming(11, com = 
    comment.model.tmp))
  # load model rds
  model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, 
    com=comment.tmp))
  # VarImportance PDF
  varimp.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(13, com = 
    comment.model.tmp))
  # prediction rds
  pred.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(15, 
    com=comment.tmp))
  pred.oob.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(16, 
    com=comment.tmp))
  ###########################################################################

  ### load model
  if (TRUE) {
    printf("Loading model %s..", model[model.id])
    gbm.cox <- readRDS(model.path)
  
    # estimates the optimal number of boosting iterations 
    best.iter <- gbm.perf(gbm.cox, method="OOB", plot.it = FALSE, oobag.curve = TRUE)  # returns out-of-bag estimated best number of trees
    # print(best.iter)

    # using estimated best number of trees
    all.coef.value <- summary(gbm.cox, n.trees = best.iter, plotit = FALSE) 
    top.coef.value <- all.coef.value[which(all.coef.value[, 2] > 0), ]
    # reduce top variables to 20
    if (length(top.coef.value) > 20) top.coef.value <- top.coef.value[1:top.vars, ]
  }

  ### predict GBM on Training
  if (!file.exists(pred.path)) {
    pred.total <- predict(gbm.cox, 
      newdat = as.data.frame(df.total[-oob, -c(status.idx, ddpr.idx, time.idx, sample.idx)]), n.trees = best.iter)
    saveRDS(pred.total, pred.path)
    printf("pred.total saved as %s.", pred.path)
  } else {
    printf("File exists already: %s", pred.path)
    pred.total <- readRDS(pred.path)
  }

  ### predict GBM on OOB
  if (!file.exists(pred.oob.path)) {
    pred.oob <- predict(gbm.cox, 
      newdat = as.data.frame(df.total[oob, -c(status.idx, ddpr.idx, time.idx, sample.idx)]), n.trees = best.iter)
    saveRDS(pred.oob, pred.oob.path)
    printf("pred.oob saved as %s.", pred.oob.path)
  } else {
    printf("File exists already: %s", pred.oob.path)
    pred.oob <- readRDS(pred.oob.path)
  }

  ### --------------- calculat threshold for cutoff #####################
  # "low Risk" < treshold > "high risk"
  ### scale RR between 0 -> 1
  pred.total.shift <- pred.total - min(pred.total)
  scale.pred.total <- pred.total.shift / diff(range(pred.total))
  scale.pred.oob <- pred.oob / diff(range(pred.total))

  #performance test for set of thresholds
  threshold <- seq(0, 1, 0.01)
  predictions.roc <- data.frame()
  #safe all prediction cutoffs as df  
  for (it in 1:length(threshold)){
    newline <- findInterval(scale.pred.total, threshold[it])
    predictions.roc <- bind_cols(as.data.frame(newline), predictions.roc)
    names(predictions.roc)[1] <- threshold[it]
  }
  #mirror df to set index right 
  predictions.roc <- predictions.roc[, ncol(predictions.roc):1]

  ### calculate error
  #true prositive = prediction 1 & relapse status 1
  #calculate ssensitivity and false positiv rate 
  all.spec <- all.sens <- all.FP.rate <- AUC <- vector()
  for (i in 1:length(threshold)){
    TP <- length(which(predictions.roc[, i] == 1 & df.total[-oob, status.idx] == 1))
    TN <- length(which(predictions.roc[, i] == 0 & df.total[-oob, status.idx] == 0))
    FP <- length(which(predictions.roc[, i] == 1 & df.total[-oob, status.idx] == 0))
    FN <- length(which(predictions.roc[, i] == 0 & df.total[-oob, status.idx] == 1))
    ALLP <- length(which(df.total[-oob, status.idx] == 1))
    ALLN <- length(which(df.total[-oob, status.idx] == 0))
    SENS <- TP / ALLP
    SPEC <- TN / ALLN
    all.spec <- c(all.spec, SPEC)
    all.sens <- c(all.sens, SENS)
    FP.rate <- FP / ALLN
    all.FP.rate <- c(all.FP.rate, FP.rate)
    AUC <- c(AUC, auc(roc(unlist(df.total[-oob, status.idx]), predictions.roc[, i], quiet=TRUE)))
  }


  thresh.alert <- FALSE
  if (pVal.opt) {
    ###########################################################################
    # #treshold by log - rang test. find "optimal" p - value for log rank
    p.value <- c()
    for (i in 1:ncol(predictions.roc)) {
      df.new <- data.frame(
        Survivaltime = df.total[-oob, time.idx], 
        RelapseStatus = df.total[-oob, status.idx], 
        Prediction = predictions.roc[, i])
      
      # calculate estimate of survival curve for censored data using Kaplan-Meier 
      km.type <- survfit(
        Surv(Survivaltime, RelapseStatus) ~ Prediction,
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
    ## get lowest FP rate
    FP.idx <- which(all.FP.rate <= FP.value)
    ### SENS = 1 and FP lowest
    sens.idx <- which(all.sens == max(all.sens[FP.idx]))
    if (length(FP.idx) == 0 | length(sens.idx) == 0){
      thresh.alert <- TRUE
      # get hightest AUC rate, first position
      sens.idx <- which(AUC == max(AUC))[1]
      FP.idx <-  which(AUC == max(AUC))[1]
    } 
    op.thresholds <- intersect(FP.idx, sens.idx)
    if (thresh.alert) printf("!!!! ALERT !!!! No threshold found for above FP/SENSE setting. Take first position of highest AUC!")
    printf("Positions that have been selected as possible thresholds: %s", paste(op.thresholds, collapse = " "))
    ## selected threshold via position local.thresh
    if (selected.thresh == "last") {
      local.thresh <- length(op.thresholds)
    }
    op.thresh <- threshold[op.thresholds]
    op.thresh <- op.thresh[local.thresh]
  }
  printf("Varied threshold Values: %s", op.thresh)
  printf("RR over %0.3f are interpret as Status 1 (High Risk) and under as Status 0(Low Risk)", op.thresh * max(pred.total))

  # prediction RR Values aprox Relapsstatus with calculates Threshold
  pred.total.class <- findInterval(scale.pred.total, op.thresh)
  pred.oob.class <- findInterval(scale.pred.oob, op.thresh)
  CM.total <- table(pred.total.class, unlist(df.total[-oob, status.idx]))
  # confusion matrix
  print(CM.total)

  ### collect error count
  error.pred <- sum(CM.total) - sum(diag(CM.total))
  error.oob <- unlist(df.total[oob, status.idx]) - as.numeric(as.character(pred.oob.class))

  ### collect RMSE
  RMSE.training <- sqrt(mean((scale.pred.total - unlist(df.total[-oob, status.idx])) ^ 2))
  RMSE.OOB <- sqrt(mean((scale.pred.oob - unlist(df.total[oob, status.idx])) ^ 2))
  RMSE.total <- sqrt(mean((c(scale.pred.total, scale.pred.oob) - c(unlist(df.total[-oob, status.idx]), unlist(df.total[oob, status.idx]))) ^ 2))

  ### collect for every run
  global.error.pred <- c(global.error.pred, error.pred)
  global.error.oob <- c(global.error.oob, error.oob)
  global.RMSE.training <- c(global.RMSE.training, RMSE.training)
  global.RMSE.OOB <- c(global.RMSE.OOB, RMSE.OOB)
  global.RMSE.total <- c(global.RMSE.total, RMSE.total)

  ### ------------------- saving summary as txt -------------------- ###
  sink(file = summary.path, append = TRUE)
  
  printf("### Start %s with %s on dataset=%s comment.in=%s subset=%s n.trees=%s shrinkage=%s ###", timeStart, model[model.id], dataset[dataset.id], comment.in, subset[subset.id], ntree.val, shrink.val)
    printf("Sampling size: %s with seed ID: %s", sampling.size, seed.vec[it.total + 1])
    printf("Random sampling of condition: %s", randomize.label)
    printf("Introduce Spike Ins: %s", spikeIns)
    printf("Selected Threshold: %s", selected.thresh)
    if (pVal.opt) {
      printf("optimal p-Value search via log-rank: %s", pVal.opt)
    } else {
      printf("Selected Parameters for FP rate <= %s and for SENS >= %s", FP.value, SENS.value)
      if (thresh.alert) printf("!!!! ALERT !!!! No threshold found for above FP/SENSE setting. Take first position of highest AUC!")
      printf("Possible thresholds for this configuration: %s", paste(op.thresholds, collapse= " "))
    }
    printf("For this threshold: id=%s with RR value %s", local.thresh, op.thresh)
    print(table(pred.total.class, unlist(df.total[-oob, status.idx])))
    printf("Training prediction error=%s/%s rate=%0.3f", error.pred, sum(CM.total), (error.pred / (sum(CM.total))))
    printf("OOB prediction error=%s/%s rate=%0.3f", error.oob, length(error.oob), (abs(error.oob) / length(error.oob)))
    printf("RMSE(training/OOB/total)=%0.3f/%0.3f/%0.3f", RMSE.training, RMSE.OOB, RMSE.total)
    print("###")
    if (oob == nrow(df.total)) {
      print("Global results:")
      printf("Global training prediction error=%s/%s rate=%0.3f", sum(global.error.pred), sum(CM.total) * oob, (sum(global.error.pred) / (sum(CM.total) * oob)))
      printf("Global OOB prediction error=%s/%s rate=%0.3f", sum(abs(global.error.oob)), oob, (sum(abs(global.error.oob)) / oob))
      printf("Global total prediction error=%s/%s rate=%0.3f", sum(global.error.pred) + sum(abs(global.error.oob)), (sum(CM.total) + 1) * oob, (sum(abs(global.error.pred)) + sum(abs(global.error.oob))) / ((sum(CM.total) + 1) * oob))
      print("###")
      printf("Global mean(RMSE) training/OOB/total=%0.3f/%0.3f/%0.3f", mean(global.RMSE.training), mean(global.RMSE.OOB), mean(global.RMSE.total))
      print("###")
    }
  sink()

  printf("Done with %s OOB=%s.", model[model.id], oob)
  printf("Training prediction error=%s/%s rate=%0.3f", error.pred, sum(CM.total), (error.pred / (sum(CM.total))))
  printf("OOB prediction error=%s/%s rate=%0.3f", error.oob, length(error.oob), (abs(error.oob) / length(error.oob)))
  printf("RMSE(training/OOB/total)=%0.3f/%0.3f/%0.3f", RMSE.training, RMSE.OOB, RMSE.total)
  ###########################################################################


  if (FALSE) {
    ###########################################################################
    ### ------------------------- writing section ------------------------- ###


    ####################################################################
    ### find error in Prediction.oob
    # for DDPR
    error.ddpr <- df.total[oob, 2] - df.total[oob, 3]
    error.ddpr[error.ddpr == 1] <- "FN"
    error.ddpr[error.ddpr == - 1] <- "FP"
    error.ddpr[error.ddpr == 0] <- ""
    # for PRI model
    error.model <- df.total[oob, 2] - pred.oob.class
    error.model[error.model == 1] <- "FN"
    error.model[error.model == - 1] <- "FP"
    error.model[error.model == 0] <- ""
    # compare RealStatus / DDPRStatus / ModelStatus 
    compred.oob <- data.frame(
      PatientID = rownames(df.total)[oob], 
      RealStatus = df.total[oob, 2], 
      DDPRStatus = df.total[oob, 3], 
      DDPRError = error.ddpr, 
      ModelStatus = pred.oob.class, 
      ModelError = error.model)
    ### find error in Prediction.oob
    # for DDPR
    error.ddpr <- df.total[-oob, 2] - df.total[-oob, 3]
    error.ddpr[error.ddpr == 1] <- "FN"
    error.ddpr[error.ddpr == - 1] <- "FP"
    error.ddpr[error.ddpr == 0] <- ""
    # for PRI model
    error.model <- df.total[-oob, 2] - pred.total.class
    error.model[error.model == 1] <- "FN"
    error.model[error.model == - 1] <- "FP"
    error.model[error.model == 0] <- ""
    # compare RealStatus / DDPRStatus / ModelStatus 
    compred.total <- data.frame(
      PatientID = rownames(df.total)[-oob], 
      RealStatus = df.total[-oob, 2], 
      DDPRStatus = df.total[-oob, 3], 
      DDPRError = error.ddpr, 
      ModelStatus = pred.total.class, 
      ModelError = error.model)
    ### bind training and validation error
    compred.combined <- bind_rows(compred.total, compred.oob)
    ### save 
    write.xlsx(compred.combined, comp.total.path)


    # saving all FP - rates, SENS, SPEC and associated thresholds
    if (output.model){
      comp.all <- data.frame(data.frame(
        threshold = threshold * max(pred.total)), 
        SENS = data.frame(all.sens), 
        SPEC = data.frame(all.spec), 
        'FP-rate' = data.frame(all.FP.rate))
      write.xlsx(comp.all, comp.SESPFP.path)
    }

    if (output.model) {
      ### get active coeffizients!
      all.active.names <- top.coef.value$var
      all.active.names.split <- strsplit(as.character(all.active.names)[1:length(all.active.names)], ".", fixed = TRUE)
      # get single elements in variables
      x <- y <- z <- modus <- quad <- active.idx <- relinfo <- vector()
      i <- 0
      while (i < length(all.active.names)){
        i <- i + 1
        active.idx <- c(active.idx, which(colnames(df.total) == all.active.names[i]))
        x <- c(x, all.active.names.split[[i]][1])
        y <- c(y, all.active.names.split[[i]][2])
        z <- c(z, all.active.names.split[[i]][3])
        modus <- c(modus, all.active.names.split[[i]][4])
        quad <- c(quad, all.active.names.split[[i]][5])
        relinfo <- c(relinfo, top.coef.value$rel.inf[i])
      }
      ### ------------ safe coef names and values for op.fit as excel
      df.active <- df.total[, which(colnames(df.total) %in% all.active.names)]
      rm(all.active.names, all.coef.value)

      #Relapsstatus as Vector Red / blue for Heatmap
      relapse.status <- as.numeric(df.total$RelapseStatus)
      relapse.status[relapse.status == 1] <- "red"
      relapse.status[relapse.status == 0] <- "blue"

      #calculate mean and var for training and validation 
      all.relapsed <- cohort$PatientID[which(cohort$RelapseStatus == "Yes")]
      all.relapsefree <- cohort$PatientID[which(cohort$RelapseStatus == "No")]

      all.mean.range.relapsed <- all.var.range.relapsed <- all.mean.range.relapsefree <- all.var.range.relapsefree <- vector()
      it <- 0
      for (i in 1:ncol(df.active)){
        it <- it  + 1
        all.mean.range.relapsed <- c(all.mean.range.relapsed, mean(df.active[which(row.names(df.active) %in% all.relapsed), it]))
        all.var.range.relapsed <- c(all.var.range.relapsed, var(df.active[which(row.names(df.active) %in% all.relapsed), it]))
        all.mean.range.relapsefree <- c(all.mean.range.relapsefree, mean(df.active[which(row.names(df.active) %in% all.relapsefree), it]))
        all.var.range.relapsefree <- c(all.var.range.relapsefree, var(df.active[which(row.names(df.active) %in% all.relapsefree), it]))
      }

      ### safe all real absRange values 
      #transponse df.active 
      df.active.t <- t(df.active)
      # colnames(df.active.t) <- row.names(df.active)
      # row.names(df.active.t) <- c(1:dim(df.active.t)[1])

      ### safe df as excel
      result.data <- data.frame(active.idx, x, y, z, modus, quad, relinfo, 
        all.mean.range.relapsed, all.var.range.relapsed, 
        all.mean.range.relapsefree, all.var.range.relapsefree, df.active.t)
      colnames(result.data) <- c("Variable.ID", "X", "Y", "Z", "Feature", "Quadrant", "relInfo", "Mean(Relapsed)", "Var(Relapsed)", "Mean(Relapsefree)", "Var(Relapsefree)", colnames(df.active.t))
      rownames(result.data) <- top.coef.value[, 1]
      write.xlsx(result.data, variables.path) 
    }
    ###########################################################################

    ### ------------------------- plotting section ------------------------- ###
    # Plot relative influence of each variable
    if (output.model) {
      pdf(varimp.path, width = 12)
      # par(mfrow = c(1, 2), oma = 4, 1, 1, 1)
      par(mfrow = c(1, 1), oma = c(4, 1, 1, 1), mar=c(3, 3, 4, 2))
      # summary(gbm.cox, n.trees = 1, cBars = 20) # using first tree
      summary(gbm.cox, n.trees = best.iter, cBars = 20) # using estimated best 
      gbm.perf(gbm.cox, method="OOB", plot.it = TRUE, oobag.curve = TRUE)
      dev.off()
    }
    # save ROC plot as pdf
    if (output.model){
      pdf(roc.pdf.path)
        par(mgp = c(2, 0.5, 0), pty = "s", mar = c(4, 4, 3, 1), cex.lab=1.2)
        plot(all.FP.rate, all.sens, type = "l", ylab = "TP Rate (Sensitivity)", 
        xlab = "FP Rate (1-Specificity)", ylim = c(0, 1), xlim = c(0, 1), main = "ROC Curve")
      dev.off()
    }
    #create kaplan-meier-survival Kurve
    if (output.model){
      time <- c(df.total[-oob, 1], df.total[oob, 1])
      status <- c(df.total[-oob, 2], df.total[oob, 2])
      df.new <- data.frame(time, status, compred.combined[, 5])
      km.type = survfit(Surv(df.new[, 1], df.new[, 2]) ~ df.new[, 3], data = df.new, type = "kaplan-meier")
      ggsurvplot(km.type, conf.int = TRUE,
        legend.labs = c("low Risk", "high Risk"),
        ggtheme = theme_minimal(),
        pval = TRUE,
        pval.method = TRUE,
        risk.table = FALSE)
      ggsave(kaplan.path)
    }
    ###########################################################################
  }
}
stopCluster(cl)

RMSE <- list(
  RMSE.training = global.RMSE.training, 
  RMSE.OOB = global.RMSE.OOB, 
  RMSE.total = global.RMSE.total)
saveRDS(RMSE, RMSElist.path)
ERROR <- list(
  error.training = global.error.pred,
  error.OOB = global.error.oob,
  error.total = global.error.pred + abs(global.error.oob)
)
saveRDS(ERROR, errorlist.path)
print("###")
printf("Done with model=%s with n.trees=%s, shrinkage=%s.", model[model.id], ntree.val, shrink.val)
print(Sys.time() - timeStart)
print("Global results:")
printf("Global training prediction error=%s/%s rate=%0.3f", sum(global.error.pred), sum(CM.total) * oob, (sum(global.error.pred) / (sum(CM.total) * oob)))
printf("Global OOB prediction error=%s/%s rate=%0.3f", sum(abs(global.error.oob)), oob, (sum(global.error.oob) / oob))
printf("Global total prediction error=%s/%s rate=%0.3f", sum(global.error.pred) + sum(abs(global.error.oob)), (sum(CM.total) + 1) * oob, (sum(abs(global.error.pred)) + sum(abs(global.error.oob))) / ((sum(CM.total) + 1) * oob))
print("###")
printf("Global mean(RMSE) training/OOB/total=%0.3f/%0.3f/%0.3f", mean(global.RMSE.training), mean(global.RMSE.OOB), mean(global.RMSE.total))
print("###")