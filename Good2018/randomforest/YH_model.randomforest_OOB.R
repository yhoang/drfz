#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

# add.libraries <- c("caret", "randomForest", "doParallel", "e1071", "ggplot2")
add.libraries <- c("randomForest", "doParallel", "ggplot2")
lapply(add.libraries, require, character.only = TRUE)
### no need at OOB
if (exists("df.training")) rm(df.training)
if (exists("df.validation")) rm(df.validation)

print("Running random forest model with OOB on df.total..")

options(expressions = 5e5, digits=4)
memory.limit(size = 1e8)
### in command
# cd "C:\Program Files\R\R-3.6.1\bin"
# ./R.exe  --max-ppsize 500000
# source(file.path("D:", "drfz", "Good2018", "YH_run.R"))

### ------------------ random forest model 
it.total <- 0

timeStart <- Sys.time()
ptm <- proc.time()
printf("### Start %s with %s on dataset=%s comment.in=%s subset=%s  random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], comment.in, subset[subset.id], randomize.label, spikeIns)

set.seed(seed.vec[it.total + 1])

if (FALSE) {
  ############################################################################
  #### FIND BEST ntree/mtry
  # tuning cross-validation
  # creating stratified/balanced indices
  cv.folds <- createMultiFolds(df.total[, 2], k = 4, times = 3)
  # creating folds 
  trainctrl <- trainControl(method = "cv", number = 4, repeats = 3, index = cv.folds, allowParallel=TRUE)

  ### Create a Random Forest model with default parameters
  # ntree: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
  # mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 
  # default: ntree=500, mtry=floor(sqrt(ncol(df.total)))

  # initiate accuracy collection
  max.acc <- list()
  best.mtry <- vector()
  acc.it <- 0

  # model tuning matrix
  rf.grid <- expand.grid(mtry = c(2, 10, 50))
  trees <- c(5000, 7000)
  ### adding cores
  cluster.size <- 3
  cl <- makeCluster(cluster.size)
  registerDoParallel(cl)
  for (it in 1:length(trees)) {
    acc.it <- acc.it + 1
    printf("%s/%s:: do rf with ntree=%s and mtry=%s. acc.it=%s", it, length(trees), trees[it], paste(rf.grid, collapse=" "), acc.it)
    rf_tree <- train(form = as.factor(RelapseStatus) ~ . - Survivaltime - DDPRStatus - PatientID, data = df.total, method = "rf", ntree = trees[it], trControl = trainctrl, tuneGrid = rf.grid, do.trace = 500, proximity = TRUE)
    
    curr.max.mtry <- rf_tree$bestTune
    best.mtry <- c(best.mtry, curr.max.mtry)

    curr.max.acc <- max(rf_tree$results$Accuracy)

    printf("Best mtry=%s acc=%s", curr.max.mtry, curr.max.acc)
    max.acc[[acc.it]] <- c(Acc=curr.max.acc, ntree=trees[it], curr.max.mtry)
  }
  stopCluster(cl)
  saveRDS(max.acc, paste0(Output.path, "max.acc.list.rds"))

  print("mtries used in models:")
  print(table(unlist(best.mtry)))

  best.acc <- best.acc.id <- 0
  for (i in 1:acc.it) {
    print(max.acc[[i]]$Acc)
    if (max.acc[[i]]$Acc > best.acc){
      best.acc <- max.acc[[i]]$Acc
      best.acc.id <- i
      # print(best.acc)
      # print(best.acc.id)
    }
  }
  printf("Highest Acc=%#.4f at ntree=%s and mtry=%s", max.acc[[best.acc.id]]$Acc, max.acc[[best.acc.id]]$ntree, max.acc[[best.acc.id]]$mtry)

  rf.ntree1250.mtry5 <- rf_tree
  rf.ntree1250.mtry10 <- rf_tree
  rf.ntree3000.mtry5 <- rf_tree
  rf.ntree7000.mtry10 <- rf_tree

  ### comparison of models
  modelcomp.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(14, com=comment.tmp))
  resamps <- resamples(list(
    rf.ntree1250.mtry5 = rf.ntree1250.mtry5,
    rf.ntree1250.mtry10 = rf.ntree1250.mtry10,
    rf.ntree3000.mtry5 = rf.ntree3000.mtry5, 
    rf.ntree7000.mtry10 = rf.ntree7000.mtry10))
  summary(resamps)
  pdf(modelcomp.path)
    dotplot(resamps, main = "model comparison")
  dev.off()

  
  # not working
  # rf_tree <- train(form = as.factor(RelapseStatus) ~ ., data = df.training[, -c(1,3)], method = "ranger", trControl = trainctrl, num.trees = 500, tuneGrid = rf.grid, importance = "impurity", do.trace = 1000, proximity = TRUE)
  # is working
  # rf_tree2 <- train(x = df.training[, - c(1:3)], y = as.factor(df.training[, 2]), method = "rf", trControl = trainctrl, tuneGrid = rf.grid, do.trace = 1000, proximity = TRUE)
  # is working
  # rf_tree <- train(x = df.training[, - c(1:3)], y = as.factor(df.training[, 2]), method = "rf", trControl = trainctrl, do.trace = 1000)
}

### find/set best setting
# best.acc <- sapply(max.acc, max)
ntree.val <- 500
mtry.val <- floor(sqrt(ncol(df.total)))
# model tuning matrix
rf.grid <- expand.grid(mtry = mtry.val)

status.idx <- which(colnames(df.total) == "RelapseStatus")
ddpr.idx <- which(colnames(df.total) == "DDPRStatus")
time.idx <- which(colnames(df.total) == "Survivaltime")
sample.idx <- which(colnames(df.total) == "PatientID")

###------------------ run random forest with fixed parameters -------------###
if (FALSE) {
  cluster.size <- 3
  cl <- makeCluster(cluster.size)
  registerDoParallel(cl)
  ptm <- proc.time()

  # for (oob in 1:nrow(df.total)) {
  
  printf("Run model=%s with n.trees=%s, mtry=%s and oob=%s..", model[model.id], ntree.val, mtry.val, oob)
  comment.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("mtry", mtry.val), paste0("oob", oob), sep=".")
  model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, com=comment.tmp))

  if (!file.exists(model.path)) {

    if (TRUE) {
      # run model and save
      df.train <- df.total[-oob, - c(status.idx, time.idx, ddpr.idx, sample.idx)]
      y.train <- df.total$RelapseStatus[-oob]
      df.oob <- df.total[oob, - c(status.idx, time.idx, ddpr.idx, sample.idx)]
      y.oob <- df.total$RelapseStatus[oob]
      
      rf.tree <- randomForest(x = df.train, y = as.factor(y.train), ntree = ntree.val, mtry = mtry.val)
      
    } else {
      # tuning cross-validation
      # cv.folds <- createMultiFolds(df.total[-oob, 2], k = 4, times = 3)
      # # creating folds 
      # trainctrl <- trainControl(method = "cv", number = 4, repeats = 3, index = cv.folds, allowParallel=TRUE)
      trainctrl <- trainControl(allowParallel = TRUE)

      rf.tree <- train(form = as.factor(RelapseStatus) ~ . - Survivaltime - DDPRStatus - PatientID, data = as.data.frame(df.total[-oob,]), method = "rf", ntree = ntree.val, trControl = trainctrl, tuneGrid = rf.grid, do.trace = 500, proximity = TRUE)
      #tuneGrid = rf.grid, proximity = TRUE)
    }
    saveRDS(rf.tree, model.path)
    # rm(rf.tree)
    # garbage collection
    gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  } else {
    printf("File exists already: %s", model.path)
  }
# }
printf("Done with model=%s with n.trees=%s, mtry=%s and oob=%s..", model[model.id], ntree.val, mtry.val, oob)
print(proc.time() - ptm)
stopCluster(cl)
}

### ---------------- run for prediction -------------------------------- ###
comment.model.pred <- paste(comment.out, paste0("ntree", ntree.val), paste0("mtry", mtry.val), sep = ".")
errorlist.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(17, 
  com=comment.model.pred))
global.error.pred <- global.error.oob <- global.error.pred2 <- global.error.oob2 <- global.RMSE.training <- global.RMSE.OOB <- global.RMSE.total <- vector()
output.thresholds <- output.model <- FALSE
for (oob in 1:nrow(df.total)) {
  # load data
  printf("Predict model=%s with n.trees=%s, mtry=%s and oob=%s..", model[model.id], ntree.val, mtry.val, oob)
  comment.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("mtry", mtry.val), paste0("oob", oob), sep=".")
  model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, com = comment.tmp))


  df.train <- df.total[-oob, -c(status.idx, time.idx, ddpr.idx, sample.idx)]
  y.train <- df.total$RelapseStatus[-oob]
  df.oob <- df.total[oob, - c(status.idx, time.idx, ddpr.idx, sample.idx)]
  y.oob <- df.total$RelapseStatus[oob]

  # rf.tree <- readRDS(model.path)
  rf.tree <- randomForest(x = df.train, y = as.factor(y.train), ntree = ntree.val, mtry = mtry.val)
      
  ### predict class with model
  # pred.total2 <- predict(rf.tree, newdata = df.total[-oob, ])
  # pred.oob2 <- predict(rf.tree, newdata = df.total[oob, ])
  pred.total2 <- predict(rf.tree, df.train)
  pred.oob2 <- predict(rf.tree, df.oob)
  pred.oob <- predict(rf.tree, df.oob, type="prob")
  pred.total <- predict(rf.tree, df.train, type="prob")
    
  # MSE <- mean(((as.numeric(pred.oob) - 1) - df.total[oob, 2]) ^ 2)
  # https://stackoverflow.com/questions/38621556/advice-on-usage-of-dplyr-do-vs-purrr-map-tidynest-for-predictions

  ###########################################################################
  ### Create Confusion Matrix
  # print(confusionMatrix(data = pred.total2,  
  #                 reference = as.factor(df.total[-oob, 2]),
  #                 positive = "1"))
  # print(confusionMatrix(data = pred.oob2,  
  #                 reference = factor(df.total[oob, 2], levels = c(0,1)),
  #                 positive = "1"))
  ###########################################################################
  
  ### determine the misclassification rate. 
  # First, build a confusion matrix. Each column of the matrix represents the number of predictions of each class, while each row represents the instances in the actual class.
  CM.total <- table(pred.total2, y.train)
  print(CM.total)
  # collect error count
  error.pred <- sum(CM.total) - sum(diag(CM.total))
  error.oob <- y.oob - as.numeric(as.character(pred.oob2))

  ### collect RMSE
  RMSE.training <- sqrt(mean((pred.total[, 2] - y.train) ^ 2))
  RMSE.OOB <- sqrt(mean((pred.oob[, 2] - y.oob) ^ 2))
  RMSE.total <- sqrt(mean((c(pred.total[, 2], pred.oob[, 2]) - c(y.train, y.oob)) ^ 2))

  ### collect for every run
  global.error.pred <- c(global.error.pred, error.pred)
  global.error.oob <- c(global.error.oob, error.oob)
  global.error.pred2 <- c(global.error.pred2, pred.total)
  global.error.oob2 <- c(global.error.oob2, pred.oob[, 1])
  global.RMSE.training <- c(global.RMSE.training, RMSE.training)
  global.RMSE.OOB <- c(global.RMSE.OOB, RMSE.OOB)
  global.RMSE.total <- c(global.RMSE.total, RMSE.total)
}

#############################################################################
if (FALSE) {
  ### run with plotting and writing table
  error.pred <- error.oob <- vector()
  for (oob in 1:nrow(df.total)) {
    # load data
    printf("Predict model=%s with n.trees=%s, mtry=%s and oob=%s with plotting..", model[model.id], ntree.val, mtry.val, oob)
    comment.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("mtry", mtry.val), paste0("oob", oob), sep=".")

    model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, com=comment.tmp))
    rf.tree <- readRDS(model.path)

    ### predict class with model
    pred.total2 <- predict(rf.tree, newdata = df.total[-oob, ])
    pred.oob2 <- predict(rf.tree, newdata = df.total[oob, ])
    
    ### determine the misclassification rate. 
    # First, build a confusion matrix. Each column of the matrix represents the number of predictions of each class, while each row represents the instances in the actual class.
    CM.total <- table(pred.total2, df.total[-oob, 2])
    # collect error count
    error.pred <- c(error.pred, sum(CM.total)-sum(diag(CM.total)))
    error.oob <- c(error.oob, df.total[oob, 2] - as.numeric(as.character(pred.oob2)))

    ### predict probability with model
    pred.total <- predict(rf.tree, newdata = df.total[-oob, ], type ="prob")
    pred.oob <- predict(rf.tree, newdata = df.total[oob, ], type ="prob")

    if (output.thresholds){
      threshold <- seq(0, 1, 0.01)
      predictions.roc <- data.frame()
      #safe all prediction cutoffs as df  
      for (it in 1:length(threshold)){
        newline <- findInterval(pred.total[, 2], threshold[it])
        predictions.roc <- bind_cols(as.data.frame(newline), predictions.roc)
        names(predictions.roc)[1]= threshold[it]
      }
      #mirror df to set index right 
      predictions.roc <- predictions.roc[, ncol(predictions.roc):1]

      ###calculate error
      #true prositive = prediction 1 & relapse status 1
      #calculate ssensitivity and false positiv rate 
      all.spec <- all.sens <- all.FP.rate <- AUC <- vector()
      for (i in 1:length(threshold)){
        TP <- length(which(predictions.roc[, i] == 1 & df.total[-oob, 2] == 1))
        TN <- length(which(predictions.roc[, i] == 0 & df.total[-oob, 2] == 0))
        FP <- length(which(predictions.roc[, i] == 1 & df.total[-oob, 2] == 0))
        ALLP <- length(which(df.total[-oob, 2] == 1))
        ALLN <- length(which(df.total[-oob, 2] == 0))
        SENS <- TP / ALLP
        SPEC <- TN / ALLN
        all.spec <- c(all.spec, SPEC)
        all.sens <- c(all.sens, SENS)
        FP.rate <- FP / ALLN
        all.FP.rate <- c(all.FP.rate, FP.rate)
        AUC <- c(AUC, auc(roc(df.total[-oob, 2], predictions.roc[, i], quiet=TRUE)))
      }

      # saving all FP - rates, SENS, SPEC and associated thresholds
      comp.all <- data.frame(
        threshold = threshold,
        SENS = data.frame(all.sens), 
        SPEC = data.frame(all.spec), 
        'FP-rate' = data.frame(all.FP.rate))
      write.xlsx(comp.all, comp.SESPFP.path)


      if (pVal.opt) {
        ##########################################################################
        # #treshold by log - rang test. find "optimal" p - value for log rank
        p.value <- c()
        for (i in 1:ncol(predictions.roc)) {
          df.new <- data.frame(
            Survivaltime = df.total[-oob, 1], 
            RelapseStatus = df.total[-oob, 2], 
            Prediction = predictions.roc[, i])
          
          # calculate estimate of survival curve for censored data using Kaplan-Meier 
          km.type = survfit(
            Surv(Survivaltime, RelapseStatus) ~ Prediction,
            data = df.new,
            type = "kaplan-meier")
          # or try coxph( Surv(df.new$Time, df.new$Status) ~ df.new$Prediction, data = df.new)
          # calculate log rank test and extract p-value
          p.value <- c(p.value, surv_pvalue(km.type, method = "1")$pval)
        }
        p.value[which(is.na(p.value) == TRUE)] <- 1
        op.thresh <- threshold[which(p.value == min(p.value))[1]]
        ##########################################################################
      } else {
        op.thresh = 0.5
      }
      printf("RR over %0.3f are interpret as Status 1 (High Risk) and under as Status 0 (Low Risk)", op.thresh)
      pred.oob.class <- findInterval(pred.oob[, 2], op.thresh)
      pred.total.class <- findInterval(pred.total[, 2], op.thresh)

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

      ###----------------------- output file names --------------------------###
      #### variables xlsx
      variables.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(7, com = comment.tmp))
      # VarImportance PDF
      varimp.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(13, com = comment.tmp))
      # summary txt
      summary.path <- sprintf("%s/%s.log", Output.path, output_filenaming(11, com = comment.tmp))
      ##########################################################################

      ### save 
      write.xlsx(compred.combined, comp.total.path)
      ##########################################################################
      ### write summary
      sink(file = summary.path, append = TRUE)
        printf("### Start %s_OOB with %s on dataset=%s and subset=%s n.trees=%s, mtry=%s, random=%s spikes=%s OOB=%s###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], ntree.val, mtry.val, randomize.label, spikeIns, oob)
        printf("Sampling size: %s with seed ID: %s", sampling.size, seed.vec[it.total + 1])
        printf("Random sampling of condition: %s", randomize.label)
        printf("Introduce Spike Ins: %s", spikeIns)
        if (pVal.opt) {
          printf("Optimal p-Value search via log-rank: %s", pVal.opt)
        }
        printf("RR for this threshold: %s", op.thresh)
        print("###")
      sink()
      ##########################################################################
          
      ### sort variables by 3 MeanDecreaseAccuracy or 4 MeanDecreaseGini
      #and select > 0
      imp.idx <- 4
      rf.sort.idx <- order(rf$importance[, imp.idx], decreasing = TRUE)
      rf.imp.val <- rf$importance[rf.sort.idx, ]
      rf.imp.idx <- 1:20
      rf.imp.val <- rf.imp.val[rf.imp.idx, ]
      rf.imp.val.split <- strsplit(rownames(rf.imp.val), ".", fixed = TRUE)

      ###############
      #get active coefficients!
      x <- y <- z <- modus <- quad <- vector()
      i <- 0
      while (i < nrow(rf.imp.val)){
        i <- i + 1
        x <- c(x, rf.imp.val.split[[i]][1])
        y <- c(y, rf.imp.val.split[[i]][2])
        z <- c(z, rf.imp.val.split[[i]][3])
        modus <- c(modus, rf.imp.val.split[[i]][4])
        quad <- c(quad, rf.imp.val.split[[i]][5])
      }
      result.data <- cbind(rf.imp.val[, 4], x, y, z, modus, quad)
      write.xlsx(result.data, variables.path) 
      
    }

    #################################################################
    ### plot section
    if (output.model) {

      pdf(varimp.path, width=20)
      varImpPlot(rf, n.var = 20, main = "Top 20 Variable Importance")
      dev.off()

      #create AUC
      pdf(auc.pdf.path)
        plot(roc(compred.combined[, 2], compred.combined[, 5]))
      dev.off()
      AUC <- auc(roc(compred.combined[, 2], compred.combined[, 5]))
      print(AUC)

      # create iAUC
      iAUC <- AUC.uno(Surv(df.train[, 1], df.train[, 2]), Surv(df.validation[, 1], df.validation[, 2]), pred.validation, times = seq(10, 1000, 10), savesensspec = TRUE)

      # save  ROC plot of training as pdf
      pdf(roc.pdf.path)
        par(mgp = c(2, 0.5, 0), pty = "s", mar = c(4, 4, 3, 1), cex.lab=1.2)
        plot(all.FP.rate, all.sens, type = "l", ylab = "TP Rate (Sensitivity)", 
        xlab = "FP Rate", ylim = c(0, 1), xlim = c(0, 1), main = "ROC Curve Training")
      dev.off()

      #create kaplan - meier - survival Kurve
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

  }
#############################################################################
}

### write summary
sink(file = summary.path, append = TRUE)
  printf("### Start %s_OOB with %s on dataset=%s and subset=%s n.trees=%s, mtry=%s, random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], ntree.val, mtry.val, randomize.label, spikeIns)
  printf("Sampling size: %s with seed ID: %s", sampling.size, seed.vec[it.total])
  printf("Random sampling of condition: %s", randomize.label)
  printf("Introduce Spike Ins: %s", spikeIns)
  printf("Done with %s OOB.", model[model.id])
  
  print("###")
  print("Global results:")
  printf("Global training prediction error=%s/%s rate=%0.3f", sum(global.error.pred), sum(CM.total) * oob, (sum(global.error.pred) / (sum(CM.total) * oob)))
  printf("Global OOB prediction error=%s/%s rate=%0.3f", sum(abs(global.error.oob)), oob, (sum(abs(global.error.oob)) / oob))
printf("Global total prediction error=%s/%s rate=%0.3f", sum(global.error.pred) + sum(abs(global.error.oob)), (sum(CM.total) + 1) * oob, (sum(abs(global.error.pred)) + sum(abs(global.error.oob))) / ((sum(CM.total) + 1) * oob))
  print("###")
  printf("Global mean(RMSE) training/OOB/total=%0.3f/%0.3f/%0.3f", mean(global.RMSE.training), mean(global.RMSE.OOB), mean(global.RMSE.total))
  print("###")
  print(Sys.time() - timeStart)
  print("###")
sink()
###########################################################################

printf("Done with %s OOB.", model[model.id])
print("Global results:")
  printf("Global training prediction error=%s/%s rate=%0.3f", sum(global.error.pred), sum(CM.total) * oob, (sum(global.error.pred) / (sum(CM.total) * oob)))
  printf("Global OOB prediction error=%s/%s rate=%0.3f", sum(abs(global.error.oob)), oob, (sum(abs(global.error.oob)) / oob))
printf("Global total prediction error=%s/%s rate=%0.3f", sum(global.error.pred) + sum(abs(global.error.oob)), (sum(CM.total) + 1) * oob, (sum(abs(global.error.pred)) + sum(abs(global.error.oob))) / ((sum(CM.total) + 1) * oob))
  print("###")
  printf("Global mean(RMSE) training/OOB/total=%0.3f/%0.3f/%0.3f", mean(global.RMSE.training), mean(global.RMSE.OOB), mean(global.RMSE.total))
  print("###")
print(Sys.time() - timeStart)


ERROR <- list(
  error.training = global.error.pred,
  error.OOB = global.error.oob,
  error.total = global.error.pred + abs(global.error.oob),
  prob.training = global.error.pred2,
  prob.oob = global.error.oob2,
  RMSE.training = global.RMSE.training, 
  RMSE.OOB = global.RMSE.OOB, 
  RMSE.total = global.RMSE.total
)
saveRDS(ERROR, errorlist.path)

RMSE.oob.zscore <- (ERROR$RMSE.OOB - mean(ERROR$RMSE.OOB)) / sd(ERROR$RMSE.OOB)

hist(RMSE.oob.zscore)