#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

add.libraries <- c("plsRcox", "ggplot2")
lapply(add.libraries, require, character.only = TRUE)
### no need at OOB
if (exists("df.training")) rm(df.training)
if (exists("df.validation")) rm(df.validation)

print("Running plsRcox model with OOB on df.total..")

options(expressions = 5e5, scipen = 100, digits = 6)
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

status.idx <- which(colnames(df.total) == "RelapseStatus")
ddpr.idx <- which(colnames(df.total) == "DDPRStatus")
time.idx <- which(colnames(df.total) == "Survivaltime")
sample.idx <- which(colnames(df.total) == "PatientID")

ntree.val <- 1000
shrink.val <- 0.0001

### ---------------- run for prediction -------------------------------- ###
comment.model.pred <- paste(comment.out, paste0("ntree", ntree.val), paste0("shrinkage", shrink.val), sep = ".")
errorlist.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(17, 
  com=comment.model.pred))
global.error.pred <- global.error.oob <- global.RMSE.training <- global.RMSE.OOB <- global.RMSE.total <- vector()
output.thresholds <- output.model <- FALSE

for (oob in 1:nrow(df.total)) {
  # load data
  printf("Predict model=%s ntree=%s shrink=%s and oob=%s..", model[model.id], ntree.val, shrink.val, oob)
  comment.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("shrinkage", shrink.val), paste0("oob", oob), sep=".")
  model.path <- sprintf("%s/%s.rds", Output.path, output_filenaming(12, com = comment.tmp))


  df.train <- df.total[-oob, -c(status.idx, time.idx, ddpr.idx, sample.idx)]
  y.train <- df.total$RelapseStatus[-oob]
  time.train <- df.total$Survivaltime[-oob]

  df.oob <- df.total[oob, - c(status.idx, time.idx, ddpr.idx, sample.idx)]
  y.oob <- df.total$RelapseStatus[oob]
  time.oob <- df.total$Survivaltime[oob]

  # doi:10.1093/bioinformatics/btu660
  # did not work, too computational intense
  # coxpls2   Fitting a Cox-Model on PLSR components using R package pls
  # coxpls2DR Fitting a Direct Kernel PLS model on the (Deviance) Residuals
  # pls.cox1 <- coxpls2DR(Xplan = df.train, time = time.train, time2 = y.train, ncomp = 6, plot=FALSE, validation = "CV", allres = FALSE, shrinkage = shrink.val, n.trees =ntree.val)

  # df <- as.matrix(df.train)
  # does work
  # coxpls3   Fitting a Cox-Model on PLSR components using R package plsRglm
  if (dataset.id == 6) {
    pls.cox2 <- coxpls3(Xplan = df.train, time = time.train, time2 = y.train, nt = 6, validation = "CV", allres = TRUE, shrinkage = shrink.val, n.trees = ntree.val)
  } else {
    pls.cox2 <- coxpls3(~. - RelapseStatus - DDPRStatus - Survivaltime - PatientID, time = as.numeric(df.total$Survivaltime[-oob]), time2 = as.numeric(df.total$RelapseStatus[-oob]), data = df.total[-oob, ], nt = 6, validation = "CV", allres = TRUE)
  }
  # coxpls3(~. - RelapseStatus - DDPRStatus - Survivaltime - PatientID, time = df.total.pri$Survivaltime[-oob], time2 = df.total.pri$RelapseStatus[-oob], data = df.total.pri[-oob, ])

  # , plot = TRUE
  # PCs
  pls.cox2.comp <- pls.cox2$tt_pls3
  print(summary(pls.cox2))

  # does work
  # coxDKpls2DR   Fitting a PLSR model on the (Deviance) Residuals using plsRglm
  if (dataset.id == 6) {
    pls.cox3 <- coxDKpls2DR(Xplan = df.train, time = time.train, time2 = y.train, shrinkage = shrink.val, n.trees = ntree.val, nt = 6, allres=TRUE)
  } else {
    pls.cox3 <- coxDKpls2DR(~. - RelapseStatus - DDPRStatus - Survivaltime - PatientID, time = as.numeric(df.total$Survivaltime[-oob]), time2 = as.numeric(df.total$RelapseStatus[-oob]), data = df.total[-oob, ], nt = 6, allres=TRUE)
  }


  predict(pls.cox2$cox_pls3, type = "survival", comps = 5)
  # type = "lp", "risk", "expected" (class), "terms", "survival"
  pred.pls.cox3 <- round(predict(pls.cox3, type = "expected"))
  CM.total <- table(pred.pls.cox3, as.numeric(df.total[-oob, status.idx]))
  # confusion matrix
  print(CM.total)

  ### another function
  pls.cox4 <- coxDKpls2DR(~. - RelapseStatus - DDPRStatus - Survivaltime - PatientID, time = as.numeric(df.total$Survivaltime), time2 = as.numeric(df.total$RelapseStatus), data = df.total, nt = 6, allres = TRUE, validation = "LOO")
  
  pred.pls.cox4 <- round(predict(pls.cox4$cox_DKpls2DR, type = "expected"))
  CM.total <- table(pred.pls.cox4, as.numeric(df.total[, status.idx]))
  # confusion matrix
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
    printf("Predict model=%s with n.trees=%s, shrinkage=%s and oob=%s with plotting..", model[model.id], ntree.val, shrink.val, oob)
    comment.tmp <- paste(comment.out, paste0("ntree", ntree.val), paste0("shrinkage", shrink.val), paste0("oob", oob), sep=".")

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
        printf("### Start %s_OOB with %s on dataset=%s and subset=%s n.trees=%s, shrinkage=%s, random=%s spikes=%s OOB=%s###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], ntree.val, shrink.val, randomize.label, spikeIns, oob)
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
  printf("### Start %s_OOB with %s on dataset=%s and subset=%s n.trees=%s, shrinkage=%s, random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], ntree.val, shrink.val, randomize.label, spikeIns)
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
  
  RMSE.training = global.RMSE.training, 
  RMSE.OOB = global.RMSE.OOB, 
  RMSE.total = global.RMSE.total)
)
saveRDS(ERROR, errorlist.path)