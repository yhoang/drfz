#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

add.libraries <- c("glmnet", "doParallel")
lapply(add.libraries, require, character.only = TRUE)

print("Running cox hazard model..")
### ------------------ cox hazard model in cross validation
### initiate variables for cross validation / cox model

it.total <- 0
all.activ.Index <- all.coef.value <- vector()
variabl.count <- vector()
min.cvm <- 100

cl <- makeCluster(cluster.size)
registerDoParallel(cl)

timeStart <- Sys.time()
ptm <- proc.time()
printf("### Start %s with %s on dataset=%s and subset=%s lambda=%s random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], lambda[lambda.id], randomize.label, spikeIns)


while (it.total < sampling.size) {
  it.total <- it.total + 1 

  set.seed(seed.vec[it.total])

  ### --------------------- create fold values
  # Basal has 13 relapsed and 31 relapse free patients
  # set folds for cross validation manual because of imbalance data
  # set folds that 1 fold contains at least 3 relapsed patients
  # num.zero <- table(df.training[, 2])[1]
  # num.one <- table(df.training[, 2])[1]
  zero.id <- one.id <- vector()
  fold.id <- rep(0, nrow(df.training))
  one.id <- c(rep(sample(1:4), 3), 3)
  zero.id <- rep(sample(1:4), 8)
  length(one.id)
  it.set <- it.null <- it.one <- 1

  #create 4 folds with at least 3 relapsed patients
  while (it.set < nrow(df.training) + 1){
      if (df.training[it.set, 2] == 0) {
      # if sample relapse status is false / 0
        fold.id[it.set] <- zero.id[it.null]
        it.null <- it.null + 1
    } else {
      fold.id[it.set] <- one.id[it.one]
      it.one <- it.one  + 1
    }
    it.set <- it.set  + 1
  }
  ##############
  
  ### ----------- create model with cross validation
  cv.fit <- cv.glmnet(df.training[, -c(1:3)],
                      Surv(df.training[, 1], df.training[, 2]), 
                      family = "cox",
                      alpha = set.alpha, 
                      foldid = fold.id, 
                      parallel = TRUE)
  #plot(cv.fit)

  ### collect all active (!=0) coef from fited model
  # lambda.min <- model with min cross validation error
  # lambda.1se <- model with min cross validation error + 1 standard error
  if (lambda.id == 1) {
    Coefficients <- coef(cv.fit, s=cv.fit$lambda.min)
    coef.value <- coef(cv.fit, s=cv.fit$lambda.min)
  } else {
    Coefficients <- coef(cv.fit, s=cv.fit$lambda.1se)
    coef.value <- coef(cv.fit, s=cv.fit$lambda.1se)
  }
  Active.Index <- which(Coefficients != 0)
  all.coef.value <- c(all.coef.value, coef.value)
  variabl.count <- c(variabl.count, length(Active.Index))
  printf("%s PRI features chosen.", length(Active.Index))
  tmp <- min(cv.fit$cvm)

  ### collect best model (model with min cross validation error) 
  # if variabl.count is not empty!
  # safe best model as op.fit
  # safe active coef
  if ((min.cvm > tmp) & (length(Active.Index) > 0)){
    min.cvm <- tmp
    op.fit <- cv.fit
    op.variabl.count <- length(Active.Index)
    op.index <- Active.Index

    }
  if (it.total %% 10 == 0) {
    printf("At Work IT:%s ", it.total)
    print(Sys.time() - timeStart)
    print(proc.time() - ptm)
  }
}
printf("Done with %s.", model[model.id])
print(Sys.time() - timeStart)

stopCluster(cl)

print(table(variabl.count))

if (all(variabl.count == 0)) {
  sink(file = summary.path, append = TRUE)
    print(timeStart)
    printf("Sampling size: %s with seed ID: %s", sampling.size, seed.vec[it.total])
    print("Variable distribution in models:")
    print(table(variabl.count))
    printf("Random sampling of condition: %s", randomize.label)
    printf("Introduce Spike Ins: %s", spikeIns)
    printf("optimal p-Value search via log-rank: %s", pVal.opt)
    print("No variables found. Stop.")
  sink()
  next

} else {

  printf("Best model >0: %s features", op.variabl.count)
  printf("Cross validation error: %s", min.cvm)

  #prediction for training and validation set based on op.fit
  #prediction for type cox model = relativ risk (RR)
  if (lambda.id == 1) {
    pred.training <- predict(op.fit, newx = df.training[, -c(1:3)], s="lambda.min", type="response")
    pred.validation <- predict(op.fit, newx = df.validation[, -c(1:3)], s="lambda.min", type="response")
  } else {
    pred.training <- predict(op.fit, newx = df.training[, -c(1:3)], s="lambda.1se", type="response")
    pred.validation <- predict(op.fit, newx = df.validation[, -c(1:3)], s="lambda.1se", type="response")
  }


  ### --------------- calculat threshold for cutoff #####################
  # "low Risk" < treshold > "high risk"
  ### scale RR between 0 -> 1
  scale.prediction <- (pred.training - min(0)) / (max(pred.training) - min(0))
  scale.prediction[which(scale.prediction < 0)] <- 0 
  #performance test for set of thresholds
  threshold <- seq(0, 1, 0.01)
  predictions.roc <- data.frame()
  #safe all prediction cutoffs as df  
  for (it in 1:length(threshold)){
    newline <- findInterval(scale.prediction, threshold[it])
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
    printf("Sampling size: %s with seed ID:%s", sampling.size, seed.vec[it.total])
    print("Variable distribution in models:")
    print(table(variabl.count))
    printf("Best model >0: %s features", op.variabl.count)
    printf("Random sampling of condition: %s", randomize.label)
    printf("Introduce Spike Ins: %s", spikeIns)
    if (pVal.opt) {
      printf("optimal p-Value search via log-rank: %s", pVal.opt)
    } else {
      printf("Selected Parameters for FP rate <= %s and for SENS  >=  %s", FP.value, SENS.value)
      printf("Possible thresholds for this configuration: %s", op.thresholds)
    }
    printf("Specifically selected Threshold for this configuration: %s at position: %s", op.thresholds[selected.thresh], selected.thresh)
    printf("RR for this threshold: %s", op.thresh * max(pred.training))
    sink()
  }
  #pediction RR Values aprox Relapsstaus with calculates Threshold
  predicted.validation=findInterval(pred.validation / max(pred.training), op.thresh)
  #print(predicted.validation)
  predicted.training=findInterval(pred.training / max(pred.training), op.thresh)
  #print(predicted.training)

  ####################################################################
  # find error in Prediction.validation
  # for DDPR
  error.ddpr <- df.validation[, 2] - df.validation[, 3]
  error.ddpr[error.ddpr == 1] <- "FN"
  error.ddpr[error.ddpr == - 1] <- "FP"
  error.ddpr[error.ddpr == 0] <- ""
  # for PRI model
  error.model <- df.validation[, 2] - predicted.validation
  error.model[error.model == 1] <- "FN"
  error.model[error.model == - 1] <- "FP"
  error.model[error.model == 0] <- ""
  # compare RealStatus / DDPRStatus / ModelStatus 
  compred.validation <- data.frame(rownames.validation, df.validation[, 2], df.validation[, 3], error.ddpr, predicted.validation, error.model)
  names(compred.validation)[1] <- "Patient ID"
  names(compred.validation)[2] <- "RealStatus"
  names(compred.validation)[3] <- "DDPRStatus"
  names(compred.validation)[4] <- "DDPRErrors"
  names(compred.validation)[5] <- "ModelStatus"
  names(compred.validation)[6] <- "ModelErrors"

  ### find error in Prediction.training 
  # for DDPR
  error.ddpr <- df.training[, 2] - df.training[, 3]
  error.ddpr[error.ddpr == 1] <- "FN"
  error.ddpr[error.ddpr == - 1] <- "FP"
  error.ddpr[error.ddpr == 0] <- ""
  # for PRI model
  error.model <- df.training[, 2] - predicted.training
  error.model[error.model == 1] <- "FN"
  error.model[error.model == - 1] <- "FP"
  error.model[error.model == 0] <- ""
  # compare RealStatus / DDPRStatus / ModelStatus
  compred.training <- data.frame(rownames.training, df.training[, 2], df.training[, 3], error.ddpr, predicted.training, error.model)
  names(compred.training)[1] <- "Patient ID"
  names(compred.training)[2] <- "RealStatus"
  names(compred.training)[3] <- "DDPRStatus"
  names(compred.training)[4] <- "DDPRErrors"
  names(compred.training)[5] <- "ModelStatus"
  names(compred.training)[6] <- "ModelErrors"

  ### bind training and validation error
  comp.total <- bind_rows(compred.training, compred.validation)
  ### save
  if (output.thresholds) write.xlsx(comp.total, comp.total.path)

  ##############################################################################
  #create AUC
  if (output.thresholds){
    pdf(auc.pdf.path)
    plot(roc(comp.total[, 2], comp.total[, 5]))
    dev.off()
  }
  AUC <- auc(roc(comp.total[, 2], comp.total[, 5]))
  print(AUC)

  ########################################################################
  # create iAUC
  iAUC <- AUC.uno(Surv(df.training[, 1], df.training[, 2]), Surv(df.validation[, 1], df.validation[, 2]), pred.validation, times = seq(10, 1000, 10), savesensspec = TRUE)

  if (output.thresholds){
    pdf(iauc.pdf.path)
    names(iAUC)
    iAUC$iauc
    plot(iAUC)
    dev.off()

  }


  ##############################################################################
  #get active coeffizients!
  Coef.names <- colnames(df.training[, -c(1:3)])

  #get names from coef.names for all active index in op.fit
  all.activ.names <- colnames(df.training[, -c(1:3)])[op.index]
  all.activ.names.split <- strsplit(all.activ.names[1:length(all.activ.names)], ".", fixed = TRUE)

  x <- y <- z <- modus <- quad <- vector()
  i <- 0
  while (i < length(all.activ.names)){
    i <- i + 1
    x <- c(x, all.activ.names.split[[i]][1])
    y <- c(y, all.activ.names.split[[i]][2])
    z <- c(z, all.activ.names.split[[i]][3])
    modus <- c(modus, all.activ.names.split[[i]][4])
    quad <- c(quad, all.activ.names.split[[i]][5])
  }

  ### ------------ safe coef names and values for op.fit as excel
  df.training <- as.data.frame(df.training)
  df.validation <- as.data.frame(df.validation)

  df.training.active <- df.training[which(colnames(df.training) %in% all.activ.names)]
  df.validation.active <- df.validation[which(colnames(df.validation) %in% all.activ.names)]
  df.active <- bind_rows(df.training.active, df.validation.active)
  row.names(df.active) <- c(rownames.training, rownames.validation)

  #Relapsstatus as Vector Red / blue for Heatmap
  relapse.status <- as.numeric(c(df.training$RelapseStatus, df.validation$RelapseStatus))
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
  result.data <- data.frame(op.index, x, y, z, feature[feature.id], modus, 
                          all.mean.range.relapsed, all.var.range.relapsed, 
                          all.mean.range.relapsefree, all.var.range.relapsefree, 
                          df.active.t)
  names(result.data) <- c("Variable.ID", "X", "Y", "Z", "Feature", "Quadrant", 
    "Mean(Relapsed)", "Var(Relapsed)", "Mean(Relapsefree)", "Var(Relapsefree)", colnames(df.active.t))
  if (output.model) write.xlsx(result.data, variables.path) 

  ##########################################################################
  #create heatmap for all active coef.
  # df.training <- as.data.frame(df.training)
  # df.validation <- as.data.frame(df.validation)
  # df.training.active <- df.training[which(colnames(df.training) %in% names(df.training[, - c(1:3)])[op.index])]
  # df.validation.active <- df.validation[which(colnames(df.validation) %in% names(df.validation[, - c(1:3)])[op.index])]
  # df.active <- bind_rows(df.training.active, df.validation.active)
  # row.names(df.active) <- c(rownames.training, rownames.validation)
  # df.active <- bind_cols(as.data.frame(relapse), df.active)
  # row.names(df.active) <- c(rownames.training, rownames.validation)
  # df.active <- df.active[order(df.active$relapse), ]

  if (output.model) {
    ### for heatmap at least 2 features need to be present
    if (nrow(df.active.t) == 1) {
      df.active.t <- rbind(df.active.t, df.active.t)
    }
    pdf(heat.path, width = 12)
    heatmap(as.matrix(df.active.t),
    # heatmap(as.matrix(rbind(df.active.t, df.active.t)),
      scale = "none", Colv = NA, ColSideColors = relapse.status[order(relapse.status)],
      margins = c(5, 10), cexRow = 1.7, cexCol = 1.2)
    dev.off()
  }


  ########################################################################

  #create kaplan - meier - survival Kurve
  if (output.model) {
    time <- c(df.training[, 1], df.validation[, 1])
    status <- c(df.training[, 2], df.validation[, 2])
    df.new <- data.frame(time, status, comp.total[, 5])
    km.type = survfit(Surv(df.new[, 1], df.new[, 2]) ~ df.new[, 3], data = df.new, type = "kaplan-meier")
    ggsurvplot(km.type, conf.int = TRUE,
      legend.labs = c("low Risk", "high Risk"),
      ggtheme = theme_minimal(),
      pval = TRUE,
      pval.method = TRUE,
      risk.table = FALSE)
    ggsave(kaplan.path)
  }


  printf("AUC Value: %s", auc(roc(comp.total[, 2], comp.total[, 5])))
  printf("iAUC Value: %s", iAUC$iauc)

  sink(file = summary.path, append = TRUE)
  printf("AUC Value: %#.3f", AUC)
  printf("Err(train): %#.3f Err(val): %#.3f Err(total): %#.3f",
    1 - auc(roc(compred.training[, 2], compred.training[, 5])),
    1 - auc(roc(compred.validation[, 2], compred.validation[, 5])),
    1 - AUC)
  printf("iAUC Value: %s", iAUC$iauc)
  sink()

}
printf("Err(train): %#.3f Err(val): %#.3f Err(total): %#.3f",
    1 - auc(roc(compred.training[, 2], compred.training[, 5])),
    1 - auc(roc(compred.validation[, 2], compred.validation[, 5])),
    1 - AUC)
printf("### Done %s :: %s with on dataset=%s and subset=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id])