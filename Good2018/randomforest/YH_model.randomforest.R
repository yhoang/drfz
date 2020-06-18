#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

add.libraries <- c("randomForest")
lapply(add.libraries, require, character.only = TRUE)

print("Running random forest model..")
ntree.val = 10000
mtry.val = 1000

comment.out <- paste(comment.out, paste0("ntree", ntree.val), paste0("mtry", mtry.val), sep=".")

#### variables xlsx
variables.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(7))
# VarImportance PDF
varimp.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(13))
# summary txt
summary.path <- sprintf("%s/%s.log", Output.path, output_filenaming(11))

### ------------------ random forest model 
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
  
  ### Create a Random Forest model with default parameters
  # We can tune the random forest model by changing the number of trees (ntree) and the number of variables randomly sampled at each stage (mtry). According to Random Forest package description:
  # ntree: Number of trees to grow. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times.
  # mtry: Number of variables randomly sampled as candidates at each split. Note that the default values are different for classification (sqrt(p) where p is number of variables in x) and regression (p/3) 
  # default: ntree=500, mtry=2
  # options(expressions = 5e5)
  # memory.limit(size = 5000000)
  ### in command
  # cd "C:\Program Files\R\R-3.6.1\bin"
  #./R.exe  --max-ppsize 500000
  #source("YH_run.R")
  rf <- randomForest(x = df.training[, - c(1:3)], y = as.factor(df.training[, 2]), 
    ntree = ntree.val, mtry = mtry.val, importance = TRUE, do.trace = 1000, proximity = TRUE)

if (FALSE) {
trctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  allowParallel = TRUE)

tgctrl <- expand.grid(
    .mtry = 7,
    .splitrule = "gini",
    .min.node.size = 10)

fit.rf <- train(Survived ~ 
                    Sex + 
                    title.class + 
                    family.size.simple + 
                    Fare.grps + 
                    Pclass + 
                    Age.grps + 
                    Embarked + 
                    ticket.prefix + 
                    race,
             data = train, 
             trControl = trctrl,
             metric = "Accuracy",
             importance = "impurity",
             tuneGrid = tgctrl,
             num.trees = 2000,
             method = "ranger")
}
  # rf3 <- randomForest(x = df.training[, - c(1:3)], y = as.factor(df.training[, 2]), 
  #   xtest = df.validation[, - c(1:3)], ytest = as.factor(df.validation[, 2]), 
  #   ntree = 1000, importance = TRUE, do.trace = 100, proximity = TRUE)

  if (output.model) {
    pdf(varimp.path, width=20)
    varImpPlot(rf, n.var = 20, main = "Top 20 Variable Importance")
    dev.off()
  }

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
  

  pred.training <- predict(rf, newdata = df.training[, -c(1:3)])
  pred.validation <- predict(rf, newdata = df.validation[, -c(1:3)])

  # show prediction vs. reality
  table(pred.training, df.training[, 2])
  table(pred.validation, df.validation[, 2])
  MSE <- mean(((as.numeric(pred.validation) - 1) - df.validation[, 2]) ^ 2)

  # plot to see the margin, positive or negative, if positif it means correct classification
  # pdf("test.pdf")
  #   plot(margin(rf, as.factor(df.training[, 2])))
  # dev.off()

  # pdf("MDS_rf.pdf")
  #   MDSplot(rf, df.training[, 2], k=2, palette = 2, pch=df.training[, 2])
  # dev.off()
      
  # pdf("rf.pdf")
  #   plot(rf)
  #   # red = 0, green = 1
  # dev.off()

  # determine the misclassification rate. First, build a confusion matrix. Each column of the matrix represents the number of predictions of each class, while each row represents the instances in the actual class.
  CM.training = table(pred.training, df.training[, 2])
  CM.validation = table(pred.validation, df.validation[, 2])
  
  # Second, build a diagonal mark quality prediction. Applying the diag function to this table then selects the diagonal elements, i.e., the number of points where random forest agrees with the true classification, and the sum command simply adds these values up.
  accuracy.train <- (sum(diag(CM.training))) / sum(CM.training)
  accuracy.val <- (sum(diag(CM.validation))) / sum(CM.validation)
  accuracy.total <- (sum(diag(CM.training) + diag(CM.validation))) / (sum(CM.training) + sum(CM.validation))
  
  printf("AUC(train): %s AUC(val): %s AUC(total): %s", accuracy.train, accuracy.val, accuracy.total)
  printf("Err(train): %s Err(val): %s Err(total): %s", 1 - accuracy.train, 1 - accuracy.val, 1 - accuracy.total)
  printf("MSE.val=%s", MSE)
  if (it.total %% 10 == 0) {
    printf("At Work IT:%s ", it.total)
    print(Sys.time() - timeStart)
    print(proc.time() - ptm)
  }
  
# Create Confusion Matrix
print(confusionMatrix(data = pred.training,  
                reference = df.training[, 2],
                positive = 1))
  
  
  sink(file = summary.path, append = TRUE)
    printf("### Start %s with %s on dataset=%s and subset=%s lambda=%s random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], lambda[lambda.id], randomize.label, spikeIns)
    printf("Sampling size: %s with seed ID: %s", sampling.size, seed.vec[it.total])
    printf("Random sampling of condition: %s", randomize.label)
    printf("Introduce Spike Ins: %s", spikeIns)
    printf("AUC(train): %#.3f AUC(val): %#.3f AUC(total): %#.3f", accuracy.train, accuracy.val, accuracy.total)
    printf("Err(train): %#.3f Err(val): %#.3f Err(total): %#.3f", 1 - accuracy.train, 1 - accuracy.val, 1 - accuracy.total)
  sink()

  pred.training <- predict(rf, newdata = df.training[, -c(1:3)], type ="probe")
  pred.validation <- predict(rf, newdata = df.validation[, -c(1:3)], type ="probe")

}
printf("Done with %s.", model[model.id])
print(Sys.time() - timeStart)
