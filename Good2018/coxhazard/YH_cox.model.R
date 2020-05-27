#!/usr/bin/R
# Authors: Eric Urbansky, Felix Lohrke and Yen Hoang

### - - - - - - - - - -  input
# working.station   FL, YH; different path settings
# conditional   TRUE/FALSE, if apply Ria's conditions
# dataset.id    1:5, for Basal, BCR, IL7, Pervanadate, TSLP
# subgroup      "Training", "Validation", "Total"
# subset.id     1:4 for func, func3, func6; different marker combinations
# subject.in.id 1:13 for "df", "quadrant", "ROC", "AUC", "iAUC", "CompV", "CompT", "CompTotal", "Variables", "Heat", "CompAll", "KM", "Thresholds"
# set.alpha     0:1, set alpha which was found to have lowest error rate, see find_alpha.R
# sampling.size 1:10000, number of CV iterations
# cluster.size  1:12, find out maximum cluster size with detectCores()
working.station <- "YH"
initdate <- "YH200526"
dataset.id <- 1
subset.id <- 4
subject.in.id <- 2
feature.id <- 1
comment.in <- "autoSec.cof0.2"
set.alpha <- 1
sampling.size <- 20
cofactor <- 0.2
cluster.size <- 3
conditional <- FALSE

### Parameters 
## FP == value
FP.value <- 0.1
### Sensitivity  <=  value
SENS.value <- 1
## threshold to select (should multiple apply to parameters)
selected.thresh <- 1
### Output activated?
## generate output that doesnt change in the model (collection of ROC values, ROC curve, variables, heat map)
## generally only apply once
output.model <- FALSE
## generate output for different selected thresholds (threshold selected + info, AUC)
output.thresholds <- FALSE
comment.out <- sprintf("%s.FP%s.SENS%s.Thresh%s", comment.in, FP.value, SENS.value, selected.thresh)
#############################################################






### - - - - - - - - - -  load R packages and scripts
librarys <- c("survival", "glmnet", "readxl", "dplyr", "doParallel", "xlsx", "survminer", "pROC", "stringr")
#install.packages(librarys, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(librarys, require, character.only = TRUE)
#printfunktion
#a useful print function
printf <- function(...) invisible(print(sprintf(...)))
source(file.path("D:", "drfz", "Good2018", "coxhazard", "YH_func_naming.R"))
source(file.path("D:", "drfz", "Good2018", "coxhazard", "YH_load_patientdata.R"))


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
    Output.path <- file.path(Project.path, "coxhazard", dataset[dataset.id], subset[subset.id])
}


######## INPUT FILES #####################
# patient data 
patient.data.path <- file.path(Text.path, "patient_cohort.xlsx")
# subgroup <- c("Training", "Validation", "Total")
training.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(1))
validation.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(2))
# total data as rds (if it does not exist a new will one be created at path location)
total.data.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(3))
##########################


######## OUTPUT FILES #####################
### create folder first
dir.create(Output.path, showWarnings = FALSE)
### subject = c("df", "quadrant", "ROC", "AUC", "iAUC", "Comp.Val", "Comp.Train", "Comp.Total", "Variables", "Heat", "Comp.Error", "KM", "Thresholds")
# ROC / AUC Curve pdf
roc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(3))
auc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(4))
iauc.pdf.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(5))
# comparisons of validation, training and total xlsx
comp.val.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(6))
comp.train.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(7))
comp.total.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(8))
#### variables xlsx
variables.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(9))
# heatmap pdf
heat.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(10))
# comparisons of Sensitivity, specificity and FP-rate to associated thresholds
comp.error.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(11))
# kaplan-meier curve
kaplan.path <- sprintf("%s/%s.pdf", Output.path, output_filenaming(12))
# thresholds txt
thresholds.path <- sprintf("%s/%s.xlsx", Output.path, output_filenaming(13))
################################


#safe triplot rds file as df. for training and validation set
#absRange without condition
printf("Loading RDS")
df.training <- readRDS(training.data.path)
df.validation <- readRDS(validation.data.path)
printf("RDS loaded")

printf("Dimensions of training and validation set: ")
print(dim(df.training))
print(dim(df.validation))

#safe rownames ("Patient ID")
rownames.training <- row.names(df.training)
rownames.validation <- row.names(df.validation)

#### Check if df.total already exists as RDS (binding is bottleneck)
#### RDS with df.total will be created once if it doesnt exist and then read every run
if (file.exists(total.data.path)){
  df.total <- readRDS(total.data.path)
  printf("df.total does exist...using created RDS")
} else {
  printf("df.total does not exist...creating RDS")
  df.total <- bind_rows(df.training, df.validation)
  # rownames(df.total) <- c(rownames(df.training), rownames(df.validation))
  saveRDS(df.total, file=total.data.path)
}
rownames(df.total)=c(rownames.training, rownames.validation)

printf("Num of patients:")
print(nrow(df.total))
#View(df.total)
#total 60 patients in cohort
#Basal no condition total 54 patients
#Basal condition 1.1 total 41 patients
#BCR 45 patients total, divided in 37 in training and 8 in validation
typeColNum <- 1

#relapsstatus in training set as numeric
condition <- training.set$`Relapse Status`[which(training.set$`Patient ID` %in%  row.names(df.training))]
condition[condition %in% c("Yes")] <- 1
condition[condition %in% c("No")] <- 0
condition=as.numeric(condition)

#survival time in training.set as numeric
survival=training.set$`Survival Time (Day)`[which(training.set$`Patient ID` %in% row.names(df.training))]
survival=as.numeric(survival)

#DDPR STATUS as numeric
DDPR=training.set$`DDPR Risk`[which(training.set$`Patient ID` %in% row.names(df.training))]
DDPR=as.numeric(DDPR)

#add condition/survivaltime/age to df.training
df.training=bind_cols(as.data.frame(DDPR), df.training)
names(df.training)[typeColNum] <- "DDPR Status"
df.training=bind_cols(as.data.frame(condition), df.training)
names(df.training)[typeColNum]="Relapse Status"
df.training=bind_cols(as.data.frame(survival), df.training)
names(df.training)[typeColNum] <- "Survivaltime (Day)"



#do same for validation set
condition <- validation.set$`Relapse Status`[which(validation.set$`Patient ID` %in%  row.names(df.validation))]
condition[condition %in% c("Yes")]=1
condition[condition %in% c("No")]=0
condition=as.numeric(condition)
survival <- validation.set$`Survival Time (Day)`[which(validation.set$`Patient ID` %in% row.names(df.validation))]
survival <- as.numeric(survival)
DDPR=validation.set$`DDPR Risk`[which(validation.set$`Patient ID` %in% row.names(df.validation))]
DDPR=as.numeric(DDPR)

df.validation=bind_cols(as.data.frame(DDPR), df.validation)
names(df.validation)[typeColNum] <- "DDPR Status"
df.validation=bind_cols(as.data.frame(condition), df.validation)
names(df.validation)[typeColNum]="Relapse Status"
df.validation=bind_cols(as.data.frame(survival), df.validation)
names(df.validation)[typeColNum] <- "Survivaltime (Day)"

df.training.strat <- df.training
rownames(df.training.strat) <- rownames.training

# glmnet need input matrix for model 
# safe df as matrix
df.training=as.matrix(df.training)
df.validation=as.matrix(df.validation)


### convert NaN/+ - Inf to -0.01 in df.training
if (any(is.nan(df.training))) df.training[is.nan(df.training) | is.infinite(df.training)] <- -0.01
if (any(is.infinite(df.training))) df.training[is.infinite(df.training)] <- -0.01
### convert NAs to sample group mean
if (any(is.na(df.training))) {
  for (i in 2:ncol(df.training)) {
    NA.idx <- which(is.na(df.training[, i]))
    
    for (j in NA.idx) {
      tmp <- df.training[which(df.training[j, 1] == df.training[, 1]), i]
      tmp <- tmp[-which(is.na(tmp))]
      tmp <- round(sample(seq(mean(tmp) - sd(tmp), mean(tmp) + sd(tmp), by=0.01), 1), 2)
      
      df.training[j, i] <- tmp
    }
  }
}
### convert NaN/+ - Inf to -0.01 in df.validation
if (any(is.nan(df.validation))) df.validation[is.nan(df.validation) | is.infinite(df.validation)] <- -0.01
if (any(is.infinite(df.validation))) df.validation[is.infinite(df.validation)] <- -0.01
### convert NAs to sample group mean
if (any(is.na(df.validation))) {
  for (i in 2:ncol(df.validation)) {
    NA.idx <- which(is.na(df.validation[, i]))
    
    for (j in NA.idx) {
      tmp <- df.validation[which(df.validation[j, 1] == df.validation[, 1]), i]
      tmp <- tmp[-which(is.na(tmp))]
      tmp <- round(sample(seq(mean(tmp) - sd(tmp), mean(tmp) + sd(tmp), by=0.01), 1), 2)
      
      df.validation[j, i] <- tmp
    }
  }
}

### initiate variables
# set seed for reproduction 
seed.vec <- sample(sampling.size)
it.total <- 0
cluster.size <- 1
all.activ.Index <- all.coef.value <- vector()
variabl.count <- vector()
min.cvm <- 100


cl <- makeCluster(cluster.size)
registerDoParallel(cl)

timeStart <- Sys.time()
ptm <- proc.time()
printf("### Start %s on dataset=%s and subset=%s ###", timeStart, dataset[dataset.id], subset[subset.id])

while (it.total < sampling.size) {
  it.total <- it.total + 1 

  set.seed(seed.vec[it.total])

  ### create fold values
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

  #create Model with cross validation
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
  Coefficients <- coef(cv.fit, s=cv.fit$lambda.1se)
  Active.Index <- which(Coefficients != 0)
  coef.value <- coef(cv.fit, s=cv.fit$lambda.1se)
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
printf("End Process")
print(Sys.time() - timeStart)

stopCluster(cl)


print(table(variabl.count))
printf("Best model >0: %s features", op.variabl.count)
printf("Cross validation error: %s", min.cvm)

#prediction for training and validation set based on op.fit
#prediction for type cox model = relativ risk (RR)
p.training <- predict(op.fit, newx = df.training[, -c(1:3)], s="lambda.1se", type="response")
p.validation <- predict(op.fit, newx = df.validation[, -c(1:3)], s="lambda.1se", type="response")

#calculat threshold for cutoff
# "low Risk" < treshold > "high risk"
#scale RR between 0 -> 1
scale.prediction <- (p.training - min(0)) / (max(p.training) - min(0))
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

#callculate error
#true prositiv = prediction 1 & relapsstatus 1
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

#safe  ROC plot as pdf
if (output.model == TRUE){
  pdf(roc.pdf.path)
  par(mgp = c(2, 0.5, 0), pty = "s", mar = c(4, 4, 3, 1), cex.lab=1.2)
  #plot Roc Kurve
  plot(all.FP.rate, all.sens, type = "l", ylab = "TP Rate (Sensitivity)", 
  xlab = "FP Rate", ylim = c(0, 1), xlim = c(0, 1), main = "ROC Curve")
  dev.off()
}

###########################################################################
# #treshold by log - rang test. find "optimal" p - value for log rank
p.value <- c()
for (i in 1:length(names(predictions.roc))){
  df.new <- data.frame(df.training[, 1], df.training[, 2], predictions.roc[, i])
  names(df.new)[1] <- "TIME"
  names(df.new)[2] <- "Status"
  names(df.new)[3] <- "PREDICTION"
  km.type=survfit(Surv(df.new$TIME, df.new$Status) ~ df.new$PREDICTION, 
                  data = df.new, 
                  type="kaplan-meier")
  tmp <- surv_pvalue(km.type, method = "1")$pval
  p.value <- c(p.value, tmp)
}
p.value[which(is.na(p.value) ==  TRUE)] <- 1
op.thresh <- threshold[which(p.value == min(p.value))[1]]
##########################################################################


# saving all FP - rates, SENS, SPEC and associated thresholds
if (output.model == TRUE){
  comp.all <- data.frame(data.frame(threshold * max(p.training)), data.frame(all.sens), data.frame(all.spec), data.frame(all.FP.rate))
  names(comp.all)[1] <- "threshold"
  names(comp.all)[2] <- "SENS"
  names(comp.all)[3] <- "SPEC"
  names(comp.all)[4] <- "FP-rate"
  write.xlsx(comp.all, comp.error.path)
}

### Eric's interpretation
#finde best threshold
# temp <- which(all.sens[which(all.FP.rate == FP.value)] <= SENS.value)[1]
# op.thresh <- which(all.sens == all.sens[which(all.FP.rate == FP.value)][temp])
# op.thresholds <- op.thresh

#sens = 1 and FP lowest
sens.idx <- which(all.sens == SENS.value)
## get lowest FP
FP.min.idx <- which(all.FP.rate == min(all.FP.rate[sens.idx]))
op.thresholds <- FP.min.idx
printf("Positions that have been selected as possible thresholds: %s", op.thresholds)

## selected threshold via position selected.thresh
op.thresh <- threshold[op.thresholds][selected.thresh]

## saved all found thresholds in op.thresholds
op.thresholds <- threshold[op.thresholds]

printf("Varied threshold Values: ")
print(op.thresholds)
printf("RR over %s are interpret as Status 1(Hight Risk) and under as Status 0(Low Risk)", op.thresholds * max(p.training))
printf("Threshold that has been selected %s at position %s", op.thresholds[selected.thresh], selected.thresh)

## saving threshold data in txt
if (output.thresholds == TRUE){
  sink(file <- thresholds.path)
  printf("Selected Parameters for FP rate <= %s and for SENS  ==  %s", FP.value, SENS.value)
  printf("Possible thresholds for this configuration: %s", op.thresholds)
  printf("Specifically selected Threshold for this configuration: %s at position: %s", op.thresholds[selected.thresh], selected.thresh)
  printf("RR for this threshold: %s", op.thresh * max(p.training))
  sink()
}

#pediction RR Values aprox Relapsstaus with calculates Threshold
predicted.validation=findInterval(p.validation / max(p.training), op.thresh)
#print(predicted.validation)
predicted.training=findInterval(p.training / max(p.training), op.thresh)
#print(predicted.training)

####################################################################

#find error in Prediction.validation for DDPR
fehler.ddpr <- df.validation[, 2] - df.validation[, 3]
fehler.ddpr[fehler.ddpr == 1]="FN"
fehler.ddpr[fehler.ddpr ==  - 1]="FP"
fehler.ddpr[fehler.ddpr == 0]=""

#find error in Prediction model
fehler.model <- df.validation[, 2] - predicted.validation
fehler.model[fehler.model == 1]="FN"
fehler.model[fehler.model ==  - 1]="FP"
fehler.model[fehler.model == 0]=""

#compare Real Status / DDPR Status / Model Status 
vergleich.validation <- data.frame(rownames.validation, df.validation[, 2], df.validation[, 3], fehler.ddpr, predicted.validation, fehler.model)
names(vergleich.validation)[1] <- "Patient ID"
names(vergleich.validation)[2] <- "Real Status"
names(vergleich.validation)[3] <- "DDPR Status"
names(vergleich.validation)[4] <- "Fehler DDPR"
names(vergleich.validation)[5] <- "Predicted Status Model"
names(vergleich.validation)[6] <- "Fehler Model"
if (output.thresholds == TRUE){
  write.xlsx(vergleich.validation, comp.val.path)
}

#find error in Prediction.training for DDPR
fehler.ddpr <- df.training[, 2] - df.training[, 3]
fehler.ddpr[fehler.ddpr == 1]="FN"
fehler.ddpr[fehler.ddpr ==  - 1]="FP"
fehler.ddpr[fehler.ddpr == 0]=""

#find error in prediction model
fehler.model <- df.training[, 2] - predicted.training
fehler.model[fehler.model == 1]="FN"
fehler.model[fehler.model ==  - 1]="FP"
fehler.model[fehler.model == 0]=""

#compare Real Status / DDPR Status / Model Status
vergleich.training <- data.frame(rownames.training, df.training[, 2], df.training[, 3], fehler.ddpr, predicted.training, fehler.model)
names(vergleich.training)[1] <- "Patient ID"
names(vergleich.training)[2] <- "Real Status"
names(vergleich.training)[3] <- "DDPR Status"
names(vergleich.training)[4] <- "Fehler DDPR"
names(vergleich.training)[5] <- "Predicted Status Model"
names(vergleich.training)[6] <- "Fehler Model"
if (output.thresholds == TRUE){
  write.xlsx(vergleich.training, comp.train.path)
}

#bind training and validation error
vergleich.total <- bind_rows(vergleich.training, vergleich.validation)
if (output.thresholds == TRUE){
  write.xlsx(vergleich.total, comp_total.path)
}

## saving threshold data in txt
if (output.thresholds == TRUE){
  sink(file = thresholds_path)
  printf("Selected Parameters for SENS >= %s and alpha = %s", FP.value, SENS.value, alpha.value)
  printf("Possible thresholds for this configuration: %s", op.thresholds)
  printf("Specifically selected Threshold for this configuration: %s at position: %s", op.thresholds[selected.thresh] * max(p.training), selected.thresh)
  printf("RR for this threshold: %s", op.thresh * max(p.training))
  printf("AUC Value: %s", auc(roc(vergleich.total[, 2], vergleich.total[, 5])))
  sink()
}
##############################################################################
#AUC
if (output.thresholds == TRUE){
  pdf(auc.pdf.path)
  plot(roc(vergleich.total[, 2], vergleich.total[, 5]))
  dev.off()
}
print(auc(roc(vergleich.total[, 2], vergleich.total[, 5])))

dev.off()

##############################################################################
#get active coeffizienten !!
Coef.names <- colnames(df.training[,  - c(1:3)])

#get names from coef.names for all active index in op.fit
all.activ.names <- colnames(df.training[,  - c(1:3)])[op.index]
all.activ.names.split <- strsplit(all.activ.names[1:length(all.activ.names)], ".", fixed = TRUE)

x <- vector()
y <- vector()
z <- vector()
modus <- vector()
quat <- vector()
i=0
while (i < length(all.activ.names)){
  i=i + 1
  x= c(x, all.activ.names.split[[i]][1])
  y= c(y, all.activ.names.split[[i]][2])
  z= c(z, all.activ.names.split[[i]][3])
  modus= c(modus, all.activ.names.split[[i]][4])
  quat <- c(quat, all.activ.names.split[[i]][5])
}

########################################################################
# create iAUC

if (output.thresholds == TRUE){
  pdf(iauc.pdf.path)
  iAUC = AUC.uno(Surv(df.training[, 1], df.training[, 2]), Surv(df.validation[, 1], df.validation[, 2]), p.validation, times=seq(10, 1000, 10),  savesensspec=TRUE)
  names(iAUC)
  iAUC$iauc
  plot(iAUC)
  dev.off()
}



#safe coef names and values for op.fit as excel
df.training <- as.data.frame(df.training)
df.validation <- as.data.frame(df.validation)

df.head.training <- df.training[which(colnames(df.training) %in% all.activ.names)]
df.head.validation <- df.validation[which(colnames(df.validation) %in% all.activ.names)]
df.head <- bind_rows(df.head.training, df.head.validation)
row.names(df.head) <- c(rownames.training, rownames.validation)

#Relapsstatus as Vector Red / blue for Heatmap
relaps <- as.numeric(c(df.training$`Relaps Status`, df.validation$`Relaps Status`))
relaps[relaps == 1] <- "red"
relaps[relaps == 0] <- "blue"

#calculate mean and var for training and validation 
all.relaps <- cohort$`Patient ID`[which(cohort$`Relapse Status` ==  "Yes")]
all.relapsfree <- cohort$`Patient ID`[which(cohort$`Relapse Status` ==  "No")]
all.mean.range.relaps <- all.var.range.relaps <- all.mean.range.relapsfree <- all.var.range.relapsfree <- vector()

it <- 0
for (i in 1:length(names(df.head))){
  it <- it  + 1
  all.mean.range.relaps <- c(all.mean.range.relaps, mean(df.head[which(row.names(df.head) %in% all.relaps), it]))
  all.var.range.relaps <- c(all.var.range.relaps, var(df.head[which(row.names(df.head) %in% all.relaps), it]))
  all.mean.range.relapsfree <- c(all.mean.range.relapsfree, mean(df.head[which(row.names(df.head) %in% all.relapsfree), it]))
  all.var.range.relapsfree <- c(all.var.range.relapsfree, var(df.head[which(row.names(df.head) %in% all.relapsfree), it]))
}

#safe all real absRange values 
#transponse df.head 
new.head <- t(df.head)
colnames(new.head) <- row.names(df.head)
row.names(new.head) <- c(1:dim(new.head)[1])

#safe df as excel
all.right.index <- op.index
result.data <- data.frame(all.right.index, x, y, z, modus, quat, 
                         all.mean.range.relaps, all.var.range.relaps, 
                         all.mean.range.relapsfree, all.var.range.relapsfree, 
                         new.head)
names(result.data)[10] <- "Var(AbsRange) Relapsefree"
names(result.data)[9] <- "Mean(AbsRange) Relapsefree"
names(result.data)[8] <- "Var(AbsRange) Relapsed"
names(result.data)[7] <- "Mean(AbsRange) Relapsed"
names(result.data)[6] <- "Quadrant"
names(result.data)[5] <- "Modus"
names(result.data)[4] <- "Z Variable"
names(result.data)[3] <- "Y Variable"
names(result.data)[2] <- "X Variable"
names(result.data)[1] <- "Variable Index"
if (output.model == TRUE){
  write.xlsx(result.data, variables.path) 
}


##########################################################################
#create heatmap for all active coef.
df.training <- as.data.frame(df.training)
df.validation <- as.data.frame(df.validation)
df.head.training <- df.training[which(colnames(df.training) %in% names(df.training[,  - c(1:3)])[op.index])]
df.head.validation <- df.validation[which(colnames(df.validation) %in% names(df.validation[,  - c(1:3)])[op.index])]
df.head <- bind_rows(df.head.training, df.head.validation)
row.names(df.head) <- c(rownames.training, rownames.validation)
df.head <- bind_cols(as.data.frame(relaps), df.head)
row.names(df.head) <- c(rownames.training, rownames.validation)
df.head <- df.head[order(df.head$relaps), ]

if (output.model == TRUE){
  pdf(heat.path, width = 12)

  heatmap(t(as.matrix(df.head[, - 1])),
    scale = "none", Colv = NA, ColSideColors = relaps[order(relaps)],
    margins = c(5, 15), cexRow = 1.7, cexCol = 1.1)
  dev.off()
}


########################################################################

#create kaplan - meier - survival Kurve
if (output.model == TRUE){
  time <- c(df.training[, 1], df.validation[, 1])
  status <- c(df.training[, 2], df.validation[, 2])
  df.new <- data.frame(time, status, vergleich.total[, 5])
  km.type=survfit(Surv(df.new[, 1], df.new[, 2]) ~ df.new[, 3], data = df.new, type="kaplan-meier")
  ggsurvplot(km.type, conf.int = TRUE, legend.labs = c("low Risk", "high Risk"), ggtheme = theme_minimal(), pval = TRUE, pval.method = TRUE)
  ggsave(kaplan.path)
}

printf("AUC Value: %s", auc(roc(vergleich.total[, 2], vergleich.total[, 5])))