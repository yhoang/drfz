#!/usr/bin/R
# Authors: Eric Urbansky and Yen Hoang
# Transferred from 20200302

rm(list = ls())
options(max.print = 100)

### ---------- input
# working.station   FL, YH; different path settings
# conditional   TRUE/FALSE, if apply Ria's conditions
# project.id    1:5, for Basal, BCR, IL7, Pervanadate, TSLP
# sub.sample.name   func, func_plus3, func_plus6; different marker combinations
# set.alpha     0:1, set alpha which was found to have lowest error rate, see find_alpha.R
# sampling.size 1:10000, number of CV iterations
# cluster.size  1:12, find out maximum cluster size with detectCores()

working.station <- "YH"
conditional <- FALSE
project.id  <- 2
sub.sample.name <- "func"
set.alpha <- 1
sampling.size <- 100
cluster.size <- 3


# set path
if (working.station == "FL") {
    Project.path <- "/home/felix/AG_Baumgrass/Data/Good/Basal_for_Felix/Basal_for_Felix/"
    Rdata.path <- file.path("", "home", "felix", "AG_Baumgrass", "Scripts", "Pri", "Pri_good_established", "Rds")
    Text.path <- file.path("", "home", "felix", "AG_Baumgrass", "Data", "Good", "Marker_combinations")
} else if (working.station == "YH") {
    Project.path <- file.path("D:", "drfz", "Good2018")
    Rdata.path <- file.path("D:", "drfz", "Good2018", "Rdata")
    Text.path <- file.path("D:", "drfz", "Good2018", "tables")
}


### ---------- load R packages
librarys <- c("survival", "glmnet", "readxl", "dplyr", "doParallel", "xlsx", "survminer", "pROC")
#install.packages(librarys, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(librarys, require, character.only = TRUE)

#a useful print function
printf <- function(...) invisible(print(sprintf(...)))
project.name <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")

condi <- vector()
if (conditional) condi <- "condi1.1_"

### ---------- load and clean patient data from excel
patient_data <- read_excel(paste0(Text.path, "/patient_cohort.xlsx"))
#View(patient_data)

#set columstype from patien data str to numeric values
#colums with numeric value : Age at Diagnosis, WBC Count, Date of Diagnosis, Time to Relapse(Days), CCR (Days)
patient_data$`Time to Relapse (Days)` <- as.numeric(patient_data$`Time to Relapse (Days)`)
patient_data$`CCR (Days)` <- as.numeric(patient_data$`CCR (Days)`)
patient_data$`Time to Relapse (Days)`[is.na(patient_data$`Time to Relapse (Days)`)] <- 0
patient_data$`CCR (Days)`[is.na(patient_data$`CCR (Days)`)] <- 0
patient_data$`Age at Diagnosis` <- as.numeric(patient_data$`Age at Diagnosis`)
patient_data$`WBC Count` <- as.numeric(patient_data$`WBC Count`)
# patient_data$`Date of Diagnosis` <- as.numeric(patient_data$`Date of Diagnosis`)

#new Collum "Survival Time (Days)
patient_data$`Survival Time (Days)` <- patient_data$`Time to Relapse (Days)` +  patient_data$`CCR (Days)`
#fehlermeldung unknown column <- "Survival Time (Days) vorher erstellen?

#yes2 == yes
patient_data$`Relapse Status`[patient_data$`Relapse Status` == "Yes2"] <- "Yes"

#Set DDPR Status to nuermic binary values
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "Low"] <- 0
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "High"] <- 1

#reduce patien data to necessary colums `Patient ID`, `Relapse Status`,`MRD Risk`, `DDPR Risk`, Cohort, `Survival Time (Days)`
cohort <- patient_data[, c(1, 11, 8, 16, 15, 17)]

#divide cohort into training and validation set
training.set <- cohort[which(cohort$Cohort == "Training"), 1:6]
validation.set <- cohort[which(cohort$Cohort == "Validation"), 1:6]

#bind training and validation set to totalset
total.set <- bind_rows(training.set, validation.set)

### ---------- load PRI features as RDS
df.training <- readRDS(paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Training_", sub.sample.name, "_quadrant_absRange_", condi, "cof0.2.rds"))
df.validation <- readRDS(paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Validation_", sub.sample.name, "_quadrant_absRange_", condi, "cof0.2.rds"))


#safe rownames ("patien ID")
rownames.validation <- row.names(df.validation)
rownames.training <- row.names(df.training)

typeColNum <- 1

#relapsstatus in training set as numeric
condition.train <- training.set$`Relapse Status`[which(training.set$`Patient ID` %in% row.names(df.training))]
condition.train[condition.train %in% c("Yes")] <- 1
condition.train[condition.train %in% c("No")] <- 0
condition.train <- as.numeric(condition.train)

#survival time in training.set as numeric
survival.train <- training.set$`Survival Time (Days)`[which(training.set$`Patient ID` %in% row.names(df.training))]

#ddpr STATUS as numeric
ddpr.train <- training.set$`DDPR Risk`[which(training.set$`Patient ID` %in% row.names(df.training))]
ddpr.train <- as.numeric(ddpr.train)

#add condition/survivaltime/age to df.training
df.training <- bind_cols(as.data.frame(ddpr.train), df.training)
names(df.training)[typeColNum] <- "DDPR Status"
df.training <- bind_cols(as.data.frame(condition.train), df.training)
names(df.training)[typeColNum] <- "Relapse Status"
df.training <- bind_cols(as.data.frame(survival.train), df.training)
names(df.training)[typeColNum] <- "Survivaltime (Days)"
# glmnet need input matrix for model
df.training <- as.matrix(df.training)


#do same for validation set
condition.val <- validation.set$`Relapse Status`[which(validation.set$`Patient ID` %in% row.names(df.validation))]
condition.val[condition.val %in% c("Yes")] <- 1
condition.val[condition.val %in% c("No")] <- 0
condition.val <- as.numeric(condition.val)
survival.val <- validation.set$`Survival Time (Days)`[which(validation.set$`Patient ID` %in% row.names(df.validation))]
survival.val <- as.numeric(survival.val)
ddpr.val <- validation.set$`DDPR Risk`[which(validation.set$`Patient ID` %in% row.names(df.validation))]
ddpr.val <- as.numeric(ddpr.val)

df.validation <- bind_cols(as.data.frame(ddpr.val), df.validation)
names(df.validation)[typeColNum] <- "DDPR Status"
df.validation <- bind_cols(as.data.frame(condition.val), df.validation)
names(df.validation)[typeColNum] <- "Relapse Status"
df.validation <- bind_cols(as.data.frame(survival.val), df.validation)
names(df.validation)[typeColNum] <- "Survivaltime (Days)"
# glmnet need input matrix for model
df.validation <- as.matrix(df.validation)

### ------------ convert NaN/ +-Inf to -0.01 to not to confuse with 0 as PRI feature value
if (any(is.nan(df.training))) df.training[is.nan(df.training) | is.infinite(df.training)] <- -0.01
if (any(is.infinite(df.training))) df.training[is.infinite(df.training)] <- -0.01
if (any(is.nan(df.validation))) df.validation[is.nan(df.validation) | is.infinite(df.validation)] <- -0.01
if (any(is.infinite(df.validation))) df.validation[is.infinite(df.validation)] <- -0.01

### ------------ convert NAs to sample group mean
if (any(is.na(df.training))) {
  for (i in 2:ncol(df.training)) {
    NA.idx <- which(is.na(df.training[, i]))

    for (j in NA.idx) {
      tmp <- df.training[which(df.training[j, 1] == df.training[, 1]), i]
      tmp <- tmp[-which(is.na(tmp))]
      tmp <- round(sample(seq(mean(tmp) - sd(tmp), mean(tmp) + sd(tmp), by = 0.01), 1), 2)

      df.training[j, i] <- tmp
    }
  }
}

if (any(is.na(df.validation))) {
  for (i in 2:ncol(df.validation)) {
    NA.idx <- which(is.na(df.validation[, i]))

    for (j in NA.idx) {
      tmp <- df.validation[which(df.validation[j, 1] == df.validation[, 1]), i]
      tmp <- tmp[-which(is.na(tmp))]
      tmp <- round(sample(seq(mean(tmp) - sd(tmp), mean(tmp) + sd(tmp), by = 0.01), 1), 2)

      df.validation[j, i] <- tmp
    }
  }
}

### initialize and set seed for reproduction
seed.vec <- sample(sampling.size)
it.total <- 0
all.activ.Index <- all.coef.value <- vector()
variable.count <- vector()
min.cvm <- 100

# initiate cluster for parallelism
cl <- makeCluster(cluster.size)
registerDoParallel(cl)

### ------ cross validation
timeStart <- Sys.time()
printf("###Start %s.###", timeStart)

while (it.total < sampling.size) {
  it.total <- it.total + 1
  set.seed(seed.vec[it.total])

  ### stratify folds for cross validation because of imbalanced data
  #set folds that 1 fold contains at least 1 relapse
  # result: 4 folds with 8 patients (at least 2 relapse)
  if (conditional) {
    #fold.id for condition 1.1
    null.id <- one.id <- vector()
    fold.id <- rep(0, nrow(df.training))
    one.id <- rep(sample(1:4), 3)
    length(null.id)
    null.id <- rep(sample(1:4), 5)
    length(one.id)
    it.set <- it.null <- it.one <- 1
  } else {
    # #fold.id for no condidion
    null.id <- one.id <- vector()
    fold.id <- rep(0, nrow(df.training))
    one.id <- c(rep(sample(1:4), 3), 3)
    length(null.id)
    null.id <- rep(sample(1:4), 8)
    length(one.id)
    it.set <- it.null <- it.one <- 1
  }

  # create folds
  while (it.set < nrow(df.training) + 1) {
    if (df.training[it.set, 2] == 0) {
      fold.id[it.set] <- null.id[it.null]
      it.null <- it.null + 1
    } else {
      fold.id[it.set] <- one.id[it.one]
      it.one <- it.one + 1
    }
    it.set <- it.set + 1
  }
    
  # columns sampling
  tmp <- df.training[, - c(1:3)]
  tmp <- tmp[, sample(ncol(tmp))]

  #create Model with cross validation
  cv.fit <- cv.glmnet(tmp, Surv(df.training[, 1],
   df.training[, 2]), family = "cox",
   alpha = set.alpha,
   foldid = fold.id,
   parallel = TRUE)
  #plot(cv.fit)

  ### collect all active (>0) coefficients from fited model
  # lambda.min: model with min cross validation error
  Coefficients <- coef(cv.fit, s = cv.fit$lambda.min)
  Active.Index <- which(Coefficients != 0)
  coef.value <- coef(cv.fit, s = cv.fit$lambda.min)
  all.coef.value <- c(all.coef.value, coef.value)
  variable.count <- c(variable.count, length(Active.Index))
  print(length(Active.Index))
  tmp <- min(cv.fit$cvm)

  #collect best model (model with min cross validation error)
  #safe best model as op.fit
  #safe active coefficients
  if (min.cvm > tmp) {
    min.cvm <- tmp
    op.fit <- cv.fit
    op.variable.count <- length(Active.Index)
    op.index <- Active.Index

  }
  if (it.total %% 10 == 0) {
    printf("Run it #%s..", it.total)
  }

}
printf("End Process")
print(Sys.time() - timeStart)

stopCluster(cl)


print(table(variable.count))
printf("Best model with #features: %s", op.variable.count)

#prediction for training and validation set based on op.fit
#prediction for type cox model = relativ risk (RR)
p.training <- predict(op.fit, newx = df.training[, - c(1:3)], s = "lambda.min", type = "response")
p.validation <- predict(op.fit, newx = df.validation[, - c(1:3)], s = "lambda.min", type = "response")

#calculat threshold for cutoff
# "low Risk" < treshold > "high risk"
#scale RR between 0->1
scale.prediction <- (p.training - min(0)) / (max(p.training) - min(0))
scale.prediction[which(scale.prediction < 0)] <- 0

#performance test for set of thresholds
threshold <- seq(0, 1, 0.01)
predictions.roc <- data.frame()

#safe all prediction cutoffs as df
for (it in 1:length(threshold)) {
  newline <- findInterval(scale.prediction, threshold[it])
  predictions.roc <- bind_cols(as.data.frame(newline), predictions.roc)
  names(predictions.roc)[1] <- threshold[it]
}

#mirror df to set index right
predictions.roc <- predictions.roc[, ncol(predictions.roc):1]

#callculate error
#true prositiv = prediction 1 & relapsstatus 1
#calculate ssensitivity and false positiv rate
all.sens <- all.fp.rate <- auc <- vector()
for (i in 1:length(threshold)) {
  tp <- length(which(predictions.roc[, i] == 1 & df.training[, 2] == 1))
  fp <- length(which(predictions.roc[, i] == 1 & df.training[, 2] == 0))
  all.p <- length(which(df.training[, 2] == 1))
  all.n <- length(which(df.training[, 2] == 0))
  sens <- tp / all.p
  all.sens <- c(all.sens, sens)
  fp.rate <- fp / all.n
  all.fp.rate <- c(all.fp.rate, fp.rate)
  auc <- c(auc, auc(roc(df.training[, 2], predictions.roc[, i])))
}

#safe ROC plot as pdf
pdf(paste0(Project.path, "coxhazard/", project.name[project.id], "/HIST-ROC-Plot.pdf"))
par(mfrow = c(2, 2), pty = "s")
#plot Roc Kurve
plot(all.fp.rate, all.sens, type = "l", ylab = "Sensitivity", xlab = "False positiv Rate", ylim = c(0, 1), xlim = c(0, 1), main = "ROC Kurv")

###########################################################################
# #treshold by log-rank test. find "optimal" p-value for log rank
# p.value <- c()
# for(i in 1:length(names(predictions.roc))){
# df.new <- data.frame(df.training[, 1], df.training[, 2], predictions.roc[, i])
# names(df.new)[1] <- "TIME"
# names(df.new)[2] <- "Status"
# names(df.new)[3] <- "PREDICTION"
# km.type <- survfit(Surv(df.new$TIME, df.new$Status) ~ df.new$PREDICTION,
#   data = df.new,
#   type="kaplan-meier")
# tmp <- surv_pvalue(km.type, method = "1")$pval
# p.value <- c(p.value, tmp)
# }
# p.value[which(is.na(p.value) == TRUE)] <- 1
# op.thresh <- threshold[which(p.value == min(p.value))[1]]
##########################################################################

#finde best threshold
#fp = 0 & sens <=90
temp <- which(all.sens[which(all.fp.rate == 0)] <= 0.99)[1]

op.thresh <- which(all.sens == all.sens[which(all.fp.rate == 0)][temp])
op.thresh <- threshold[op.thresh][1]

printf("RR over %s are interpret as Status 1(Hight Risk)", op.thresh * max(p.training))
printf("RR under %s are Status 0(Low Risk)", op.thresh * max(p.training))


#pediction RR Values aprox Relapsstaus with calculates Threshold
predicted.validation <- findInterval(p.validation / max(p.training), op.thresh)
#print(predicted.validation)
predicted.training <- findInterval(p.training / max(p.training), op.thresh)
#print(predicted.training)

####################################################################

#find error in Prediction.validation for ddpr
fehler.ddpr <- df.validation[, 2] - df.validation[, 3]
fehler.ddpr[fehler.ddpr == 1] <- "FN"
fehler.ddpr[fehler.ddpr == -1] <- "FP"
fehler.ddpr[fehler.ddpr == 0] <- ""

#find error in Prediction model
fehler.model <- df.validation[, 2] - predicted.validation
fehler.model[fehler.model == 1] <- "FN"
fehler.model[fehler.model == -1] <- "FP"
fehler.model[fehler.model == 0] <- ""

#compare Real Status / ddpr Status / Model Status
vergleich.validation <- data.frame(rownames.validation, df.validation[, 2], df.validation[, 3], fehler.ddpr, predicted.validation, fehler.model)
names(vergleich.validation)[1] <- "Patien ID"
names(vergleich.validation)[2] <- "Real Status"
names(vergleich.validation)[3] <- "DDPR Status"
names(vergleich.validation)[4] <- "Fehler DDPR"
names(vergleich.validation)[5] <- "Predicted Status Model"
names(vergleich.validation)[6] <- "Error Model"
write.xlsx(vergleich.validation, paste0(Project.path, "coxhazard/", project.name[project.id], "/Prediction_DDPR_PRI_Training.xlsx"))

#find error in Prediction.training for ddpr
fehler.ddpr <- df.training[, 2] - df.training[, 3]
fehler.ddpr[fehler.ddpr == 1] <- "FN"
fehler.ddpr[fehler.ddpr == -1] <- "FP"
fehler.ddpr[fehler.ddpr == 0] <- ""

#find error in prediction model
fehler.model <- df.training[, 2] - predicted.training
fehler.model[fehler.model == 1] <- "FN"
fehler.model[fehler.model == -1] <- "FP"
fehler.model[fehler.model == 0] <- ""

#compare Real Status / ddpr Status / Model Status
vergleich.training <- data.frame(rownames.training, df.training[, 2], df.training[, 3], fehler.ddpr, predicted.training, fehler.model)
names(vergleich.training)[1] <- "Patient ID"
names(vergleich.training)[2] <- "Real Status"
names(vergleich.training)[3] <- "DDPR Status"
names(vergleich.training)[4] <- "Error DDPR"
names(vergleich.training)[5] <- "Predicted Status Model"
names(vergleich.training)[6] <- "Error Model"
write.xlsx(vergleich.training, paste0(Project.path, "coxhazard/", project.name[project.id], "/Prediction_DDPR_PRI_Validation.xlsx"))

#bind training and validation error
vergleich.total <- bind_rows(vergleich.training, vergleich.validation)
write.xlsx(vergleich.total, paste0(Project.path, "coxhazard/", project.name[project.id], "/Prediction_DDPR_PRI_Total.xlsx"))

##############################################################################
#auc
print(auc(roc(vergleich.total[, 2], vergleich.total[, 5])))
plot(roc(vergleich.total[, 2], vergleich.total[, 5]))

##############################################################################
#get active coeffizienten !!
Coef.names <- colnames(df.training[, - c(1:3)])

#get names from coef.names for all active index in op.fit
all.activ.names <- colnames(df.training[, - c(1:3)])[op.index]
all.activ.names.split <- strsplit(all.activ.names[1:length(all.activ.names)], ".", fixed = TRUE)

x <- vector()
y <- vector()
z <- vector()
modus <- vector()
quat <- vector()
i <- 0
while (i < length(all.activ.names)) {
  i <- i + 1
  x <- c(x, all.activ.names.split[[i]][1])
  y <- c(y, all.activ.names.split[[i]][2])
  z <- c(z, all.activ.names.split[[i]][3])
  modus <- c(modus, all.activ.names.split[[i]][4])
  quat <- c(quat, all.activ.names.split[[i]][5])
}

dev.off()

#safe coef names and values for op.fit as excel
df.training <- as.data.frame(df.training)
df.validation <- as.data.frame(df.validation)

#create heatmap for op.fit coef
df.head.training <- df.training[which(colnames(df.training) %in% all.activ.names)]
df.head.validation <- df.validation[which(colnames(df.validation) %in% all.activ.names)]
df.head <- bind_rows(df.head.training, df.head.validation)
row.names(df.head) <- c(rownames.training, rownames.validation)

#Relapsstatus as Vector Red/blue for Heatmap
relapse <- as.numeric(c(df.training$`Relapse Status`, df.validation$`Relapse Status`))
relapse[relapse == 1] <- "red"
relapse[relapse == 0] <- "blue"

#calculate mean and var for training and validation
all.relapse <- cohort$`Patient ID`[which(cohort$`Relapse Status` == "Yes")]
all.relapsfree <- cohort$`Patient ID`[which(cohort$`Relapse Status` == "No")]
all.mean.range.relapse <- all.var.range.relapse <- all.mean.range.relapsfree <- all.var.range.relapsfree <- vector()

it <- 0
for (i in 1:length(names(df.head))) {
  it <- it + 1
  all.mean.range.relapse <- c(all.mean.range.relapse, mean(df.head[which(row.names(df.head) %in% all.relapse), it]))
  all.var.range.relapse <- c(all.var.range.relapse, var(df.head[which(row.names(df.head) %in% all.relapse), it]))
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
    all.mean.range.relapse, all.var.range.relapse,
    all.mean.range.relapsfree, all.var.range.relapsfree,
    new.head)
names(result.data)[10] <- "Var(AbsRange) Relapsfree"
names(result.data)[9] <- "Mean(AbsRange) Relapsfree"
names(result.data)[8] <- "Var(AbsRange) Relaps"
names(result.data)[7] <- "Mean(AbsRange) Relaps"
names(result.data)[6] <- "Quadrant"
names(result.data)[5] <- "Modus"
names(result.data)[4] <- "z Variable"
names(result.data)[3] <- "y Variable"
names(result.data)[2] <- "x Variable"
names(result.data)[1] <- "Variable Index"

write.xlsx(result.data, paste0(Project.path, "coxhazard/", project.name[project.id], "/Prediction_Variables.xlsx"))

##########################################################################

#create heatmap for all active coef.
df.training <- as.data.frame(df.training)
df.validation <- as.data.frame(df.validation)
df.head.training <- df.training[which(colnames(df.training) %in% names(df.training[, - c(1:3)])[op.index])]
df.head.validation <- df.validation[which(colnames(df.validation) %in% names(df.validation[, - c(1:3)])[op.index])]
df.head <- bind_rows(df.head.training, df.head.validation)
row.names(df.head) <- c(rownames.training, rownames.validation)
df.head <- bind_cols(as.data.frame(relapse), df.head)
row.names(df.head) <- c(rownames.training, rownames.validation)
df.head <- df.head[order(df.head$relapse), ]

pdf(paste0(Project.path, "coxhazard/", project.name[project.id], "/Heatmap.pdf"), width = 10)
heatmap(t(as.matrix(df.head[, -1])), scale = "none", Colv = NA, ColSideColors = relapse[order(relapse)])
dev.off()

########################################################################

#create kaplan-meier - survival Kurve
time <- c(df.training[, 1], df.validation[, 1])
status <- c(df.training[, 2], df.validation[, 2])
df.new <- data.frame(time, status, vergleich.total[, 5])
km.type <- survfit(Surv(df.new[, 1], df.new[, 2]) ~ df.new[, 3], data = df.new, type = "kaplan-meier")
ggsurvplot(km.type, conf.int = TRUE, legend.labs = c("low Risk", "high Risk"), ggtheme = theme_minimal(), pval = TRUE, pval.method = TRUE)

