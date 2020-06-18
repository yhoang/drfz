#!/usr/bin/R

### ------------------- load PRI features and patient meta data and preprocess
### load custom R script to load patient metadata, which covers in detail
# loading RDS files
# load patient metadata and convert into numerical data
# add column "Survival Time"
# convert data frames as matrices
# convert NaNs and NAs
# create df.total and save

### load PRI features rds file as df. for training and validation set
printf("Loading RDS from data set: %s", dataset[dataset.id])
df.training <- readRDS(training.data.path)
df.validation <- readRDS(validation.data.path)
printf("Dimensions of training and validation set: ")
print(dim(df.training))
print(dim(df.validation))

#safe rownames ("Patient ID")
rownames.training <- row.names(df.training)
rownames.validation <- row.names(df.validation)

typeColNum <- 1

### ---------------------- loading patient data from excel
patient_data <- read_excel(patient.data.path)
#View(patient_data)

###### PREPARE patient data #####################
#set columstype from patien data str to numeric values
#colums with numeric value : Age at Diagnosis, WBC Count, Date of Diagnosis, Time to Relapse(Day), CCR (Day) 
patient_data$`Time to Relapse (Days)` <- as.numeric(patient_data$`Time to Relapse (Days)`)
patient_data$`CCR (Days)` <- as.numeric(patient_data$`CCR (Days)`)
patient_data$`Time to Relapse (Days)`[is.na(patient_data$`Time to Relapse (Days)`)] <- 0
patient_data$`CCR (Days)`[is.na(patient_data$`CCR (Days)`)] <- 0

#new column "SurvivalTime""
patient_data$SurvivalTime <- patient_data$`Time to Relapse (Days)` + patient_data$`CCR (Days)`

#yes2 == yes
patient_data$`Relapse Status`[patient_data$`Relapse Status` == "Yes2"] <- "Yes"

#Set DDPRStatus to nuermic binary values
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "Low"] <- 0
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "High"] <- 1

#reduce patien data to necessary colums
cohort = patient_data[, c(1, 15, 17, 11, 16, 7, 8, 9)]
#rename column names
colnames(cohort) <- c("PatientID", "Cohort", "SurvivalTime", "RelapseStatus", "DDPRStatus", "FinalRisk", "MRDRisk", "NCIRisk")

#divide cohort into training and validation set
training.set <- cohort[which(cohort$Cohort == "Training"), 1:6]
validation.set <- cohort[which(cohort$Cohort == "Validation"), 1:6]

#bind training and validation set to totalset
total.set <- bind_rows(training.set, validation.set)
#################################################


### -------------- preparation in training set 
#relapse status as numeric
cond.training <- training.set$RelapseStatus[which(training.set$PatientID %in% row.names(df.training))]
cond.training[cond.training %in% c("Yes")] <- 1
cond.training[cond.training %in% c("No")] <- 0
cond.training <- as.numeric(cond.training)
#survival time in training.set as numeric
survival.training <- training.set$SurvivalTime[which(training.set$PatientID %in% row.names(df.training))]
survival.training <- as.numeric(survival.training)
#DDPR STATUS as numeric
DDPR.training <- training.set$DDPRStatus[which(training.set$PatientID %in% row.names(df.training))]
DDPR.training <- as.numeric(DDPR.training)

### -------------- preparation in validation set 
cond.validation <- validation.set$RelapseStatus[which(validation.set$PatientID %in% row.names(df.validation))]
cond.validation[cond.validation %in% c("Yes")] <- 1
cond.validation[cond.validation %in% c("No")] <- 0
cond.validation <- as.numeric(cond.validation)
survival.validation <- validation.set$SurvivalTime[which(validation.set$PatientID %in% row.names(df.validation))]
survival.validation <- as.numeric(survival.validation)
DDPR.validation <- validation.set$DDPRStatus[which(validation.set$PatientID %in% row.names(df.validation))]
DDPR.validation <- as.numeric(DDPR.validation)

### introduce randomized cond.trainings
if (randomize.label) {
  # set.seed(seed.vec[4])
  cond.training <- sample(cond.training)
  cond.validation <- sample(cond.validation)
}
### introduce spike-ins
if (spikeIns > 0) {
  spike1 <- cond.training * 1
  df.training <- bind_cols(df.training, as.data.frame(spike1)
  )
  spike1 <- cond.validation * 1
  df.validation <- bind_cols(df.validation, as.data.frame(spike1)
  )
}
if (spikeIns > 1) {
  spike2 <- cond.training * 4
  df.training <- bind_cols(df.training, as.data.frame(spike2)
  )
  spike2 <- cond.validation * 4
  df.validation <- bind_cols(df.validation, as.data.frame(spike2)
  )
}

### add cond.training/survivaltime/age to df.training
df.training <- bind_cols(as.data.frame(DDPR.training), df.training)
colnames(df.training)[typeColNum] <- "DDPRStatus"
df.training <- bind_cols(as.data.frame(cond.training), df.training)
colnames(df.training)[typeColNum] <- "RelapseStatus"
df.training <- bind_cols(as.data.frame(survival.training), df.training)
colnames(df.training)[typeColNum] <- "Survivaltime"
#do same for validation set
df.validation <- bind_cols(as.data.frame(DDPR.validation), df.validation)
colnames(df.validation)[typeColNum] <- "DDPRStatus"
df.validation <- bind_cols(as.data.frame(cond.validation), df.validation)
colnames(df.validation)[typeColNum] <- "RelapseStatus"
df.validation <- bind_cols(as.data.frame(survival.validation), df.validation)
colnames(df.validation)[typeColNum] <- "Survivaltime"


### glmnet need input matrix for model 
# safe df as matrix
df.training=as.matrix(df.training)
df.validation = as.matrix(df.validation)

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

#### Check if df.total already exists as RDS (binding is bottleneck)
#### RDS with df.total will be created once if it doesnt exist and then read every run
if (file.exists(total.data.path)) {
  df.total <- readRDS(total.data.path)
  printf("df.total already exists...using created RDS.")
} else if (!randomize.label & spikeIns == 0) {
  printf("df.total does NOT exist...creating RDS. Takes a while..")
  require("data.table")
  # names(df.training) <- colnames(df.training)
  # names(df.validation) <- colnames(df.validation)
  df1 <- as.data.table(df.training)
  df2 <- as.data.table(df.validation)
  df.total <- rbind(df1, df2)
  df.total <- as.data.frame(df.total)
  rownames(df.total) <- c(rownames(df.training), rownames(df.validation))
  saveRDS(df.total, file=total.data.path)
  rm(df1, df2)
}
rownames(df.total) <- c(rownames.training, rownames.validation)

printf("Number of patients in total: %s", nrow(df.total))
#View(df.total)
#total 60 patients in cohort
#Basal no condition total 54 patients
#Basal condition 1.1 total 41 patients
#BCR 45 patients total, divided in 37 in training and 8 in validation


printf("Randomizing patient outcomes: %s.", randomize.label)
printf("Introducing two spikeIn columns: %s.", spikeIns)
print("Done preprocessing.")