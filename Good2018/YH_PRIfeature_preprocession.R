#!/usr/bin/R

### ------------------- load PRI features and patient meta data and preprocess
### load custom R script to load patient metadata, which covers in detail
# loading RDS files
# load patient metadata and convert into numerical data
# add column "Survival Time"
# convert data frames as matrices
# convert NaNs and NAs
# create df.total and save

###### PREPARE patient data #####################
### ---------------------- loading patient data from excel
patient_data <- read_excel(patient.data.path)
#View(patient_data)

#set columstype from patien data str to numeric values
#colums with numeric value : Age at Diagnosis, WBC Count, Date of Diagnosis, Time to Relapse(Day), CCR (Day) 
patient_data$`Time to Relapse (Days)` <- as.numeric(patient_data$`Time to Relapse (Days)`)
patient_data$`CCR (Days)` <- as.numeric(patient_data$`CCR (Days)`)
patient_data$`Time to Relapse (Days)`[is.na(patient_data$`Time to Relapse (Days)`)] <- 0
patient_data$`CCR (Days)`[is.na(patient_data$`CCR (Days)`)] <- 0

#new column "Survivaltime""
patient_data$Survivaltime <- patient_data$`Time to Relapse (Days)` + patient_data$`CCR (Days)`

#yes2 == yes
patient_data$`Relapse Status`[patient_data$`Relapse Status` == "Yes2"] <- "Yes"

#Set DDPRStatus to nuermic binary values
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "Low"] <- 0
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "High"] <- 1

#reduce patient data to necessary colums
cohort <- patient_data[, c(1, 15, 17, 11, 16, 7, 8, 9, 3, 4)]
#rename column names
colnames(cohort) <- c("PatientID", "Cohort", "Survivaltime", "RelapseStatus", "DDPRStatus", "FinalRisk", "MRDRisk", "NCIRisk", "Age", "WBC")
cohort <- cohort[-which(cohort$Cohort == "NA"), ]
cohort$DDPRStatus <- as.numeric(cohort$DDPRStatus)

cohort$RelapseStatus[cohort$RelapseStatus %in% "No"] <- 0
cohort$RelapseStatus[cohort$RelapseStatus %in% "Yes"] <- 1

typeColNum <- 1
if (dataset.id == 6) {
  printf("Loading table from data set: %s", dataset[dataset.id])
  # skip first two lines
  # use third line as header
  # read only 78 lines, after that are comments and description
  df.total <- read.table(file.path(Rdata.path, "Good_suppl2.csv"), skip = 2, sep = ",", header = TRUE, nrows = 78)
  colnames(df.total)[1] <- "PatientID"
  df.total <- as_tibble(df.total)

  ### Features in DDPR
  # For each of 5 expanded populations (Pre-Pro-B, Pro-BI, Pro-BII, Pre-BI, Pre-BII, Early Progenitors) --> Pop5-Pop9, Mixed2
  # 1 Percent of cells in population
  # 2 Mean arsinh-transformed expression of 24 phenotypic markers
  # 3 For each of 9 phosphorylated proteins,
  # % positive cells in Basal condition and in response to 4 perturbations (Pervanadate, IL-7, BCR Crosslink, TSLP)
  # Plus 2 clinical features: age and white blood cell count at diagnosis (see Supplementary Table 1).
  # Total: 5 x (1 + 24 + 9 x 5) + 2 = 352

  ### FIRST change label of pPLCg1_2 to pPLCg1.2
  colnames(df.total) <- sub("^pPLCg1_2", "pPLCg1.2", colnames(df.total))

  ### SECOND get column indices
  df.split <- strsplit(colnames(df.total), "_")
  # extract with sapply the 2nd column from df.split
  ### 1: get percent of cells in population
  # n = 1
  pop.idx <- which(sapply(df.split, "[", 2) == "P")
  ### 3: get phosphorylated proteins (9) % positive cells in Basal condition and in response to 4 perturbations 
  # n = 9 * 5 = 45
  dp.idx <- grep("dP", sapply(df.split, "[", 2))

  ddpr.idx <- vector()
  for (i in 5:9) {
    # get the values for the population i
    grep.label <- paste0("Pop", i)
    grep.idx <- grep(grep.label, colnames(df.total))

    ### 2: get MSI expressions of 24 phenotypic markers
    # n = 24
    msi.idx <- grep(grep.label, sapply(df.split, "[", 2))

    tmp.idx <- intersect(grep.idx, c(pop.idx, msi.idx, dp.idx))
    # combine into DDPR.idx
    ddpr.idx <- c(ddpr.idx, tmp.idx)
  }
  df.total <- df.total[, c(1, ddpr.idx)]
  
  df.total$PatientID <- as.character(df.total$PatientID)


  ### reduce to 54 samples
  if (TRUE) {
    pat.idx <- match(cohort$PatientID, df.total$PatientID)
    # df.total <- df.total[pat.idx, ]

    # df.total$Age <- cohort$Age
    # df.total$WBC <- as.numeric(cohort$WBC)
    # df.total$DDPRStatus <- cohort$DDPRStatus
    # df.total$RelapseStatus <- cohort$RelapseStatus
    # df.total$RelapseStatus[df.total$RelapseStatus %in% "No"] <- 0
    # df.total$RelapseStatus[df.total$RelapseStatus %in% "Yes"] <- 1
    # df.total$RelapseStatus <- as.numeric(df.total$RelapseStatus)
    # df.total$Survivaltime <- cohort$SurvivalTime

    df.total <- df.total %>%
      slice(pat.idx) %>%
      add_column(Age = cohort$Age,
                 WBC = as.numeric(cohort$WBC), 
                 Survivaltime = cohort$Survivaltime,
                 RelapseStatus = as.numeric(cohort$RelapseStatus),
                 DDPRStatus = cohort$DDPRStatus)


  } else {
    
    ### THIRD get row indices
    df.split <- strsplit(df.total$PatientID, "-")
    pat.idx <- match(sapply(df.split, "[", 1), cohort$PatientID)
    ### 4: add age and white blood cell count at diagnosis from patient_table 
    df.total$Age <- df.total$WBC <- df.total$DDPRStatus <- df.total$RelapseStatus <- df.total$Survivaltime <- NA
    df.total$Age <- cohort$Age[pat.idx]
    df.total$WBC <- as.numeric(cohort$WBC)[pat.idx]

    ### add condition/survivaltime/age to df.total
    df.total$DDPRStatus <- cohort$DDPRStatus[pat.idx]
    df.total$RelapseStatus <- cohort$RelapseStatus[pat.idx]
    df.total$Survivaltime <- cohort$Survivaltime[pat.idx]

    #### check for NAs
    if (FALSE) {
      isNA.row <- list()
      for (row in 1:nrow(df.total)) {
        isNA.row[[row]] <- which(is.na(df.total[row,]))
      }
      ### reminder
      # 38 out of 78 rows have NAs
      # 20 out of 78 rows have >90 NAs
      # lenNA <- sapply(isNA.row, function(x) length(x))
      # NA.idx <- which(lenNA > 10)
    }
    print("Removing NA positions in condition RelapseStatus..")
    df.total$RelapseStatus[which(df.total$RelapseStatus == "NA")] <- NA
    df.total <- df.total[-which(is.na(df.total$RelapseStatus)), ]
    # df.total$RelapseStatus <- as.factor(df.total$RelapseStatus)

    print("Setting NAs in WBC (UPN90) to patient UPN22 with similar RelapseStatus, Age and Survivaltime..")
    WBC.NA.idx <- which(rownames(df.total) == "UPN90")
    WBC.time.idx <- which(df.total$Survivaltime < df.total$Survivaltime[WBC.NA.idx] + 300 & df.total$Survivaltime > df.total$Survivaltime[WBC.NA.idx] - 300)
    WBC.age.idx <- which(df.total$Age < df.total$Age[WBC.NA.idx] + 3 & df.total$Age > df.total$Age[WBC.NA.idx] - 3)
    WBC.status.idx <- which(df.total$RelapseStatus == df.total$RelapseStatus[WBC.NA.idx])
    # intersect
    WBC.idx <- Reduce(intersect, list(WBC.age.idx, WBC.time.idx, WBC.status.idx))
    WBC.NA.idx <- grep("UPN90", rownames(df.total))
    df.total$WBC[WBC.NA.idx] <- mean(df.total$WBC[WBC.idx[which(WBC.idx != WBC.NA.idx)]])
    
    # change Yes and No to 1/0
    df.total$RelapseStatus[df.total$RelapseStatus %in% "Yes"] <- 1
    df.total$RelapseStatus[df.total$RelapseStatus %in% "No"] <- 0
    df.total$RelapseStatus <- as.numeric(df.total$RelapseStatus)
  }

  print("Convert column NAs to sample group mean")
  status.idx <- which(colnames(df.total) == "RelapseStatus")
  # ptm <- proc.time()
  ### convert NAs to sample group mean
  for (i in 1:ncol(df.total)) {
    NA.idx <- which(is.na(df.total[, i]))
    for (j in NA.idx) {
      # select group
      # tmp <- df.total[which(df.total[[j, status.idx]] == df.total[, status.idx]), i]
      # remove NA
      # tmp <- tmp %>% na.omit()
      # select group and remove NA
      tmp <- df.total %>% filter(RelapseStatus == df.total[[j, status.idx]]) %>% select(all_of(i)) %>% drop_na()
      # sample a number between mean +- sd
      tmp <- unlist(tmp)

      df.total[j, i] <- round(sample(seq(mean(tmp) - sd(tmp), mean(tmp) + sd(tmp), by = 0.01), 1), 2)
    }
  }
  # print(proc.time())

  # }

} else {
  ### load PRI features rds file as df. for training and validation set
  printf("Loading RDS from data set: %s", dataset[dataset.id])
  df.training <- readRDS(training.data.path)
  df.validation <- readRDS(validation.data.path)
  printf("Dimensions of training and validation set: ")
  print(dim(df.training))
  print(dim(df.validation))

  pat.idx.training <- match(rownames(df.training), cohort$PatientID)
  pat.idx.validation <- match(rownames(df.validation), cohort$PatientID)

  #divide cohort into training and validation set
  training.set <- cohort[which(cohort$Cohort == "Training"), 1:6]
  validation.set <- cohort[which(cohort$Cohort == "Validation"), 1:6]

  #################################################

  ### -------------- preparation in training set 
  #relapse status as numeric
  cond.training <- training.set$RelapseStatus[which(training.set$PatientID %in% row.names(df.training))]
  cond.training[cond.training %in% c("Yes")] <- 1
  cond.training[cond.training %in% c("No")] <- 0
  cond.training <- as.numeric(cond.training)
  #survival time in training.set as numeric
  survival.training <- training.set$Survivaltime[which(training.set$PatientID %in% row.names(df.training))]
  survival.training <- as.numeric(survival.training)
  #DDPR STATUS as numeric
  DDPR.training <- training.set$DDPRStatus[which(training.set$PatientID %in% row.names(df.training))]
  DDPR.training <- as.numeric(DDPR.training)

  ### -------------- preparation in validation set 
  cond.validation <- validation.set$RelapseStatus[which(validation.set$PatientID %in% row.names(df.validation))]
  cond.validation[cond.validation %in% c("Yes")] <- 1
  cond.validation[cond.validation %in% c("No")] <- 0
  cond.validation <- as.numeric(cond.validation)
  survival.validation <- validation.set$Survivaltime[which(validation.set$PatientID %in% row.names(df.validation))]
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
  df.training <- df.training %>%
    add_column(PatientID = cohort$PatientID[pat.idx.training],
               Survivaltime = cohort$Survivaltime[pat.idx.training],
               RelapseStatus = as.numeric(cohort$RelapseStatus[pat.idx.training]),
               DDPRStatus = cohort$DDPRStatus[pat.idx.training])

  ### add cond.training/survivaltime/age to df.validation
  df.validation <- df.validation %>%
    add_column(PatientID = cohort$PatientID[pat.idx.validation],
               Survivaltime = cohort$Survivaltime[pat.idx.validation],
               RelapseStatus = as.numeric(cohort$RelapseStatus[pat.idx.validation]),
               DDPRStatus = cohort$DDPRStatus[pat.idx.validation])

  ### glmnet need input matrix for model 
  # safe df as matrix
  df.training <- as.matrix(df.training)
  df.validation <- as.matrix(df.validation)

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
    df.total <- rbindlist(list(df1, df2))
    df.total <- as.data.frame(df.total)
    # rownames(df.total) <- c(rownames(df.training), rownames(df.validation))
    saveRDS(df.total, file=total.data.path)
    rm(df1, df2)
  }

  printf("Dimension df.total should be 24292: %s", dim(df.total))
  #View(df.total)
  #total 60 patients in cohort
  #Basal no condition total 54 patients
  #Basal condition 1.1 total 41 patients
  #BCR 45 patients total, divided in 37 in training and 8 in validation
}


status.idx <- which(colnames(df.total) == "RelapseStatus")
ddpr.idx <- which(colnames(df.total) == "DDPRStatus")
time.idx <- which(colnames(df.total) == "Survivaltime")
sample.idx <- which(colnames(df.total) == "PatientID")

printf("Randomizing patient outcomes: %s.", randomize.label)
printf("Introducing two spikeIn columns: %s.", spikeIns)


df.tmp <- df.total[, -c(status.idx, ddpr.idx, time.idx, sample.idx)]
df.mean <- mean(unlist(df.tmp))
df.sd <- sd(unlist(df.tmp))
### collect Z score
df.zscore <- (df.tmp - df.mean) / df.sd
# hist(unlist(df.zscore), xlim=c(-0.06,0.05), n=10000)
# hist(unlist(df.zscore), xlim=c(-4,4), n=1000)
  
print("Done preprocessing.")


df.relapsed <- df.total %>%
  select(RelapseStatus, Survivaltime) %>%
  filter(RelapseStatus == 1) %>%
  arrange(Survivaltime)

df.unrelapsed <- df.total %>%
  select(RelapseStatus, Survivaltime) %>%
  filter(RelapseStatus == 0) %>%
  arrange(Survivaltime)


if (FALSE) {
  # Change histogram plot line colors by groups
  ggplot(df.total[, 354:355], aes(x=Survivaltime, color=as.factor(RelapseStatus))) +
    geom_histogram(fill="white", binwidth=500, alpha=0.5, position="dodge")


}