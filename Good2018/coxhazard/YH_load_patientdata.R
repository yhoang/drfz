#!/usr/bin/R


#loading patient data from excel
patient_data <- read_excel(patient.data.path)
#View(patient_data)

###### PREPARE patient data #####################
#set columstype from patien data str to numeric values
#colums with numeric value : Age at Diagnosis, WBC Count, Date of Diagnosis, Time to Relapse(Day), CCR (Day) 
patient_data$`Time to Relapse (Days)` <- as.numeric(patient_data$`Time to Relapse (Days)`)
patient_data$`CCR (Days)` <- as.numeric(patient_data$`CCR (Days)`)
patient_data$`Time to Relapse (Days)`[is.na(patient_data$`Time to Relapse (Days)`)] <- 0
patient_data$`CCR (Days)`[is.na(patient_data$`CCR (Days)`)] <- 0

#new column "Survival Time (Day)""
patient_data$`Survival Time (Day)` <- patient_data$`Time to Relapse (Days)` + patient_data$`CCR (Days)`

#yes2 == yes
patient_data$`Relapse Status`[patient_data$`Relapse Status` == "Yes2"] <- "Yes"

#Set DDPR Status to nuermic binary values
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "Low"] <- 0
patient_data$`DDPR Risk`[patient_data$`DDPR Risk` == "High"] <- 1

#reduce patien data to necessary colums
cohort=patient_data[, c(1, 11, 8, 16, 15, 17)]

#divide cohort into training and validation set
training.set=cohort[which(cohort$Cohort == "Training"), 1:6]
validation.set=cohort[which(cohort$Cohort == "Validation"), 1:6]

#bind training and validation set to totalset
total.set <- bind_rows(training.set, validation.set)
#################################################
