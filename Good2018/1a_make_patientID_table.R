#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019 - 2020
# Goood2018

rm(list = ls())
options(max.print = 100)

### ----------- initiate
subgroups <- c("Training", "Validation")
cofactor <- 0.2
coverage <- "full"
stat.id <- 1
### set paths
folder.path <- file.path("D:", "drfz", "Good2018")
db.path <- file.path("D:", "DB")
################


### --------- load libraries
# RSQLite     :  interact with database
load.libraries <- c("RSQLite", "dplyr")
#install.packages(libraries, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(load.libraries, require, character.only = TRUE)
### load functions
source(file.path(folder.path, "PRI_funs.R"))

setwd(folder.path)
project.name <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")
project.name.long <- c("Basal", "BCR-Crosslink", "IL-7", "Pervanadate", "TSLP")
stat.info <- c("absRange", "variance", "freq_green", "mean")

### metatable file
cohort_full <- read.csv(sprintf("%s/tables/patient_cohort.csv", folder.path))
cohort <- cohort_full[, c(1, 7:9, 16, 15, 11, 12, 14)]
# "Patient.ID" 
# "Final.Risk"    
# "MRD.Risk"
# "NCI.Rome.Risk"
# "DDPR.Risk"
# "Cohort" 
# "Relapse.Status"
# "Time to Relapse (Days)"
# "CCR (Days)"
sub.set <- cohort[which(cohort$Cohort == subgroups[1] | cohort$Cohort == subgroups[2]), ]
sub.set$Relapse.Status[which(sub.set$Relapse.Status == "Yes2")] <- "Yes"
sub.set$Relapse.Status <- factor(sub.set$Relapse.Status)
sub.set$Patient.ID <- factor(sub.set$Patient.ID)
# combine column 6+7 (Time to Relapse and CCR)
sub.set$Time.to.Relapse..Days. <- as.numeric(sub.set$Time.to.Relapse..Days.)
sub.set$Time.to.Relapse..Days.[is.na(sub.set$Time.to.Relapse..Days.)] <- 0
sub.set$CCR..Days. <- as.numeric(sub.set$CCR..Days.)
sub.set$CCR..Days.[is.na(sub.set$CCR..Days.)] <- 0
sub.set$Survival.Time <- sub.set$Time.to.Relapse..Days. + sub.set$CCR..Days.
sub.set <- sub.set[, - c(8, 9)]
# sort to Cohort Training and Validation
sub.set <- sub.set[order(sub.set$Cohort), ]


### ------ load files from database
db.name <- "RB_20191002_Good2018.sqlite3"
fcs$connectDb(file.path(db.path, db.name))

for (project.id in 1:length(project.name)) {
  ### get files from DB
  files.in.DB <- fcs$getDFtable(paste0(project.name[project.id], "_fileIdentity"))[2]

  ### add column for TRUE/FALSE if file is found in each project in DB
  sub.set[, paste0(project.name[project.id], "_DB")] <- sapply(1:nrow(sub.set), function(x) {
    search.file.in.DB <- paste0(sub.set[x, 1], "_", project.name.long[project.id], ".fcs")

    search.file.in.DB2 <- paste0(sub.set[x, 1], "_", tolower(project.name.long[project.id]), ".fcs")

    if (length(which(files.in.DB == search.file.in.DB)) > 0 | length(which(files.in.DB == search.file.in.DB2)) > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  
  ### load quad tables in temp.data.all
  # training set
  quad.table.train <- sprintf("%s/Rdata/%s/%s_%s_%s_quadrant_%s_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroups[1], coverage, stat.info[stat.id], cofactor)
  temp.data.train <- readRDS(quad.table.train)[, 1:3]
  # validation set
  quad.table.val <- sprintf("%s/Rdata/%s/%s_%s_%s_quadrant_%s_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroups[2], coverage, stat.info[stat.id], cofactor)
  temp.data.val <- readRDS(quad.table.val)[, 1:3]
  ### combine both sets
  temp.data.all <- rbind(temp.data.train, temp.data.val)

  ### add column for TRUE/FALSE if file is found in each project in PRI
  sub.set[, paste0(project.name[project.id], "_PRIquad")] <- sapply(1:nrow(sub.set), function(x) {
    if (length(which(rownames(temp.data.all) == as.character(sub.set[x, 1]))) > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
}
dbDisconnect(fcs$conn)

### save overview table
write.csv(sub.set, sprintf("%s/tables/patient_cohort_overview.csv", folder.path, row.names=FALSE))
