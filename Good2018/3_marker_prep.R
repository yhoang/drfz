#!/usr/bin/R
# Authors: Felix Lohrke and Yen Hoang
# function that checks the differences between two marker lists and then removes the quadrants which incorporates these markers
# difference from dataframe read ins (columns)

rm(list = ls())
options(max.print = 100)

### ---------- input
# working.station   FL, YH; different path settings
# dataset.id        1:5; for Basal, BCR, IL7, Pervanadate, TSLP
# subset.id         1:4 for full, func, func3, func6; 
# subject.in.id     1:13 for "df", "quadrant", "ROC", "AUC", "iAUC",  "CompTotal", "Variables", "Heat", "CompAll", "KM", "Summary"
# feature.id        1:4 for "absRange", "variance", "freq.green", "mean"

# working.station <- "FL"
working.station <- "YH"
initdate <- "YH200526"
dataset.id <- 3
subset.id <- 4
subject.in.id <- 2
feature.id <- 1
comment.in <- "autoSec.cof0.2"
comment.train <- comment.val <- comment.in
comment.out <- comment.in

### load custom R script to generate standardized file names
source(file.path("D:", "drfz", "Good2018", "YH_func_naming.R"))

if (working.station == "FL") {
    Rdata.path <- file.path("", "home", "felix", "AG_Baumgrass", "Scripts", "Pri", "Pri_good_established", "Rds")
    Text.path <- file.path("", "home", "felix", "AG_Baumgrass", "Data", "Good", "Marker_combinations")
} else if (working.station == "YH") {
    Rdata.path <- file.path("D:", "drfz", "Good2018", "Rdata", dataset[dataset.id])
    Text.path <- file.path("D:", "drfz", "Good2018", "tables")
}

#a useful print function
printf <- function(...) invisible(print(sprintf(...)))



######## INPUT FILES FULL #####################
# subgroup <- c("Training", "Validation", "Total")
training.in.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(1, com = comment.train, subs = subset[1]))
validation.in.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(2, com = comment.val, subs = subset[1]))
##########################

######## OUTPUT FILES func/func3/func6 #####################
# subgroup <- c("Training", "Validation", "Total")
training.out.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(1, com = comment.train))
validation.out.path <- sprintf("%s/%s.rds", Rdata.path, input_filenaming(2, com = comment.val))
##########################

### load PRI features rds file as df. for training and validation set
printf("Loading full RDS from data set: %s", dataset[dataset.id])
df.training.full <- readRDS(training.in.path)
df.validation.full <- readRDS(validation.in.path)
printf("Dimensions of training and validation set: ")
print(dim(df.training.full))
print(dim(df.validation.full))


# list names of all markers 
markers.full <- scan(paste0(Text.path, "/columns_full.txt"), what="character", sep="\n")

# list names of markers to extract (by removing difference from full)
markers.func <- scan(paste0(Text.path, "/columns_", subset[subset.id], ".txt"), what="character", sep="\n")

# list names of markers that are in full but not func (used to remove them from dataframe)
markers.notfunc <- markers.full[!(markers.full %in% markers.func)]

# comparing marker inputs
printf("markers full: %s", length(markers.full))
print(markers.full)
print("######################################")

printf("markers func: %s", length(markers.func))
print(markers.func)
print("######################################")

printf("markers not in func: %s", length(markers.notfunc))
print(markers.notfunc)
print("######################################")

# saving as new df from read ins to change
df.training.func <- as.data.frame(df.training.full)
df.validation.func <- as.data.frame(df.validation.full)
# saving rownames
rownames.training <- row.names(df.training.full)
rownames.validation <- row.names(df.validation.full)
# saving marker names of dataframes
names.df.training <- names(df.training.full)
names.df.validation <- names(df.validation.full)

if (TRUE) {
    print("Check if all markers are present in full and subset..")
    ### test if all markers are present in full and subsample
    strings <- unlist(strsplit(names.df.training, "[.]"))
    strings.tab <- table(strings)
    strings.tab <- strings.tab[-which(names(strings.tab) %in% c("absRange", "Q1", "Q2", "Q3", "Q4"))]
    markers.in.rds <- names(strings.tab)
    # should be true
    all(markers.in.rds == markers.full)
    all(markers.in.rds[which(markers.in.rds %in% markers.notfunc)] == markers.notfunc)
}


### training
print("#### Working on training set ####")
for (i in 1:length(markers.notfunc)) {
    # index pos of markers to be removed
    printf("%st/%s marker selected to be removed: %s", i, length(markers.notfunc), markers.notfunc[i])
    ind.train <- grep(paste0(markers.notfunc[i], "\\."), names.df.training)
    
    # only execute if to removed markers are still in the dataframe
    if (length(ind.train) != 0){
        df.training.func <- df.training.func[, -ind.train]
        names.df.training <- names(df.training.func)
        # printf("PRI features left: %s", ncol(df.training.func))
    }  else {
        printf("Marker not found: %s", markers.notfunc[i])
    }
}
printf("PRI features left: %s", ncol(df.training.func))


### validation
print("#### Working on validation set ####")
for (i in 1:length(markers.notfunc)) {
    # index pos of markers to be removed
    printf("%st/%s marker selected to be removed: %s", i, length(markers.notfunc), markers.notfunc[i])
    ind.val <- grep(paste0(markers.notfunc[i], "\\."), names.df.validation)
    
    # only execute if to removed markers are still in the dataframe
    if (length(ind.val) != 0){
        df.validation.func <- df.validation.func[, -ind.val]
        names.df.validation <- names(df.validation.func)
        # printf("PRI features left: %s", ncol(df.validation.func))
    } else {
        printf("Marker not found: %s", markers.notfunc[i])
    }

}
printf("PRI features left: %s", ncol(df.validation.func))

# showing new dimensions of dataframes
print("Removal succesfull.")
printf("Old dimension of Training set: %s", paste(dim(df.training.full), collapse = "  "))
printf("New dimension of Training set: %s", paste(dim(df.training.func), collapse = "  "))
printf("Old dimension of Validation set: %s", paste(dim(df.validation.full), collapse = "  "))
printf("New dimension of Validation set: %s", paste(dim(df.validation.func), collapse = "  "))

len.marker <- length(markers.func)
# m * (m - 1) * (m - 2) * 0.5 * 4
len.features <- len.marker * (len.marker - 1) * (len.marker - 2) * 0.5 * 4
printf("This data should have %s markers resulting to %s features.", len.marker, len.features)

# adding rownames
row.names(df.training.func) <- rownames.training
row.names(df.validation.func) <- rownames.validation

# saving modified training and validation dataframes
if (file.exists(training.out.path) | file.exists(validation.out.path)) {
    print("FILE EXISTS ALREADY. CHECK YOUR PATH!")
    print(training.out.path)
    print(validation.out.path)
} else {
    saveRDS(df.training.func, training.out.path)
    saveRDS(df.validation.func, validation.out.path)
    printf("Saved new df as RDS in %s", training.out.path)
    printf("Saved new df as RDS in %s", validation.out.path)
}

