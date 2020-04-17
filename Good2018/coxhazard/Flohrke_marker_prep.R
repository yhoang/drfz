#!/usr/bin/R
# Authors: Felix Lohrke and Yen Hoang
# function that checks the differences between two marker lists and then removes the quadrants which incorporates these markers
# difference from dataframe read ins (columns)

rm(list = ls())
options(max.print = 100)

### ---------- input
# working.station   FL, YH; different path settings
# project.id        1:5; for Basal, BCR, IL7, Pervanadate, TSLP
# sub.sample.name   func, func_plus3, func_plus6; different marker combinations

working.station <- "FL"
working.station <- "YH"

sub.sample.name <- "func_plus3"
project.id  <- 2


if (working.station == "FL") {
    Rdata.path <- file.path("", "home", "felix", "AG_Baumgrass", "Scripts", "Pri", "Pri_good_established", "Rds")
    Text.path <- file.path("", "home", "felix", "AG_Baumgrass", "Data", "Good", "Marker_combinations")
} else if (working.station == "YH") {
    Rdata.path <- file.path("D:", "drfz", "Good2018", "Rdata")
    Text.path <- file.path("D:", "drfz", "Good2018", "tables")
}


#a useful print function
printf <- function(...) invisible(print(sprintf(...)))
project.name <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")


# define output path and filename
outfile.training <- paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Training_", sub.sample.name, "_quadrant_absRange_cof0.2.rds")
outfile.validation <- paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Validation_", sub.sample.name, "_quadrant_absRange_cof0.2.rds")

# list names of all markers 
markers.full <- scan(paste0(Text.path, "/columns_full.txt"), what="character", sep="\n")

# list names of markers to extract (by removing difference from full)
markers.func <- scan(paste0(Text.path, "/columns_", sub.sample.name, ".txt"), what="character", sep="\n")

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

# Training and validation sets from which markers not in sub sample are to be removed
df.training.full <- readRDS(paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Training_full_quadrant_absRange_cof0.2.rds"))
df.validation.full <- readRDS(paste0(Rdata.path, "/", project.name[project.id], "/", project.name[project.id], "_Validation_full_quadrant_absRange_cof0.2.rds"))

# saving rownames
rownames.training <- row.names(df.training.full)
rownames.validation <- row.names(df.validation.full)

# saving marker names of dataframes
names.df.training <- names(df.training.full)
names.df.validation <- names(df.validation.full)

# saving as new df from read ins to change
df.training.func <- as.data.frame(df.training.full)
df.validation.func <- as.data.frame(df.validation.full)

### training
print("#### Working on training set ####")
for (i in 1:length(markers.notfunc)) {
    # index pos of markers to be removed
    printf("Marker that is selected to be removed: %s", markers.notfunc[i])
    ind.train <- grep(markers.notfunc[i], names.df.training)
    
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
    printf("Marker that is selected to remove: %s", markers.notfunc[i])
    ind.val <- grep(markers.notfunc[i], names.df.validation)
    
    # only execute if to removed markers are still in the dataframe
    if (length(ind.val) != 0){
        df.validation.func <- df.validation.func[, -ind.val]
        names.df.validation <- names(df.validation.func)
        # printf("PRI features left: %s", ncol(df.validation.func))
    } else {
        print("marker not found")
    }

}
printf("PRI features left: %s", ncol(df.validation.func))

# showing new dimensions of dataframes
print("Removal succesful.")
printf("Old dimension of Training set: %s", paste(dim(df.training.full), collapse = "  "))
printf("New dimension of Training set: %s", paste(dim(df.training.func), collapse = "  "))
printf("Old dimension of Validation set: %s", paste(dim(df.validation.full), collapse = "  "))
printf("New dimension of Validation set: %s", paste(dim(df.validation.func), collapse = "  "))

# adding rownames
row.names(df.training.func) <- rownames.training
row.names(df.validation.func) <- rownames.validation

# saving modified training and validation dataframes
saveRDS(df.training.func, outfile.training)
saveRDS(df.validation.func, outfile.validation)
printf("Saved new df as RDS in %s", outfile.training)
printf("Saved new df as RDS in %s", outfile.validation)

