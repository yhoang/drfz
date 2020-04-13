#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019 - 2020
# Goood2018

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high risk

rm(list = ls())
options(max.print = 100)

### ----------- initiate
# cofactor      1, 0.2, 0.1 [so far at 0.2]
# mincells      2,5,10,20 [so far at 5]
# stat.id       1:4 [absRange, variance, freq_green, mean]
# cluster.size  1:12, detect with detectCores()
# project.id    1:5 for Basal, BCR, IL7, Pervanadate, TSLP
# subgroup      "Validation", "Training"
# coverage      full, func, func-plus3, func-plus6 [func=18, func-plus3=21, func-plus6=24]
# check.cutoffs TRUE/FALSE, check if cutoffs are set in all files
cofactor <- 0.2
mincells <- 5
stat.id <- 1
cluster.size <- 3
project.id <- 4
subgroup <- "Training"
# subgroup <- "Validation"
coverage <- "full"
check.cutoffs <- TRUE
### set paths
folder.path <- file.path("D:", "drfz", "Good2018")
db.path <- file.path("D:", "DB")
################




### --------- load libraries
# xlsx      :  Library for Excel reading and creating
# RSQLite     :  interact with database
# reshape2    :  function melt()
# dplyr      :  faster binding of columns bind_cols() and rows bind_rows()
# foreach     :  allows for each loops
# doParallel   :  use several clusters for parallele calculations
load.libraries <- c("RSQLite", "dplyr", "reshape2", "foreach", "doParallel")
#install.packages(libraries, lib = "C:/Program Files/R/R-3.6.1/library")
lapply(load.libraries, require, character.only = TRUE)
### load functions
source(file.path(folder.path, "PRI_funs.R"))
printf <- function(...) invisible(print(sprintf(...)))


setwd(folder.path)
project.name <- c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")
project.name.long <- c("Basal", "BCR-Crosslink", "IL-7", "Pervanadate", "TSLP")
stat.info <- c("absRange", "variance", "freq_green", "mean")

### metatable file
# cohort_full <- read.csv(sprintf("%s/tables/patient_cohort.csv", folder.path))
# sub.set <- cohort_full[, c(1, 5, 8, 11, 16, 15)]
sub.set <- read.csv(sprintf("%s/tables/%s_%s_patient_cohort.csv", folder.path, project.name[project.id], subgroup))
sub.set$Relapse.Status <- factor(sub.set$Relapse.Status)
sub.set$Patient.ID <- factor(sub.set$Patient.ID)
col.vec.func <- as.vector(unlist(read.table(file=paste0(folder.path, "/tables/columns_", coverage, ".txt"))))

### ------ load data from database
db.name <- "RB_20191002_Good2018.sqlite3"

fcs$connectDb(file.path(db.path, db.name))
### select project
project.idx <- which(dbListTables(fcs$conn) == project.name[project.id])
fileID <- fcs$getDFtable(paste0(project.name[project.id], "_fileIdentity"))
stainID <- fcs$getDFtable(paste0(project.name[project.id], "_markerIdentity"))

### ------ check if all cutoffs are set for each sampleID
if (check.cutoffs) {
    stopping <- FALSE
    stainID2 <- fcs$getDFtable(paste0(project.name[1], "_markerIdentity"))
    vars <- unique(stainID$shortname)
    for (i in 1:length(vars)) {
        if (unique(stainID[which(stainID$shortname == vars[i]), 5]) != unique(stainID2[which(stainID2$shortname == vars[i]), 5])) {
            printf("cutoff(%s) is not set the same.")
        }
        if (length(unique(stainID[which(stainID$shortname == vars[i]), 5])) > 1) {
        printf("%s has different cutoffs. Set it in the DB first before continueing!", vars[i])
        print(unique(stainID[which(stainID$shortname == vars[i]), 5]))
        printf("1=%s 29=%s",
            stainID[which(stainID$shortname == vars[i])[1], 5],
            stainID[which(stainID$shortname == vars[i])[29], 5])
        stopping <- TRUE
        }
    }
    try(if (stopping) stop("Check your cutoffs!"))
}
dbDisconnect(fcs$conn)

### initiate input dataframe name RDS - file
data.table.name <- sprintf("%s/Rdata/%s/%s_%s_df_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroup, cofactor)
### initiate output dataframe name RDS - file
quad.table.name <- sprintf("%s/Rdata/%s/%s_%s_%s_quadrant_%s_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroup, coverage, stat.info[stat.id], cofactor)

### load matrix temp.data.all
temp.data.all <- readRDS(data.table.name)
print(dim(temp.data.all))
#[1] 3697393      40    for project.id=1 training
#[1] 1943896      40    for project.id=2 training


### - - - - triplots quadrants, cache=TRUE - - - - - - - - - - - - - - - - - -
col.func.idx <- which(names(temp.data.all) %in% col.vec.func)
temp.data.all <- temp.data.all[, c(1, col.func.idx)]
len.var <- ncol(temp.data.all) - 1
colvec <- colnames(temp.data.all)[2:ncol(temp.data.all)]
len.col <- length(colvec)

### ---------- initate
it <- 0
quad.df <- data.frame(matrix(,
nrow <- nrow(sub.set),
ncol <- (len.col * (len.col - 1) * (len.col - 2) * 0.5 * 4 - 1),
#ncol <- 9792, m=18, with ncol <- m * (m - 1) * (m - 2) * 0.5 * 4
byrow <- TRUE
), stringsAsFactors <- FALSE)
quad.sample_id <- vector()
fileConn <- sprintf("tables/2_create_PRI_quadrant_%s_%s.log", project.name[project.id], subgroup)

if (cluster.size > 1) {
    cl <- makeCluster(cluster.size)
    registerDoParallel(cl)
}

ptm <- proc.time()
for (i in 1:nrow(sub.set)) {
# for (i in c(5:7, 28)) {
# for (i in 1:1) {
    quad.sample_id <- c(quad.sample_id, as.character(sub.set$Patient.ID[i]))

    file.name <- paste0(sub.set$Patient.ID[i], "_", project.name.long[project.id], ".fcs")
    file.idx <- fileID$file_ID[which(fileID$filename == file.name)]
    if (length(file.idx) == 0) {
        file.name <- paste0(sub.set$Patient.ID[i], "_", tolower(project.name.long[project.id]), ".fcs") # file names vary..
        file.idx <- fileID$file_ID[which(fileID$filename == file.name)]
    }
    cutoffs <- stainID$file_savedCutoffs[which(stainID$file_ID == file.idx)]
    names(cutoffs) <- stainID$shortname[which(stainID$file_ID == file.idx)]
    # set to same order as in col.func.idx
    cutoffs <- cutoffs[col.func.idx - 1]

    quad.file <- vector()
    for (v1 in 1:(len.col - 1)) {
    # for (v1 in 1:1) {
        it <- it + 1

        quad.oper <- foreach(v2=(v1 + 1):len.col, .combine=cbind) %dopar% {
            #   quad.oper <- foreach(v2=11:14, .combine=cbind) %dopar% {
            quadrant.vec <- vector()

            for (v3 in 1:len.col) {
            #    for (v3 in 13:13) {
                if (all(v3 != c(v1, v2))) {
                sampl.data <- temp.data.all[which(temp.data.all$file_id == as.vector(sub.set$Patient.ID[i])), c(colvec[v1], colvec[v2], colvec[v3])]
                ### NEW::ONLY rows where used if v1 or v2 are >0
                # sampl.data <- sampl.data[which(sampl.data[,1]>=0 & sampl.data[,2]>=0),]

                ### --------- calculate triplot quadrants
                quad.results <- fcs$calc_triplot_quadrant(temp.data = sampl.data, calc.meth = stat.info[stat.id],  min.cells = mincells, prod.cutoff = cutoffs[c(v1, v2)])
                quadrant.vec <- c(quadrant.vec, quad.results)
                }
            }

            return(quadrant.vec)
        }

        quad.file <- c(quad.file, as.vector(quad.oper))

        ### ------ create label vector in the last step
        if (i == nrow(sub.set)) {
            label.file <- vector()
            for (v1 in 1:(len.col - 1)) {
            # for (v1 in 1:1) {

                label.oper <- foreach(v2=(v1 + 1):len.col, .combine=cbind) %dopar% {
                # label.oper <- foreach (v2=11:14,.combine=cbind) %dopar% {
                    label.vec <- vector()

                    for (v3 in 1:len.col) {
                    # for (v3 in 13:13) {
                        if (all(v3 != c(v1, v2))) {

                            if (stat.info[stat.id] == "zRange") {
                                label.vec <- c(label.vec, paste0(colvec[v1], ".", colvec[v2], ".", colvec[v3]))
                            } else {
                                label.vec <- c(label.vec,
                                    paste0(colvec[v1], ".", colvec[v2], ".", colvec[v3], ".", stat.info[stat.id ], ".", "Q1"),
                                    paste0(colvec[v1], ".", colvec[v2], ".", colvec[v3], ".", stat.info[stat.id ], ".", "Q2"),
                                    paste0(colvec[v1], ".", colvec[v2], ".", colvec[v3], ".", stat.info[stat.id ], ".", "Q3"),
                                    paste0(colvec[v1], ".", colvec[v2], ".", colvec[v3], ".", stat.info[stat.id ], ".", "Q4")
                                    )
                            }
                        }
                    }
                    return(label.vec)
                }
                label.file <- c(label.file, as.vector(label.oper))
            }
        }

        printf("%s::%s::%s::%s::%s::v1=%s[%s/%s] ready (it=%s)", i, sub.set$Patient.ID[i], project.name[project.id], subgroup, stat.info[stat.id], colvec[v1], v1, len.col, it)
    
        print(proc.time() - ptm)
    }
    quad.df[i, ] <- quad.file

    printf("File %s ready (it=%s)", i, it)
    print(proc.time() - ptm)

    sink(fileConn, append = TRUE)
    cat(paste0(date(), "\n"), append = TRUE)
    cat(sprintf("%s::%s::%s::%s::%s::cluster=%s \n", i, sub.set$Patient.ID[i], project.name[project.id], subgroup, coverage, cluster.size), append = TRUE)
    cat(paste0(proc.time() - ptm, collapse = " "), append = TRUE)
    cat("\n", append = TRUE)
    if (i == nrow(sub.set)) cat(paste0("Done. ", quad.table.name, " saved!"))
    sink()
}
print("Total time:")
print(proc.time() - ptm)

colnames(quad.df) <- label.file
rownames(quad.df) <- quad.sample_id

saveRDS(quad.df, file=quad.table.name)
printf("%s saved!", quad.table.name)

### this followed because some files were missing in training data
if (FALSE) {
    print("Still need to bind_rows with other training data!")
    ### fastest way to bind tables with rbindlist()
    library(data.table)
    quad.df2 <- quad.df
    quad.df1 <- readRDS(sprintf("%s/Rdata/%s/%s_%s_%s_quadrant_%s_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroup, coverage, stat.info[stat.id], cofactor))
    quad.df.all <- rbindlist(list(quad.df1, quad.df2))
    saveRDS(quad.df.all, file = sprintf("%s/Rdata/%s/%s_%s_%s_quadrant.all_%s_cof%s.rds", folder.path, project.name[project.id], project.name[project.id], subgroup, coverage, stat.info[stat.id], cofactor))
    printf(sprintf("%s/Rdata/%s/%s_%s_%s_quadrant.all_%s_cof%s.rds saved!", folder.path, project.name[project.id], project.name[project.id], subgroup, coverage, stat.info[stat.id], cofactor))
}

if (cluster.size > 1) stopCluster(cl)