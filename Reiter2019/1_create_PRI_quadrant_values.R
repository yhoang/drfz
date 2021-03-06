#!/usr/bin/R
# Authors: Yen Hoang, Felix Lohrke and Max Salomon
# DRFZ 2020-2021
# Reiter2019


rm(list = ls())
options(max.print = 100)

### ----------- initiate
# initials      "YH", "FL", "MS" 
# work.station  set correct paths
# load.from.DB  TRUE/FALSE, to load data from Sqlite db or from RDS file (merge from script 1+2)
# cofactor      1, 0.2, 0.1 (so far at 0.2)
# mincells      2,5,10,20 (so far at 5)
# stat.id       1:4 for absRange, variance, freq.green, mean
# cluster.size  1:12, detect with detectCores()
# db.id         1:3 for VIE_Routine, BUE_Dura, BLN_Dura
# subgroup      "Total", "Validation", "Training"
# subset        full
# comment       "autoSec" for automatic cutoffs, "manSec" for manual cutoffs set in database
# check.cutoffs TRUE/FALSE, check if cutoffs are set in all files
initials <- "YH"
work.station <- "delta"
load.from.DB <- TRUE
cofactor <- 0.2
mincells <- 5
cluster.size <- 1
db.id <- 1
subgroup <- "Total"

subset <- "full"
subject <- "quadrant"
stat.id <- 1
comment <- "autoSec"
check.cutoffs <- FALSE

### set paths
if (work.station == "asus-vivid") {
  library.dir = file.path("C:", "Program Files", "R", "R-3.6.1", "library")
  main.dir <- file.path("D:", "drfz", "Reiter2019")
  db.dir <- file.path("D:", "DB")
} else if (work.station == "delta") {
  library.dir <- file.path("usr","lib","R","library")
  main.dir <- file.path("/scratch", "drfz", "Reiter2019")
  db.dir <- file.path("/data", "databases")
}
setwd(main.dir)
################



### --------- load libraries
# xlsx      :  Library for Excel reading and creating
# RSQLite     :  interact with database
# reshape2    :  function melt()
# dplyr      :  faster binding of columns bind_cols() and rows bind_rows()
# foreach     :  allows for each loops
# doParallel   :  use several clusters for parallele calculations
.libPaths(library.dir)
load.libraries <- c("RSQLite", "dplyr", "reshape2", "foreach", "doParallel", "stringr")
#install.packages(libraries, lib = library.dir)
lapply(load.libraries, require, character.only = TRUE)
### load functions
source(file.path(main.dir, "PRI_funs.R"))
printf <- function(...) invisible(print(sprintf(...)))

dataset.name <- c("VIE_Routine", "BUE_Dura", "BLN_Dura")
stat.info <- c("absRange", "variance", "freq.green", "mean")
today <- paste0(initials, substring(str_replace_all(Sys.Date(), "-", ""), 3))

### NO METADATA YET -----------------------------------------
### NOTE: db VIE has one marker less - might affect downstream analysis!
col.vec.func <- as.vector(unlist(read.table(file=paste0(main.dir, "/tables/columns_", dataset.name[db.id], "_", subset, ".txt"))))
len.col <- length(col.vec.func)

####### load data base -----------------------------------------
db.name <- sprintf("FL_20201112_Reiter-%s.sqlite3", dataset.name[db.id])
fcs$connectDb(file.path(db.dir, db.name))
### select project
db.idx <- which(dbListTables(fcs$conn) == dataset.name[db.id])
fileID <- fcs$getDFtable(paste0(dataset.name[db.id], "_fileIdentity"))
stainID <- fcs$getDFtable(paste0(dataset.name[db.id], "_markerIdentity"))
### create folders if not existent
sub.dir <- c("Rdata", "tables", "log")
dir.create(file.path(main.dir, sub.dir[1]), showWarnings = FALSE)
dir.create(file.path(main.dir, sub.dir[1], dataset.name[db.id]), showWarnings = FALSE)
dir.create(file.path(main.dir, sub.dir[2]), showWarnings = FALSE)
dir.create(file.path(main.dir, sub.dir[3]), showWarnings = FALSE)

### ------ check if all cutoffs are set for each sampleID
if (check.cutoffs) {
    stopping <- FALSE
    stainID2 <- fcs$getDFtable(paste0(dataset.name[db.id], "_markerIdentity"))
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

############### SET OUTPUT FILE NAMES
# dataframe name RDS - file
quad.table.name <- sprintf("%s/%s/%s/%s_%s_%s_%s_%s_%s.cof%s_%s.rds", main.dir, sub.dir[1], dataset.name[db.id], dataset.name[db.id], subgroup, subset, subject, stat.info[stat.id], comment, cofactor, today)
# log file
log.File <- sprintf("%s/%s/2_create_PRI_%s_%s_%s_%s.log", main.dir, sub.dir[3], subject, dataset.name[db.id], subgroup, today)
################

if (!load.from.DB) {
    ### INPUT dataframe name RDS - file
    data.table.name <- sprintf("%s/%s/%s/%s_%s_df_cof%s.rds", main.dir, sub.dir[1], dataset.name[db.id], dataset.name[db.id], subgroup, cofactor)
    ### load matrix temp.data.all
    temp.data.all <- readRDS(data.table.name)
    temp.data.all <- as.data.frame(temp.data.all)
    print(dim(temp.data.all))

    ### removing rows with NAs
    if (sum(is.na(temp.data.all)) > 0) {
        printf("%s NAs in all data from %s. Removing..", sum(is.na(temp.data.all)), dataset.name[db.id])
        sink(log.File, append = TRUE)
        cat(sprintf("%s NAs in temp.data.all. Removing..", sum(is.na(temp.data.all))), append = TRUE)
        sink()
    }
    temp.data.all <- na.omit(temp.data.all)
}

### - - - - triplots quadrants, cache=TRUE - - - - - - - - - - - - - - - - - -

### ---------- initate
sub.set <- fileID[,2]
it <- 0
quad.df <- data.frame(matrix(,
    nrow <- length(sub.set),
    ncol <- (len.col * (len.col - 1) * (len.col - 2) * 0.5 * 4 - 1),
    byrow <- TRUE
    ), stringsAsFactors <- FALSE
)
quad.sample_id <- vector()

### turn on clusters
cl <- makeCluster(cluster.size)
registerDoParallel(cl)
ptm <- proc.time()

for (i in 1:length(sub.set)) {
# for (i in 143:143) {

    if (load.from.DB) {
        ### access data from Sqlite database
        pat.id <- file.name <- sub.set[i]
        if (db.id == 1) {pat.id <- substr(file.name, 1, 10)
        } else { pat.id <- substr(file.name, 1, 6) }
        printf("%s/%s::Looking for file %s..", i, length(sub.set), file.name)
        file.idx <- fileID$file_ID[which(fileID$filename == file.name)]
        ### get data with cofactor
        temp.data.all <- fcs$getData(table = dataset.name[db.id], fileidx = file.idx, cofactor=cofactor)
        
        ### set negative expressions to zero
        temp.data.all[temp.data.all < 0] <- 0
        ### removing rows with NAs
        if (sum(is.na(temp.data.all)) > 0) {
            printf("%s NAs in %s. Removing..", sum(is.na(temp.data.all)), pat.id)
            sink(log.File, append = TRUE)
            cat(sprintf("%s NAs in %s. Removing..", sum(is.na(temp.data.all)), pat.id), append = TRUE)
            sink()
        }
        temp.data.all <- na.omit(temp.data.all)
        
        ### get short names
        colnames(temp.data.all) <- stainID$shortname[which(stainID$file_ID == file.idx)]
        
        col.fileid <- rep(unlist(pat.id), nrow(temp.data.all))
        col.func.idx <- which(names(temp.data.all) %in% col.vec.func)
        temp.data.all <- as.data.frame(bind_cols(file_id=col.fileid, temp.data.all[, col.func.idx]))
        len.var <- ncol(temp.data.all) - 1
        colvec <- colnames(temp.data.all)[2:ncol(temp.data.all)]
    }

    quad.sample_id <- c(quad.sample_id, as.character(sub.set[i]))

    if (comment != "autoSec") {
        file.idx <- fileID$file_ID[which(fileID$filename == sub.set[i])]
        cutoffs <- stainID$file_savedCutoffs[which(stainID$file_ID == file.idx)]
        names(cutoffs) <- stainID$shortname[which(stainID$file_ID == file.idx)]
        # set to same order as in col.func.idx
        cutoffs <- cutoffs[col.func.idx - 1]  
    }

    quad.file <- vector()
    for (v1 in 1:(len.col - 1)) {
    # for (v1 in 1:1) {
        it <- it + 1

        quad.oper <- foreach(v2=(v1 + 1):len.col, .combine=cbind) %dopar% {
        # quad.oper <- foreach(v2=(v1 + 1), .combine=cbind) %dopar% {
            quadrant.vec <- vector()

            for (v3 in 1:len.col) {
            # for (v3 in 1:1) {
                if (all(v3 != c(v1, v2))) {

                    sampl.data <- temp.data.all[ ,c(colvec[v1], colvec[v2], colvec[v3])]

                    ### --------- calculate triplot quadrants
                    if (comment == "autoSec") {
                        quad.results <- fcs$calc_triplot_quadrant(temp.data = sampl.data, calc.meth = stat.info[stat.id], min.cells = mincells)
                    } else {
                      quad.results <- fcs$calc_triplot_quadrant(temp.data = sampl.data, calc.meth = stat.info[stat.id], min.cells = mincells, prod.cutoff = cutoffs[c(v1, v2)])
                    }
                    quadrant.vec <- c(quadrant.vec, quad.results)
                }
            }

            return(quadrant.vec)
        }

        quad.file <- c(quad.file, as.vector(quad.oper))

        ### ------ create label vector in the last step
        if (i == length(sub.set)) {
            label.file <- vector()
            for (v1 in 1:(len.col - 1)) {

                label.oper <- foreach(v2=(v1 + 1):len.col, .combine=cbind) %dopar% {
                    label.vec <- vector()

                    for (v3 in 1:len.col) {
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

        if (it %% 10 == 0) {
            printf("%s::%s::%s::%s::%s::%s::v1=%s[%s/%s] ready (it=%s)", i, sub.set[i], dataset.name[db.id], subgroup, stat.info[stat.id], comment, colvec[v1], v1, len.col, it)
            print(proc.time() - ptm)
        }
    }
    quad.df[i, ] <- quad.file

    printf("File %s ready (it=%s)", i, it)
    print(proc.time() - ptm)

    sink(log.File, append = TRUE)
    cat(paste0(date(), "\n"), append = TRUE)
    cat(sprintf("%s::%s::%s::%s::%s::cluster=%s \n", i, sub.set[i], dataset.name[db.id], subgroup, subset, cluster.size), append = TRUE)
    cat(paste0(proc.time() - ptm, collapse = " "), append = TRUE)
    cat("\n", append = TRUE)
    if (i == length(sub.set)) cat(paste0("Done. ", quad.table.name, " saved!"))
    sink()
}
print("Total time:")
print(proc.time() - ptm)

colnames(quad.df) <- label.file
rownames(quad.df) <- quad.sample_id

saveRDS(quad.df, file=quad.table.name)
printf("%s saved!", quad.table.name)

dbDisconnect(fcs$conn)

stopCluster(cl)

