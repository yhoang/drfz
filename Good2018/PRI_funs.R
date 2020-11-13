#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019-2020
# Goood2018

fcs <- new.env()

### helpful print function
printf <- function(...) invisible(print(sprintf(...)))

is.nan.data.frame <- function(x) {
  do.call(cbind, lapply(x, is.nan))
}


fcs$connectDb <- function(fname) {
  this <- fcs
  this$conn <- dbConnect(SQLite(), dbname <- fname)
  print(paste("Database opened:", fname))
}

fcs$disconnectDb <- function() {
  disconnectDb(fcs$conn)
}

### get data from file index in database
fcs$getData <- function (table, fileidx, columns=NA, stain=NA, cofactor = NA) {
  this <- fcs

  data <- dbGetQuery(this$conn, paste("SELECT * FROM ", table, " WHERE file_ID == '", fileidx, "'", sep=""))

  # and ignore columns with NAs
  col.NA <- NA
  for (i in 1:ncol(data)) {
    if (any(is.na(data[, i])))  col.NA <- c(col.NA, i)
  }
  if (!is.na(col.NA)) data <- data[, -col.NA]

  ### change column names if provided
  column.names <- colnames(data)
  if (!is.na(columns)) {
    data <- data[, columns + 1]
  }

  if (!is.na(stain) & !is.na(columns)) colnames(data) <- stain[columns]
  if (is.na(stain) & !is.na(columns)) colnames(data) <- this$selected.vars[columns]
  if (!is.na(stain) & is.na(columns)) colnames(data) <- c("file_ID", stain)
  if (is.na(stain) & is.na(columns)) colnames(data) <- c("file_ID", this$selected.vars)

  # set asinh cofactor to 1 if not set
  if (is.na(cofactor)) cofactor <- 1

  if (!is.na(columns)) {
    data <- asinh(data / cofactor)
  } else {
    data <- asinh(data[, (2:dim(data)[2])] / cofactor)
  }

  printf("w: do getData(%s) from table='%s' with fileidx=%s and asinh cofactor=%s", nrow(data), table, fileidx, cofactor)

  data
}

### get full table from database
fcs$getDFtable <- function (table) {
  this <- fcs

  table.df <- dbGetQuery(this$conn, paste("SELECT * FROM ", table))

  return(table.df)
}


######## LOAD FUNCTIONS WITHOUT ACTUALLY CALLING THEM ##############
##### @param
# project.idx   project index
# file.idx      file index
##### returns data table
fcs$loaddata <- function(
  project.idx,
  file.idx){
  this <- fcs

  ### list all tables in database
  this$table.lis <- dbListTables(this$conn)

  ### get every index which has metadata
  metadata.idx <- grep("markerIdentity|colnameIndex|fileIdentity|fileIndex|UserHistory|Classification|equipmentInfo|fileComments|SPILL", this$table.list)

  this$metadata.list <- vector()
  df.num <- 0
  for (i in metadata.idx){
    df.num <- df.num + 1
    this$metadata.list[df.num] <- this$table.list[i]
  }

  ### now only the table names are listed which have our data (fluorescence intensity for each feature)
  this$project.list <- this$table.list[-metadata.idx]
  ### if there are several projects in the database, choose one
  this$current.project <- this$project.list[project.idx]

  this$current.filetable <- this$getDFtable(paste0(this$current.project, "_fileIdentity"))
  this$file.list <- this$current.filetable[,  2]

  this$current.file <- this$file.list[file.idx]
  this$current.staintable <- this$getDFtable(paste0(this$current.project, "_markerIdentity"))
  this$current.vars <- this$getVariables(index=file.idx)

  data <- this$getFile(this$current.project, file.idx)

  data
}


######## CREATE TRIPLOT CALCULATION fcs$my.calc w/o PLOTTING ##############
##### @param
# data              data matrix to analyse
# stainTable        table from DB, information of column names and cutoffs
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# calc              calculation method
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(1,12,1,12)
fcs$bintriploT_table <- function(
  data,
  stainTable,
  feat.X,
  feat.Y,
  feat.Z1,
  calc,
  binsize,
  mincells,
  plot.range=c(1, 12, 1, 12)) {
  this <- fcs

  if (FALSE) {
    feat.X <- "PD-1"
    feat.Y <- "IFNg"
    feat.Z1 <- "IL21"
    calc <- "MFI+"
    binsize <- 0.2
    mincells <- 10
  }

  ### remove all signs and write anything with capitals
  feat.X.clean <- gsub("[^[:alnum:]]", "", toupper(feat.X))
  feat.Y.clean <- gsub("[^[:alnum:]]", "", toupper(feat.Y))
  feat.Z1.clean <- gsub("[^[:alnum:]]", "", toupper(feat.Z1))

  ### axes range
  xmin.val <- plot.range[1]
  xmax.val <- plot.range[2]
  ymin.val <- plot.range[3]
  ymax.val <- plot.range[4]

  ############ LOAD DATA AND CUTOFFS
  cutoffs <- as.numeric(stainTable[which(stainTable$file_ID == file.idx), 5])

  features <- colnames(data)

  ### remove all signs and write anything with capitals
  features.clean <- gsub("[^[:alnum:]]", "", make.unique(unlist(lapply(features, function(x) {
    len <- length(strsplit(x, "[.]")[[1]])
    y <- toupper(strsplit(x, "[.]")[[1]][1])
    paste(y, collapse=".")
  }))))

  idx.X <- which(features.clean == feat.X.clean)
  idx.Y <- which(features.clean == feat.Y.clean)
  idx.Z1 <- which(features.clean == feat.Z1.clean)

  if (is.na(idx.X) | is.na(idx.Y) | is.na(idx.Z1)) stop(sprintf("Either feat.X (%s) or feat.Y (%s) or feat.Z1 (%s) not in file index %s", feat.X, feat.Y, feat.Z1, i))

  ### cut data if calc method is MFI+
  if (calc == "MFI+") data <- data[which(data[, idx.Z1] > cutoffs[idx.Z1]), ]

  ### construct bin table with number of cells per bin
  fX <- cut(data[, idx.X], breaks=seq(xmin.val, xmax.val, by=binsize), include.lowest=TRUE, dig.lab=5)
  fY <- cut(data[, idx.Y], breaks=seq(ymin.val, ymax.val, by=binsize), include.lowest=TRUE, dig.lab=5)
  fXY <- as.factor(paste(fX, fY))

  tab <- table(fX, fY)
  colnames(tab) <- seq(ymin.val, ymax.val - binsize, by=binsize)
  rownames(tab) <- seq(xmin.val, xmax.val - binsize, by=binsize)

  this$tab <- tab

  # printf("length(bins with tab >= mincells)=%s",length(which(tab >= mincells)))

  ### set negative values of z-axis (colnum=3) to zero
  data[which(data[, idx.Z1] < 0), idx.Z1] <- 0

  if (grepl(calc, "MFI")) {
    ########## CALCULATE MFI
    my.lengths <- aggregate(data[, idx.Z1], by=list(fXY), length)
    idx.len <- which(my.lengths$x >= mincells)

    ### start with MFI calculation of Z1
    my.calc <- aggregate(data[, idx.Z1], by=list(fXY), mean)

    min.range.Z1 <- floor(min(my.calc[idx.len, "x"]) * 10) / 10
    max.range.Z1 <- ceiling(max(my.calc[idx.len, "x"]) * 10) / 10

    # get steps for Z1
    step <- round(diff(range(max.range.Z1, min.range.Z1)) / 10, 2)
    steps.Z1 <- seq(min.range.Z1, max.range.Z1, by=step)

    # bin color factor Z1
    my.calc.fac.Z1 <- cut(my.calc$x, breaks=steps.Z1, labels=1:10, include.lowest=TRUE)
    names(my.calc.fac.Z1) <- my.calc$x

    ### combine all frequencies in one table
    my.calc <- cbind(my.calc, fac.Z1=as.numeric(my.calc.fac.Z1))
    my.calc <- cbind(my.calc, ncells=my.lengths$x)

  } else if (calc == "freq") {
    ########## CALCULATE FREQUENCIES
    ### frequency of feature Z1
    my.calc <- aggregate(data[, idx.Z1], by=list(fXY), function(x) {
      y <- round(100 * length(which(x >= cutoffs[idx.Z1])) / length(x))
      return(y)
    })

    # bin color factor Z1
    my.calc.fac.Z1 <- cut(my.calc$x, breaks=seq(0, 100, by=10), labels=1:10, include.lowest=TRUE)
    names(my.calc.fac.Z1) <- my.calc$x

    my.lengths <- aggregate(data[, idx.Z1], by=list(fXY), length)

    ### combine all frequencies in one table
    my.calc <- cbind(my.calc, fac.Z1=as.numeric(my.calc.fac.Z1))
    my.calc <- cbind(my.calc, ncells=my.lengths$x)
  }
  this$my.calc <- my.calc
  this$my.calc.fac.Z1 <- my.calc.fac.Z1
  ########## DONE CALCULATE FREQUENCIES

}


######## CREATE DIPLOT CALCULATION fcs$my.calc w/o PLOTTING ##############
##### @param
# data              data matrix to analyse
# stainTable        table from DB, information of column names and cutoffs
# feat.X            name of feature X
# feat.Y1           name of feature Y1
# calc              calculation method
# binsize           bin size
# mincells          minimum amount of cells in a bin
# plot.range        plot range, first x axis, second y axis, default= c(1,12)
fcs$bindiploT_table <- function(
  data,
  stainTable,
  feat.X,
  feat.Y1,
  calc,
  binsize,
  mincells,
  plot.range=c(1, 12)){
  this <- fcs

  if (FALSE) {
    feat.X <- "PD-1"
    feat.Y <- "IFNg"
    feat.Z1 <- "IL21"
    calc <- "MFI+"
    binsize <- 0.2
    mincells <- 10
  }

  ### remove all signs and write anything with capitals
  feat.X.clean <- gsub("[^[:alnum:]]", "", toupper(feat.X))
  feat.Y1.clean <- gsub("[^[:alnum:]]", "", toupper(feat.Y1))

  ### axes range
  xmin.val <- plot.range[1]
  xmax.val <- plot.range[2]

  ############ LOAD DATA AND CUTOFFS
  cutoffs <- as.numeric(stainTable[which(stainTable$file_ID == file.idx), 5])

  features <- colnames(data)

  ### remove all signs and write anything with capitals
  features.clean <- gsub("[^[:alnum:]]", "", make.unique(unlist(lapply(features, function(x) {
    len <- length(strsplit(x, "[.]")[[1]])
    y <- toupper(strsplit(x, "[.]")[[1]][1])
    paste(y, collapse=".")
  }))))

  idx.X <- which(features.clean == feat.X.clean)
  idx.Y1 <- which(features.clean == feat.Y1.clean)

  if (is.na(idx.X) | is.na(idx.Y1)) stop(sprintf("Either feat.X (%s) or feat.Y1 (%s) not in file index %s", feat.X, feat.Y1, i))

  ### cut data if calc method is MFI+
  if (calc == "MFI+") data <- data[which(data[, idx.Y1] > cutoffs[idx.Y1]), ]

  ### construct bin table with number of cells per bin
  fX <- cut(data[, idx.X], breaks=seq(xmin.val, xmax.val, by=binsize), include.lowest=TRUE, dig.lab=5)

  ## set negative values of z-axis (colnum=3) to zero
  data[which(data[, idx.Y1] < 0), idx.Y1] <- 0

  ### calculate cells in bins
  my.lengths <- aggregate(data[, idx.Y1], by=list(fX), length)
  idx.len <- which(my.lengths$x >= mincells)

  if (grepl(calc, "MFI")) {
    ########## CALCULATE MSI
    ### MSI calculation of Y1
    my.calc <- aggregate(data[, idx.Y1], by=list(fX), mean)

    min.range.Y1 <- floor(min(my.calc[idx.len, "x"]) * 10) / 10
    max.range.Y1 <- ceiling(max(my.calc[idx.len, "x"]) * 10) / 10

    # get bin steps for Y1
    step <- round(diff(range(max.range.Y1, min.range.Y1)) / 10, 2)
    steps.Y1 <- seq(min.range.Y1, max.range.Y1, by=step)

    # bin color factor Y1
    my.calc.fac.Y1 <- cut(my.calc$x, breaks=steps.Y1, labels=1:10, include.lowest=TRUE)
    names(my.calc.fac.Y1) <- my.calc$x
    ########## DONE CALCULATE MSI

  } else if (calc == "freq") {
    ########## CALCULATE FREQUENCIES
    ### frequency of feature Y1
    my.calc <- aggregate(data[, idx.Y1], by=list(fX), function(x) {
      y <- round(100 * length(which(x >= cutoffs[idx.Y1])) / length(x))
      return(y)
    })

    # bin color factor Y1 always from 0 to 100
    my.calc.fac.Y1 <- cut(my.calc$x, breaks=seq(0, 100, by=10), labels=1:10, include.lowest=TRUE)
    names(my.calc.fac.Y1) <- my.calc$x
    ########## DONE CALCULATE FREQUENCIES
    #my.lengths <- aggregate(data[, idx.Y1], by=list(fX), length)
    ### combine all frequencies in one table
    # my.calc <- cbind(my.calc, fac.Y1=as.numeric(my.calc.fac.Y1))
    # my.calc <- cbind(my.calc, ncells=my.lengths$x)
  }
  ### combine all values in one table
  my.calc <- cbind(my.calc, fac.Y1=as.numeric(my.calc.fac.Y1))
  my.calc <- cbind(my.calc, ncells=my.lengths$x)

  ### globalize
  this$my.calc <- my.calc
  this$my.calc.fac.Y1 <- my.calc.fac.Y1
}


##### Calculate triplot quadrants function
##### @param
# temp.data     temporary data with 3 columns v1,v2,v3
# calc.meth     calculation method to use [mean,sd,median,zrange,absRange,relRange,var]
# prod.cutoff   cutoff to divide neg/pos cells, here to divide quadrants. Be aware! Changes with different cofactors for asinh transformation
# size.bin      bin size
# min.cells     minimum number of cells per bin to consider this a countable bin
# min.bin.count if quadrant does not have enough bins, then exclude. default=NA
fcs$calc_triplot_quadrant <- function (
  temp.data,
  calc.meth,
  prod.cutoff = NA,
  size.bin = 0.2,
  min.cells = 5,
  min.bin.count = NA){

  ### to start manually
  if (FALSE) {
    temp.data <- sampl.data
    calc.meth <- stat.info[stat.id]
    min.cells <- mincells
    prod.cutoff <- cutoffs[c(v1, v2)]
    size.bin <- 0.2
    min.bin.count <- NA
  }

  brackets.open <- c("(", "[")

  ### min value for x and y axis
  min.x <- floor(min(temp.data[,  1]) * 10) / 10
  min.y <- floor(min(temp.data[,  2]) * 10) / 10

  ### max value for x and y axis
  max.x <- ceiling(max(temp.data[,  1]) * 10) / 10
  max.y <- ceiling(max(temp.data[,  2]) * 10) / 10

  ### ranges need to be even!
  # only with bin size <- 0.2
  is.even <- function(x) x %% 2 == 0
  if (size.bin == 0.2) {
    if (!is.even(min.x * 10)) min.x <- min.x - 0.1
    if (!is.even(min.y * 10)) min.y <- min.y - 0.1
    if (!is.even(max.x * 10)) max.x <- max.x + 0.1
    if (!is.even(max.y * 10)) max.y <- max.y + 0.1
  }

  temp.data <- as.data.frame(temp.data)

  ### initialize bin construct
  fX <- cut(temp.data[,  1], breaks=seq(min.x, max.x, by=size.bin), include.lowest=TRUE, dig.lab=5)
  fY <- cut(temp.data[,  2], breaks=seq(min.y, max.y, by=size.bin), include.lowest=TRUE, dig.lab=5)
  tab <- table(fX, fY)
  seq.bin.y <- seq(min.y, max.y - size.bin, by=size.bin)
  seq.bin.x <- seq(min.x, max.x - size.bin, by=size.bin)
  colnames(tab) <- seq(min.y, max.y - size.bin, by=size.bin)
  rownames(tab) <- seq(min.x, max.x - size.bin, by=size.bin)

  ### ONLY if minimum bin count is set
  ### if triploT has not enough bins, quadrant values will be NA, otherwise continue
  if (!is.na(min.bin.count) & (length(which(tab > min.cells)) < min.bin.count)) {
    vec.q1 <- vec.q2 <- vec.q3 <- vec.q4 <- NA
  } else {

    ############# if cutoff is set
    if (all(!is.na(prod.cutoff))) {
      mean.bin.x <- prod.cutoff[1]
      mean.bin.y <- prod.cutoff[2]
    } else {
      ############# if cutoff is not set
      ########### mean values for x and y where bins are displayed to divide into 4 quadrants

      ### find minimum and maximum x where bins are displayed
      ### max.bin.y
      for (max.bin.y in ncol(tab):1) {
        if (any(tab[, max.bin.y] >= min.cells)) break
      }
      if (max.bin.y < length(seq.bin.y)) {
        max.bin.y <- seq.bin.y[max.bin.y + 1]
      } else {
        max.bin.y <- max.y
      }
      ### min.bin.y
      for (min.bin.y in 1:ncol(tab)) {
        if (any(tab[, min.bin.y] >= min.cells)) break
      }
      min.bin.y <- seq.bin.y[min.bin.y]
      ############ find minimum and maximum y where bins are displayed
      ### max.bin.x
      for (max.bin.x in nrow(tab):1) {
        if (any(tab[max.bin.x, ] >= min.cells)) break
      }
      if (max.bin.x < length(seq.bin.x)) {
        max.bin.x <- seq.bin.x[max.bin.x + 1]
      } else {
        max.bin.x <- max.x
      }
      ### min.bin.x
      for (min.bin.x in 1:nrow(tab)) {
        if (any(tab[min.bin.x, ] >= min.cells)) break
      }
      min.bin.x <- seq.bin.x[min.bin.x]

      # mean.bin.x <- prod.cutoff[1]
      # mean.bin.y <- prod.cutoff[2]
      mean.bin.x <- round(range(min.bin.x, max.bin.x) / 2, 1)
      mean.bin.y <- round(range(min.bin.y, max.bin.y) / 2, 1)

      ### make NEW bin construct with NEW minimum and maximum x/y where bins are displayed
      fX <- cut(temp.data[,  1], breaks=seq(min.bin.x, max.bin.x, by=size.bin), include.lowest=TRUE, dig.lab=5)
      fY <- cut(temp.data[,  2], breaks=seq(min.bin.y, max.bin.y, by=size.bin), include.lowest=TRUE, dig.lab=5)
      tab <- table(fX, fY)
      colnames(tab) <- seq(min.bin.y, max.bin.y - size.bin, by=size.bin)
      rownames(tab) <- seq(min.bin.x, max.bin.x - size.bin, by=size.bin)
    }

    ### combine X and Y position together
    fXY <- as.factor(paste(fX, fY))

    ### group cells which calls into bin together and calculate bin method
    my.lengths <- aggregate(temp.data[, 3], by=list(fXY), length)

    if (calc.meth == "median") {
      my.calc <- aggregate(temp.data[, 3], by=list(fXY), median)
    } else if (calc.meth == "sd") {
      my.calc <- aggregate(temp.data[, 3], by=list(fXY), sd)
    } else if (calc.meth == "zrange") {
      my.calc <- aggregate(temp.data[, 3], by=list(fXY), function(z) {
        return(diff(range(max(z), min(z))))
      })
    } else if (calc.meth == "variance") {
      my.calc <- aggregate(temp.data[, 3], by=list(fXY), var)
    } else {
      # if (calc.meth == "mean" | calc.meth == "range.quad" | calc.meth == "range.quad_scaled")
      my.calc <- aggregate(temp.data[, 3], by=list(fXY), mean)
    }
    my.calc <- cbind(my.calc, bin.cells <- my.lengths$x)

    ### only use with different calc.meth then variance or absRange
    if (!(calc.meth %in% c("variance", "absRange"))) {
      ### bin positions where we have at least minimum cells
      idx.min.cells <- which(my.calc$bin.cells >= min.cells)
      ### min/max of range var colvec[v3]
      min.v3.range <- floor(min(my.calc[idx.min.cells, "x"]) * 10) / 10
      max.v3.range <- ceiling(max(my.calc[idx.min.cells, "x"]) * 10) / 10
      ### get bin steps of range var colvec[v3]
      step <- round(diff(range(max.v3.range, min.v3.range)) / 10, 2)
      steps <- seq(min.v3.range, max.v3.range, by=step)
      ### bin factor
      ### e.g. mean here is scaled mean (bin factor mean from 0-9)
      my.calc.fac <- cut(my.calc$x, breaks=steps, labels=0:9, include.lowest=TRUE)

      ### combine length and bin factors
      my.calc <- cbind(my.calc, fac=as.numeric(my.calc.fac) - 1)
    }

    ### initialize quadrant vectors
    vec.q1 <- vec.q2 <- vec.q3 <- vec.q4 <- vector()

    ### divide bin structure into quadrants
    for (i in 1:length(rownames(tab))) {
      for (j in 1:length(colnames(tab))) {
        ### only use bins with enough cells
        x <- rownames(tab)[i]
        y <- colnames(tab)[j]
        if (tab[x, y] >= min.cells) {
          ### create factor
          brackets.idx.x <- brackets.idx.y <- 1
          if (i == 1) brackets.idx.x <- 2
          if (j == 1) brackets.idx.y <- 2

          fact <- as.factor(paste(brackets.open[brackets.idx.x], x, ",",
                      as.numeric(x) + size.bin, "] ",
                      brackets.open[brackets.idx.y], y, ",",
                      as.numeric(y) + size.bin, "]", sep=""))

          idx.bin <- which(as.character(fact) == as.character(my.calc$Group.1))

          if (calc.meth  %in% c("var", "absRange")) {
            fac.bin <- my.calc[idx.bin, "x"]
          } else {
            fac.bin <- my.calc[idx.bin, "fac"]
          }

          if ((as.numeric(x) < mean.bin.x) & (as.numeric(y) < mean.bin.y)) {
            ### Q1
            vec.q1 <- c(vec.q1, fac.bin)
          } else if ((as.numeric(x) < mean.bin.x) & (as.numeric(y) >= mean.bin.y)) {
            ### Q2
            vec.q2 <- c(vec.q2, fac.bin)
          } else if ((as.numeric(x) >= mean.bin.x) & (as.numeric(y) >= mean.bin.y)) {
            ### Q3
            vec.q3 <- c(vec.q3, fac.bin)
          } else if ((as.numeric(x) >= mean.bin.x) & (as.numeric(y) < mean.bin.y)) {
            ### Q4
            vec.q4 <- c(vec.q4, fac.bin)
          }
        }
      }
    }
  }

  ### return range of every quadrants listed bins with minimum cells and used calculation method
  if (calc.meth == "absRange") {
    return(c(round(diff(range(vec.q1)), 3),
              round(diff(range(vec.q2)), 3),
              round(diff(range(vec.q3)), 3),
              round(diff(range(vec.q4)), 3)
    ))
    ### return variance of every quadrants listed bins with minimum cells and used calculation method
  } else if (calc.meth == "variance") {
    return(c(round(mean(vec.q1), 3),
              round(var(vec.q2), 3),
              round(var(vec.q3), 3),
              round(var(vec.q4), 3)
    ))
    ### return mean of every quadrants listed bins with minimum cells and used calculation method
  } else {
    return(c(round(mean(vec.q1), 3),
              round(mean(vec.q2), 3),
              round(mean(vec.q3), 3),
              round(mean(vec.q4), 3)
    ))
  }

}


##### Calculate frequency in gates function
##### @param
# temp.data     temporary data with 3 columns v1,v2,v3
# prod.cutoff   cutoff to divide neg/pos cells, here to divide quadrants. Be aware! Changes with different cofactors for asinh transformation
# output        frequency output (default="green") [black,red,green]
# MSIplus       binary option, if data should be cutted first with feature C (default=FALSE)
fcs$calc_frequencies <- function (
  temp.data,
  prod.cutoff = NA,
  output = "green",
  MSIplus = FALSE){

  brackets.open <- c("(", "[")

  temp.data <- as.data.frame(temp.data)

  ### calc quadrants in total
  ncells <- ncells.total <- nrow(temp.data)
  # q1 quadrant unten links
  # q2 quadrant unten rechts
  # q3 quadrant oben rechts
  # q4 quadrant oben links

  ### count cells in quadrant
  temp.data.q1 <- temp.data[which(temp.data[,  1] < prod.cutoff[1] &  temp.data[,  2] < prod.cutoff[2]), 3]
  temp.data.q2 <- temp.data[which(temp.data[,  1] < prod.cutoff[1] &  temp.data[,  2] >= prod.cutoff[2]), 3]
  temp.data.q3 <- temp.data[which(temp.data[,  1] >= prod.cutoff[1] &  temp.data[,  2] >= prod.cutoff[2]), 3]
  temp.data.q4 <- temp.data[which(temp.data[,  1] >= prod.cutoff[1] &  temp.data[,  2] < prod.cutoff[2]), 3]

  ### cut all cells which are not producing cells of feature C
  # affects output "black" and "red"
  if (MSIplus) {
    temp.data.plus <- temp.data[which(temp.data[, 3] > prod.cutoff[3]), ]
    ncells <- nrow(temp.data.plus)
  }

  if (output == "black") {
    ### q[x].total [ink=black]
    ### percentage of cells in quadrant to total cells
    ### or in MSI(+): percentage of cells in quadrant to total positive cells
    q1 <- abs(100 * length(temp.data.q1) / ncells)
    q2 <- abs(100 * length(temp.data.q2) / ncells)
    q3 <- abs(100 * length(temp.data.q3) / ncells)
    q4 <- abs(100 - q1 - q2 - q3)
  } else if (output == "red") {
  ### number of cells which are producing cells in feature C
  # ncells <- nrow(temp.data[which(temp.data[, 3]> prod.cutoff[3]), ])

  ### q[x].prodcells [ink=red]
  ### percentage of cells which are positive for feature C in quadrant to total quadrant cells
  q1 <- 100 * length(temp.data.q1[which(temp.data.q1 >= prod.cutoff[3])]) / length(temp.data.q1)
  if (is.nan(q1)) q1 <- 0
  q2 <- 100 * length(temp.data.q2[which(temp.data.q2 >= prod.cutoff[3])]) / length(temp.data.q2)
  if (is.nan(q2)) q2 <- 0
  q3 <- 100 * length(temp.data.q3[which(temp.data.q3 >= prod.cutoff[3])]) / length(temp.data.q3)
  if (is.nan(q3)) q3 <- 0
  q4 <- 100 * length(temp.data.q4[which(temp.data.q4 >= prod.cutoff[3])]) / length(temp.data.q4)
  if (is.nan(q4)) q4 <- 0
  } else {
    ### percentage of cells which are positive for feature C to total cells
    q1 <- 100 * length(temp.data.q1[which(temp.data.q1 >= prod.cutoff[3])]) / ncells.total
    if (is.nan(q1)) q1 <- 0
    q2 <- 100 * length(temp.data.q2[which(temp.data.q2 >= prod.cutoff[3])]) / ncells.total
    if (is.nan(q2)) q2 <- 0
    q3 <- 100 * length(temp.data.q3[which(temp.data.q3 >= prod.cutoff[3])]) / ncells.total
    if (is.nan(q3)) q3 <- 0
    q4 <- 100 * length(temp.data.q4[which(temp.data.q4 >= prod.cutoff[3])]) / ncells.total
    if (is.nan(q4)) q4 <- 0
  }
  ### return
  return(c(round(q1, 3),
            round(q2, 3),
            round(q3, 3),
            round(q4, 3)
  ))

}



########################### Non Redundancy Score Calculation
NRS <- function(x, ncomp = 3){
  pr <- prcomp(x, center = TRUE, scale. = FALSE)
  score <- rowSums(outer(rep(1, ncol(x)),
                         pr$sdev[1:ncomp] ^ 2) * abs(pr$rotation[,  1:ncomp]))

  return(score)
}
########################### Non Redundancy Score Calculation END


########################### Stratification
##### source: https://gist.github.com/mrdwab/6424112#file-stratified-r
##### @param
# df        The input data.frame
# group     A character vector of the column or columns that make up the "strata".
# size      The desired sample size.
#           If size is a value less than 1, a proportionate sample is taken from each stratum.
#           If size is a single integer of 1 or more, that number of samples is taken from each stratum.
#           If size is a vector of integers, the specified number of samples is taken for each stratum.
#           It is recommended that you use a named vector. For example, if you have two strata, "A" and "B", and you wanted 5 samples from "A" and 10 from "B", you would enter size <- c(A <- 5, B <- 10).
# select    This allows you to subset the groups in the sampling process. This is a list. For instance, if your group variable was "Group", and it contained three strata, "A", "B", and "C", but you only wanted to sample from "A" and "C", you can use select <- list(Group <- c("A", "C")).
# replace   For sampling with replacement.
fcs$stratified <- function(df, group, size, select = NULL,
                           replace = FALSE, bothSets = FALSE) {
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop <- TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      # message(sQuote("size"), " vector entered as:\n\nsize <- structure(c(",
      #         paste(n, collapse <- ", "), "),\n.Names <- c(",
      #         paste(shQuote(names(n)), collapse <- ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse <- ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse <- ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits <- 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out <- length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse <- ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x <- size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace <- replace), ])
  set1 <- do.call("rbind", temp)

  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 <- set1, SET2 <- set2)
  } else {
    set1
  }
}

##### Select samples after a condition
##### @param
# df          The input data.frame
# cutoffs     cutoff vector
# condition   1.1 for (1.) %(CD34+CD38)>0.1 &  (2.) %(TdT)>0.1  & TdT, pSTAT5 and CD24 on X/Y-axes (later)
# frequencies optional, change frequencies in conditions (default=NA)
fcs$condition_approval <- function(
  df.sample,
  cutoffs,
  condition,
  freqs = NA
){

  approval <- FALSE
  if (condition == 1.1) {
    freqs #01 is for CD34/CD38 and 02 is for TdT
    if (any(is.na(freqs))) freqs <- c(0.1, 0.1)
    ### get marker indices, for cutoffs is one less (-1)
    ncells <- dim(df.sample)[1]
    sub.cells01 <- length(which(df.sample$CD34 > cutoffs["CD34"] & df.sample$CD38 > cutoffs["CD38"]))
    sub.cells02 <- length(which(df.sample$TdT > cutoffs["TdT"]))

    # (1.) %(CD34+CD38)>0.1
    freq01 <- sub.cells01 / ncells * 100
    freq02 <- sub.cells02 / ncells * 100
    # (2.) %(TdT)>0.1
    # (1.) %(CD34+CD38)>0.1 &  (2.) %(TdT)>0.1 &
    if (freq01 >= freqs[1] & freq02 >= freqs[2]) {
      approval <- TRUE

      # printf("%s for %s: freq01(%s)>a=%s freq02(%s)>b=%s",approval,df.sample[1,1],round(freq01, 3),freqs[1],round(freq02, 3),freqs[2])
    }
  }

  approval
}

