#!/usr/bin/R
# DRFZ 2019 - 2020
# Goood2018

add.libraries <- c("glmnet", "doParallel")
lapply(add.libraries, require, character.only = TRUE)

print("### Running cox hazard model with OOB on df.total..")
### ------------------ cox hazard model in cross validation
### initiate variables for cross validation / cox model


cl <- makeCluster(cluster.size)
registerDoParallel(cl)

timeStart <- Sys.time()
ptm <- proc.time()
printf("### Start %s with %s on dataset=%s and subset=%s lambda=%s random=%s spikes=%s ###", timeStart, model[model.id], dataset[dataset.id], subset[subset.id], lambda[lambda.id], randomize.label, spikeIns)

### run OOB (out of bag) through all patients one by one from df.total
glmnet.oob <- lapply( seq(1:2), function(a) {
# gbm.oob <- lapply( seq(1:nrow(df.total)), function(a) {
  df.train <- df.total[-a, ]
  df.val <- df.total[a, ]

  it.total <- 0
  while (it.total < sampling.size) {
    it.total <- it.total + 1 

    set.seed(seed.vec[it.total])

    ### --------------------- create fold values
    # Basal has 13 relapsed and 31 relapse free patients
    # set folds for cross validation manual because of imbalance data
    # set folds that 1 fold contains at least 3 relapsed patients
    fold.id <- rep(0, nrow(df.train))
    fold.n <- 4 # 4 folds
    zero.num <- length(which(df.train[, 2] == 0))
    one.num <- length(which(df.train[, 2] == 1))
    if ((zero.num %% fold.n) > 0) {
      zero.id <- rep(sample(1:fold.n), zero.num / fold.n)
      zero.id.mod <- sample(1 : (zero.num %% fold.n))
      zero.id <- c(zero.id, zero.id.mod)
    } else {
      zero.id <- rep(sample(1:fold.n), zero.num / fold.n)
    }
    if ((one.num %% fold.n) > 0) {
      one.id <- rep(sample(1:fold.n), one.num / fold.n)
      one.id.mod <- sample(1: (one.num %% fold.n))
      one.id <- c(one.id, one.id.mod)
    } else {
      one.id <- rep(sample(1:fold.n), one.num / fold.n)
    }
    it.set <- it.null <- it.one <- 1

    #create 4 folds with at least 3 relapsed patients
    while (it.set < nrow(df.train) + 1){
        if (df.train[it.set, 2] == 0) {
        # if sample relapse status is false / 0
          fold.id[it.set] <- zero.id[it.null]
          it.null <- it.null + 1
      } else {
        fold.id[it.set] <- one.id[it.one]
        it.one <- it.one  + 1
      }
      it.set <- it.set  + 1
    }
    ##############
    
    ### ----------- create model with cross validation
    cv.fit <- cv.glmnet(as.matrix(df.train[, -c(1:3)]),
                        Surv(df.train[, 1], df.train[, 2]), 
                        family = "cox",
                        alpha = set.alpha, 
                        foldid = fold.id, 
                        parallel = TRUE
                        )
  }
})
stopCluster(cl)
