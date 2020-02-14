#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high ristk

rm(list = ls())



############### initiate
cofactor = 0.2
stat.info = "freq_green"
# stat.info = "absRange"
# coverage = "func"
coverage = "full"
cluster.size = 10

### set paths
sub.path = "/scratch/drfz/Good2018/glmnet/"
setwd(sub.path)
project.name = "Basal"
Training.table.name = sprintf("%s/Rdata/%s_Training_%s_quadrant_%s_cof%s.rds",sub.path,project.name,coverage,stat.info,cofactor)
Validation.table.name = sprintf("%s/Rdata/%s_Validation_%s_quadrant_%s_cof%s.rds",sub.path,project.name,coverage,stat.info,cofactor)
################

### load libraries
# xlsx            :   Library for Excel reading and creating
# dplyr           :   faster binding of columns bind_cols() and rows bind_rows()
# glmnet          :   generalized linear regression and cross validation
# doParallel      :   use several clusters for parallele calculations
libraries = c("xlsx","dplyr","glmnet","doParallel")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/Good2018/PRI_funs.R")

### metatable file
cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort = cohort_full[,c(1,5,8,11,16,15)]
train.set = cohort[which(cohort$Cohort=="Training"),1:5]
valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
# test.set = cohort[which(cohort$Cohort=="NA"),1:5]
trainval.set = bind_rows(train.set,valid.set)

### read table
df.Training = readRDS(Training.table.name)
df.Validation = readRDS(Validation.table.name)


df.total = bind_rows(df.Training,df.Validation)
sample.size = ncol(df.total)
rownames(df.total) = c(rownames(df.Training),rownames(df.Validation))
## column index of predicted variable in dataset
typeColNum = 1

condition = trainval.set$`Relapse Status`[which(trainval.set$`Patient ID` %in% rownames(df.total))]
### change condition "Yes2" / "Yes" to 1
condition[condition %in% c("Yes","Yes2")] = 1
### change condition "No" to 0
condition[condition %in% c("No")] = 0
condition = as.numeric(condition)


df.total= bind_cols(as.data.frame(condition),df.total)
names(df.total)[typeColNum] = "condition"
df.total = as.matrix(df.total)

### convert NaN/+-Inf to -1
if (any(is.nan(df.total))) df.total[is.nan(df.total) | is.infinite(df.total)] <- -1
### convert NAs to sample group mean
if (any(is.na(df.total))) {
  for ( i in 2:ncol(df.total)) {
    NA.idx = which(is.na(df.total[,i]))
    
    for (j in NA.idx) {
      tmp = df.total[which(df.total[j,1]==df.total[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.total[j,i] = tmp
    }
  }
}

alphalist <- seq(0,1,by=0.1)
alpha.collect = lambda1se.collect = vector()
seed.vec = sample(10^3)
it.total = 0

cl <- makeCluster(cluster.size)
registerDoParallel(cl)

timeSTART = Sys.time()
ptm <- proc.time()
printf("##### Start at %s.",timeSTART)
printf("run %s::%s::cluster=%s::samplesize=%s",project.name,stat.info,cluster.size,sample.size)
while (it.total < 100) {
  it.total = it.total + 1  
  if ( it.total %% 10 == 0 ) {
    printf("find_alpha %s run #%s..",stat.info,it.total)
    print(Sys.time()-timeSTART)
    print(proc.time() - ptm)
  }
  
  set.seed(seed.vec[it.total])
  ### 10-fold with 5-6 samples in one fold
  set.foldid = sample(rep(seq((1/5)*nrow(df.total)),length=nrow(df.total)))
  
  ### no need to set nfolds if foldid is provided since observations are already distributed in folds
  elasticnet <- lapply(alphalist, function(a) {
    cv.glmnet(x = df.total[,-typeColNum], y = df.total[,typeColNum],
              alpha=a, family="binomial",
              lambda.min.ratio=.0005,
              type.measure="deviance",
              foldid = set.foldid,
              parallel = TRUE
    )
  })
  
  min.err = min.err.lambda.1se = vector()
  for (i in 1:11) {
    #print(min(elasticnet[[i]]$cvm))
    #printf("min(CVmean)=%s lambda.1se=%s",min(elasticnet[[i]]$cvm),elasticnet[[i]]$lambda.1se)
    min.err = c(min.err,min(elasticnet[[i]]$cvm))
    min.err.lambda.1se = c(min.err.lambda.1se,elasticnet[[i]]$lambda.1se)
  }
  min.err.idx = which(min.err==min(min.err))
  alpha.best = alphalist[min.err.idx]
  lambda.best = elasticnet[[min.err.idx]]$lambda.1se
  # plot(min.err, ylab= "min(deviance)",
  #      main=sprintf("a=%s,seed=%s,lambda.1se=%s",alpha.best,seed.vec[it.total],lambda.best))
  
  lambda1se.collect = c(lambda1se.collect, lambda.best)
  alpha.collect = c(alpha.collect, alpha.best)
}
stopCluster(cl)
print(Sys.time())

alpha.tab = table(alpha.collect)
alpha.max = alpha.tab[which(alpha.tab==max(alpha.tab))]


printf("Best alpha after %s iterations: a=%s", it.total,names(alpha.max))
print(alpha.tab)


### create directory and file
sub.path = file.path(sub.path, sprintf("%s_findalpha_%s_%s",project.name,coverage,stat.info))
dir.create(sub.path, showWarnings = TRUE)
current.date = format(Sys.Date(),"%y%m%d")
### write information in log file
file.log<-file(sprintf("%s/%s_bestalpha.log",sub.path,current.date))
sink(file.log)
cat(current.date)
cat(sprintf("\nTime spent: %s\n",Sys.time()-timeSTART))
cat(sprintf("coverage = %s\n",coverage))
cat(sprintf("stat.info = %s\n",stat.info))
cat(sprintf("it.total = %s \n",it.total))
cat(sprintf("alpha_max = %s\n",names(alpha.max)))
cat(sprintf("%s\n",names(alpha.tab)))
cat(sprintf("%s\n",alpha.tab))
sink()

printf("Ready for %s::%s::%s::%s.",project.name,stat.info,coverage,it.total)
print(Sys.time()-timeSTART)
print(proc.time() - ptm)











