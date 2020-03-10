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
stat.info = "absRange"
subgroup = "Training"
coverage = "func"
# coverage = "full"
alphalist <- seq(0,1,by=0.002)
cluster.size = 3
condi = 1.1

### set paths
sub.path = "/scratch/drfz/Good2018/glmnet/"
setwd(sub.path)
project.name = "Basal"
Training.table.name = sprintf("%s/Rdata/%s_Training_%s_quadrant_%s_condi%s_cof%s.rds",sub.path,project.name,coverage,stat.info,condi,cofactor)
# save.Training.table.name = sprintf("%s/Rdata/%s_Training_%s_quadrant_%s_condi%s_cof%s.rds",sub.path,project.name,coverage,stat.info,condi,cofactor)
Validation.table.name = sprintf("%s/Rdata/%s_Validation_%s_quadrant_%s_cof%s.rds",sub.path,project.name,coverage,stat.info,cofactor)
sample.cond.name = sprintf("%s/%s_%s_%s_condi%s_cof%s.txt",sub.path,project.name,subgroup,coverage,condi,cofactor)
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

## -------------- metatable file, only look at samples in condition XX -------------
cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
# "Patient ID", "Prognostic Translocations", "MRD Risk", "Relapse Status", "DDPR Risk", "Cohort"  
cohort = cohort_full[,c(1,5,8,11,16,15)]
cond.samples = as.vector(unlist(read.table(sample.cond.name)))
train.set = cohort[which(cohort$`Patient ID` %in% cond.samples),1:5]
# valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
# test.set = cohort[which(cohort$Cohort=="NA"),1:5]
# trainval.set = bind_rows(train.set,valid.set)

## ---------- read table --------
df.Training = readRDS(Training.table.name)
df.idx = which(rownames(df.Training) %in% cond.samples)
df.Training = df.Training[df.idx,]
# df.Training = readRDS(Training.table.name)
# df.Validation = readRDS(Validation.table.name)
# df.total = bind_rows(df.Training,df.Validation)
# rownames(df.total) = c(rownames(df.Training),rownames(df.Validation))
sample.size = ncol(df.Training)
## column index of predicted variable in dataset
typeColNum = 1

condition = train.set$`Relapse Status`[which(train.set$`Patient ID` %in% rownames(df.Training))]
### change condition "Yes2" / "Yes" to 1
condition[condition %in% c("Yes","Yes2")] = 1
### change condition "No" to 0
condition[condition %in% c("No")] = 0
condition = as.numeric(condition)
df.Training= bind_cols(as.data.frame(condition),df.Training)
names(df.Training)[typeColNum] = "condition"
rownames(df.Training) = cond.samples

df.Training = as.matrix(df.Training)

## ----- exclude all variables with axes without TdT, pSTAT5 and CD24 -------------
# SKIP, already done at "0a_conditions_to_select...R"
if (FALSE) {
  marker.axes = c("TdT","pSTAT5","CD24")
  marker.comb = as.vector(outer(marker.axes, marker.axes, paste, sep="."))
  marker.comb = marker.comb[-c(1,5,9)]
  
  select.marker.idx = select.marker.names = vector()
  for ( i in 1:length(marker.comb) ) {
    temp.idx = eval(parse(text=paste0("grep(\"",marker.comb[i],"\",colnames(df.Training))")))
    temp.names = colnames(df.Training[,temp.idx])
    select.marker.names = c(select.marker.names, temp.names)
    select.marker.idx = c(select.marker.idx,temp.idx)
  }
  select.marker.df = t(as.data.frame(select.marker.idx))
  colnames(select.marker.df) = select.marker.names
  ## ------ identify doubled marker combinations and remove -------------------------
  colnames(select.marker.df) = make.unique(colnames(select.marker.df),sep="__")
  select.marker.df = select.marker.df[,-grep("__",colnames(select.marker.df))]
  # should be 424 +1 for condition
  df.Training = df.Training[,c(1,select.marker.df)]
  sample.size = ncol(df.Training)-1
}

# saveRDS(df.Training,save.Training.table.name)

### ------------ convert NaN/+-Inf/NA -----------------------------------------
nan.len = inf.len = na.len = 0
if (any(is.nan(df.Training))) {
  nan.len = length(is.nan(df.Training))
  df.Training[is.nan(df.Training)] <- -1
}
### infinite means there are no bins for absRange and no cells for freq_green
if (any(is.infinite(df.Training))) {
  inf.len = length(is.infinite(df.Training))
  df.Training[is.infinite(df.Training)] <- 0
}
### convert NAs to sample group mean
if (any(is.na(df.Training))) {
  na.len = length(is.na(df.Training))
  for ( i in 2:ncol(df.Training)) {
    NA.idx = which(is.na(df.Training[,i]))
    
    for (j in NA.idx) {
      tmp = df.Training[which(df.Training[j,1]==df.Training[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.Training[j,i] = tmp
    }
  }
}


alpha.collect = lambda1se.collect = vector()
seed.vec = sample(10^3)
it.total = 0

if (cluster.size>0) {
  cl <- makeCluster(cluster.size)
  registerDoParallel(cl)
}

timeSTART = Sys.time()
ptm <- proc.time()
printf("##### Start at %s.",timeSTART)
printf("run find_alpha on %s::%s::%s::cluster=%s::samplesize=%s",project.name,stat.info,subgroup,cluster.size,sample.size)
while (it.total < 100) {
  it.total = it.total + 1  
  if ( it.total %% 10 == 0 ) {
    printf("find_alpha %s run #%s..",stat.info,it.total)
    print(Sys.time()-timeSTART)
    print(proc.time() - ptm)
  }
  
  set.seed(seed.vec[it.total])
  ### 5-6 samples in one fold
  set.foldid = sample(rep(seq((1/5)*nrow(df.Training)),length=nrow(df.Training)))
  
  ### no need to set nfolds if foldid is provided since observations are already distributed in folds
  if (cluster.size==0) {
    elasticnet <- lapply(alphalist, function(a) {
      cv.glmnet(x = df.Training[,-typeColNum], y = df.Training[,typeColNum],
                alpha=a, family="binomial",
                lambda.min.ratio=.0005,
                type.measure="deviance",
                foldid = set.foldid
      )
    })
  } else {
    elasticnet <- lapply(alphalist, function(a) {
      cv.glmnet(x = df.Training[,-typeColNum], y = df.Training[,typeColNum],
                alpha=a, family="binomial",
                lambda.min.ratio=.0005,
                type.measure="deviance",
                foldid = set.foldid,
                parallel = TRUE
      )
    })
  }
  
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
  # plot(min.err, ylab= "min(deviance)",xlab = "var_count",
       # main=sprintf("a=%s,seed=%s,lambda.1se=%s",alpha.best,seed.vec[it.total],round(lambda.best,5)))
  
  lambda1se.collect = c(lambda1se.collect, lambda.best)
  alpha.collect = c(alpha.collect, alpha.best)
}
if (cluster.size>0) stopCluster(cl)

print(Sys.time())
printf("%s NaN positions converted to -1. %s Inf positions converted to 0. %s NA positions converted to group mean.",nan.len,inf.len,na.len)


alpha.tab = table(alpha.collect)
alpha.max = alpha.tab[which(alpha.tab==max(alpha.tab))]


printf("Best alpha after %s iterations: a=%s", it.total,names(alpha.max))
print(alpha.tab)


### create directory and file
sub.path = file.path(sub.path, sprintf("%s_findalpha_%s_%s_%s",project.name,coverage,stat.info,condi))
dir.create(sub.path, showWarnings = TRUE)
current.date = format(Sys.Date(),"%y%m%d")
### write information in log file
file.log<-file(sprintf("%s/%s_bestalpha.log",sub.path,current.date))
sink(file.log)
cat(current.date)
cat(sprintf("\nTime spent: %s\n",Sys.time()-timeSTART))
cat(sprintf("coverage = %s\n",coverage))
cat(sprintf("stat.info = %s\n",stat.info))
cat(sprintf("condition = %s\n",condi))
cat(sprintf("it.total = %s \n",it.total))
cat(sprintf("alpha_max = %s\n",names(alpha.max)))
cat(sprintf("%s\n",names(alpha.tab)))
cat(sprintf("%s\n",alpha.tab))
cat("\n")
sink()

printf("Ready for %s::%s::%s::%s.",project.name,stat.info,coverage,it.total)
printf("%s NaN positions converted to -1. %s Inf positions converted to 0. %s NA positions converted to group mean.",nan.len,inf.len,na.len)
print(Sys.time()-timeSTART)
print(proc.time() - ptm)











