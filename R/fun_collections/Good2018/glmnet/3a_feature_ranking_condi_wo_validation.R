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
set.alpha = 0.01
### only at alpha = 0 set amount of quads at 100
quad.cap = 100
stat.info = "freq_green"
# stat.info = "absRange"
coverage = "func"
# coverage = "full"
condi = 1.1
iterations = 500
rmse_threshold = 0.5
with_pVal = FALSE
cluster.size = 3

### set paths
folder.path = "/scratch/drfz/Good2018/glmnet"
setwd(folder.path)
project.name = "Basal"
### quadrant file names
Training.table.name = sprintf("%s/Rdata/%s_Training_%s_quadrant_%s_condi%s_cof%s.rds",folder.path,project.name,coverage,stat.info,condi,cofactor)
Validation.table.name = sprintf("%s/Rdata/%s_Validation_%s_quadrant_%s_condi%s_cof%s.rds",folder.path,project.name,coverage,stat.info,condi,cofactor)
Training.cond.name = sprintf("%s/%s_Training_%s_condi%s_cof%s.txt",folder.path,project.name,coverage,condi,cofactor)
Validation.cond.name = sprintf("%s/%s_Validation_%s_condi%s_cof%s.txt",folder.path,project.name,coverage,condi,cofactor)
################

### load libraries
# xlsx            :   Library for Excel reading and creating
# dplyr           :   faster binding of columns bind_cols() and rows bind_rows()
# glmnet          :   generalized linear regression and cross validation
# ggpubr          :   https://www.r-bloggers.com/add-p-values-and-significance-levels-to-ggplots/
# reshape2        :   function melt()
# doParallel      :   use several clusters for parallele calculations
libraries = c("xlsx","dplyr","glmnet","ggpubr","reshape2","doParallel")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/Good2018/PRI_funs.R")

### metatable file
cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort = cohort_full[,c(1,5,8,11,16,15)]
train.cond.samples = as.vector(unlist(read.table(Training.cond.name)))
train.set = cohort[which(cohort$`Patient ID` %in% train.cond.samples),1:5]

valid.cond.samples = as.vector(unlist(read.table(Validation.cond.name)))
valid.set = cohort[which(cohort$`Patient ID` %in% valid.cond.samples),1:5]

# valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
# test.set = cohort[which(cohort$Cohort=="NA"),1:5]

# ## column index of predicted variable in dataset
typeColNum = 1

### read table
df.Training = readRDS(Training.table.name)
sample.size = ncol(df.Training)
### add condition
condition = train.set$`Relapse Status`[which(train.set$`Patient ID` %in% rownames(df.Training))]
### change condition "Yes2" / "Yes" to 1
condition[condition %in% c("Yes","Yes2")] = 1
### change condition "No" to 0
condition[condition %in% c("No")] = 0
condition = as.numeric(condition)
df.Training= bind_cols(as.data.frame(condition),df.Training)
names(df.Training)[typeColNum] = "condition"
rownames(df.Training) = train.cond.samples
df.Training = as.matrix(df.Training)

df.Validation = readRDS(Validation.table.name)
tmpnames.Val = rownames(df.Validation)

# ### condition 1/0
# condition.train = train.set$`Relapse Status`[which(train.set$`Patient ID` %in% rownames(df.Training))]
condition.val = valid.set$`Relapse Status`[which(valid.set$`Patient ID` %in% rownames(df.Validation))]
# ### change condition "Yes2" / "Yes" to 1
# condition.train[condition.train %in% c("Yes","Yes2")] = 1
condition.val[condition.val %in% c("Yes","Yes2")] = 1
# ### change condition "No" to 0
# condition.train[condition.train %in% c("No")] = 0
condition.val[condition.val %in% c("No")] = 0
# condition.train = as.numeric(condition.train)
condition.val = as.numeric(condition.val)
# 
# # df.Training = bind_cols(as.data.frame(condition.train),df.Training)
df.Validation = bind_cols(as.data.frame(condition.val),df.Validation)
# names(df.Training)[typeColNum] = "condition"
names(df.Validation)[typeColNum] = "condition"
# df.Training = as.matrix(df.Training)
df.Validation = as.matrix(df.Validation)
rownames(df.Validation) = tmpnames.Val


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

### ------------ convert NaN/+-Inf/NA df.Validation-----------------------------------------
if (any(is.nan(df.Validation))) df.Validation[is.nan(df.Validation)] <- -1
### infinite means there are no bins for absRange and no cells for freq_green
if (any(is.infinite(df.Validation))) df.Validation[is.infinite(df.Validation)] <- 0
### convert NAs to sample group mean
if (any(is.na(df.Validation))) {
  for ( i in 2:ncol(df.Validation)) {
    NA.idx = which(is.na(df.Validation[,i]))
    
    for (j in NA.idx) {
      tmp = df.Validation[which(df.Validation[j,1]==df.Validation[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.Validation[j,i] = tmp
    }
  }
}

seed.vec = sample(10^3)
it_rmse = 0
it_total = 0
conf_interval=0.95

acc_total = vector()
padj_total = 0


cl <- makeCluster(cluster.size)
registerDoParallel(cl)

current.date = format(Sys.Date(),"%y%m%d")

timeSTART = Sys.time()
ptm <- proc.time()
printf("##### Start at %s.",timeSTART)
printf("run %s::%s::%s::condi=%s::alpha=%s::samplesize=%s",project.name,stat.info,coverage,condi,set.alpha,sample.size)
while (it_total < iterations) {
  
  it_total = it_total +1
  
  set.seed(seed.vec[it_total])
  ######## OLD SAMPLING
  #set.foldid = sample(rep(seq((1/10)*nrow(df.Training)),length=nrow(df.Training)))
  ######## NEW
  ### FIND OUT THE FREQUENCY OF condition "1" and "0" once
  # d = as.data.frame(table(condition.train))
  # d$Freq[2]/d$Freq[1]&
  # 0.4193548
  ##### IMPORTANT equal stratification ten times
  df.tmp = as.data.frame(df.Training[,1:2])
  strat.df = rep(0,nrow(df.Training))
  strat.df = as.data.frame(strat.df)
  rownames(strat.df) = rownames(df.Training)
  for ( i in 1:5 ) {
    ### pick 3 OR 4 from condition = 0
    ### pick 1 OR 2 from condition = 1
    # df.strat = fcs$stratified(df.tmp,"condition",size=c(sample(c(3,4),1),sample(c(1,2),1)),replace=FALSE)
    
    ### pick 3 from condition = 0
    ### pick 1 from condition = 1
    df.strat = fcs$stratified(df.tmp,"condition",size=c(3,2),replace=FALSE)
    strat.df[which(rownames(strat.df) %in% rownames(df.strat)),] = i
    
    strat.idx = which(rownames(df.tmp) %in% rownames(df.strat))
    df.tmp = df.tmp[-strat.idx,]
  }
  # check with table(strat.df)
  strat.df[which(strat.df==0),] = 6
  strat.df = as.vector(unlist(strat.df))
  #######
  
  ### no need to set nfolds if foldid is provided since observations are already distributed in stratified folds
  if (cluster.size==0) {
    cv.out <- cv.glmnet(x = as.matrix(df.Training[,-typeColNum]), y = df.Training[,typeColNum],
                        alpha=set.alpha, family="binomial",
                        lambda.min.ratio=.0005,
                        type.measure="deviance",
                        foldid = strat.df
    )
  } else {
    cv.out <- cv.glmnet(x = as.matrix(df.Training[,-typeColNum]), y = df.Training[,typeColNum],
                        alpha=set.alpha, family="binomial",
                        lambda.min.ratio=.0005,
                        type.measure="deviance",
                        foldid = strat.df,
                        parallel = TRUE
    )
  }
  # plot result
  # par(oma=c(1,1,1,1),mar=c(1,1,1,1))
  plot(cv.out)
  
  # min value of lambda
  lambda_min <- cv.out$lambda.min
  # best value of lambda
  lambda_1se <- cv.out$lambda.1se
  
  ### prediction 
  elastic_class <- predict(cv.out,
                         newx = df.Validation[,-typeColNum],
                         s=lambda_min,
                         na.action = na.pass,
                         type="class"
  )
  elastic_prob <- predict(cv.out,newx = as.matrix(df.Validation[,-typeColNum]),s=lambda_1se,type="response")
  elastic_train <- predict(cv.out,newx = as.matrix(df.Training[,-typeColNum]),s=lambda_1se,type="response")
  
  ### accuracy
  acc = mean(elastic_class==df.Validation[,typeColNum])
  acc_total = c(acc_total,acc)
  
  ### root mean squared error RMSE
  # in sample RMSE
  rmse_train = sqrt(sum(elastic_train - mean(df.Training[,typeColNum]))^2/nrow(df.Training))
  # out of sample RMSE
  rmse_test = sqrt(sum(elastic_prob - mean(df.Validation[,typeColNum]))^2/nrow(df.Validation))
  
  # if (rmse_test>=0.05) printf("it=%s RMSE(train)=%s RMSE(test)=%s",it_total,rmse_train,rmse_test);it_rmse=it_rmse+1
  # printf("%s",rmse_test-rmse_train)
  rmse_diff = rmse_test-rmse_train
  
  
  #regression coefficients
  coeff = coef(cv.out,s=lambda_1se)
  # coeff = coef(cv.out,s=lambda_min)
  ### get coefficients who has impact
  coeff.idx = which(abs(coeff)>0)
  
  coeff.pred = coeff[coeff.idx,]
  coeff.pred = coeff.pred[order(abs(coeff.pred),decreasing=T)]
  
  ### only if ridge regression only
  if (set.alpha == 0) {
    coeff.pred[1:quad.cap]
  }
  
  if ( it_total %% 50 == 0 ) {
    printf("%s::%s::%s::a=%s::cluster=%s::acc=%s::rmse_diff=%s",it_total,stat.info,coverage,set.alpha,cluster.size,round(acc,3),round(rmse_diff,3))
    print(Sys.time()-timeSTART)
    print(proc.time() - ptm)
  }
  
  if ((rmse_diff)<rmse_threshold & length(coeff.pred)>1) {
    it_rmse = it_rmse+1
    
    ### look at differentiating variables (here = quadrants)
    intercept.idx = grep("Intercept",names(coeff.pred))
    vars.pred = names(coeff.pred)[-intercept.idx]
    vars.pred = gsub("`","",vars.pred)
    
    df.Training = as.data.frame(df.Training)
    df.predict = df.Training[order(rownames(df.Training)),c(1,which(names(df.Training)%in%vars.pred))]
    
    
    ### manipulate dataframe
    dat.m = melt(df.predict, id.vars="condition")
    
    ### unpaired two samples t-test (independent samples)
    #t.test(df.predict[which(df.predict$condition==1),9],df.predict[which(df.predict$condition==0),9])
    ### calculate t.test for each variable
    comp = compare_means(value ~ condition, data = dat.m, group.by = "variable",  
                         method="t.test", p.adjust.method = "BH"
    )
    
    if (with_pVal) {
      ### get idx where adjusted p-value is <0.05
      comp.adj.idx = which(comp$p.adj<0.05)
    } else {
      comp.adj.idx = which(comp$p.adj<1)
    }
    
    padj_total = padj_total+length(comp.adj.idx)
    
    ### get quadrants and coeffs with p.adj<0.05
    coeff.pred.adj = coeff.pred[which(names(coeff.pred) %in% comp$variable[comp.adj.idx] )]
    
    
    ### create directory and file
    sub.path = file.path(folder.path, sprintf("%s_FR_%s_%s_a%s_RMSE%s",project.name,coverage,stat.info,set.alpha,rmse_threshold))
    dir.create(sub.path, showWarnings = TRUE)
    signif.stat.file = sprintf("%s/%s_signif_stat_it%s.txt",sub.path,current.date,iterations)
    signif.vars.file =  sprintf("%s/%s_vars_of_signif_stat_it%s.txt",sub.path,current.date,iterations)
    
    ### write info to signif_stat.txt
    # create file and header
    if (it_rmse==1) write(paste("project.name","calculation_meth","sample_size","acc","prediction_RMSE(Val)","seed","padj",sep="\t"),
     file=sprintf("%s/%s_signif_stat_it%s.txt",sub.path,current.date,iterations))
    write(paste(project.name,stat.info,"stat",round(acc,4),round(rmse_test,4),seed.vec[it_total],
                paste(names(coeff.pred.adj),collapse="\t"),
                sep="\t"),file=signif.stat.file,append=TRUE)
    write(paste(project.name,stat.info,"coeff",round(acc,4),round(rmse_test,4),seed.vec[it_total],
                paste(coeff.pred.adj,collapse="\t"),
                sep="\t"),file=signif.stat.file,append=TRUE)
    write(paste(project.name,stat.info,"pVal",round(acc,4),round(rmse_test,4),seed.vec[it_total],
                paste(round(comp$p[comp.adj.idx],4),collapse="\t"),
                sep="\t"),file=signif.stat.file,append=TRUE)
    write(paste(project.name,stat.info,"pAdj",round(acc,4),round(rmse_test,4),seed.vec[it_total],
                paste(round(comp$p.adj[comp.adj.idx],4),collapse="\t"),
                sep="\t"),file=signif.stat.file,append=TRUE)
    
    
    ### get protein names    
    if (coverage=="full") {
      columns_name = "sorted"
    } else columns_name = "functional"
    col.vec.func = as.vector(unlist(read.table(file=sprintf("/scratch/drfz/Good2018/columns_%s.txt", columns_name))))
    
    ### count protein in protein combination of quadrants and calculation method
    var.list = sapply(col.vec.func, grep, comp$variable[comp.adj.idx])#,value=T)
    # var.list[which( lengths(var.list) != 0)]
    # lengths(var.list[which( lengths(var.list) != 0)])
    # unlist(sapply(col.vec.func, grepl, comp$variable[comp.adj.idx]))
    var.count = lengths(var.list)
    
    ### write only counts of proteins from significant quadrants
    # create file and header
    if (it_rmse==1) write(paste("project.name","calculation_meth","sample_size","acc","prediction_RMSE(Val)","seed","padj",
               paste(col.vec.func,collapse="\t"),sep="\t"),
               file=signif.vars.file)
    write(paste(project.name,stat.info,sample.size,round(acc,4),round(rmse_test,2),seed.vec[it_total],with_pVal,paste(var.count,collapse="\t"),sep="\t"),
          file=signif.vars.file,append=TRUE)
  }
}
stopCluster(cl)

printf("Ready for %s::%s::%s::%s::pVal=%s::cluster=%s::%s/%s with RMSE_diff<%s.",project.name,stat.info,coverage,set.alpha,with_pVal,cluster.size,it_rmse,it_total,rmse_threshold)
printf("%s NaN positions converted to -1. %s Inf positions converted to 0. %s NA positions converted to group mean.",nan.len,inf.len,na.len)
print(Sys.time()-timeSTART)
print(proc.time() - ptm)



# rm(list = ls())








