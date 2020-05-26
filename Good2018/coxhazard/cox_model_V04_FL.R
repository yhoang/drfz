#!/usr/bin/R
# Authors: Eric Urbansky, Yen Hoang and Felix Lohrke



# naming order: dataset_subgroup_subset_subject_featureinfo_comment_initalsdate
rm(list = ls())
#install.packages(c("survival","glmnet","readxl","dplyr","doParallel","xlsx","survminer","pROC","survAUC"))
#necessary r packages
libarys = c("survival","glmnet","readxl","dplyr","doParallel","xlsx","survminer","pROC","survAUC","stringr")
lapply(libarys,require,character.only = TRUE)

#a useful print function
printf <- function(...) invisible(print(sprintf(...)))

# command line parser
args = commandArgs(trailingOnly = TRUE)

############ CONFIGURATIONS #################


### data to work with
dataset.selected = 1
subset.selected = 1

# comment needs to be written as _comment
comment = ""

#  mostly static
subject.selected = 1
featureinfo.selected = 1
cof.selected = 1

### model parameters
alpha.value = 1
SENS.value = 1
FP.value = 0.2

#threshold to select (should multiple apply to parameters)
selected.thresh = 1

### remove 4 patients?
patients.removed = F

### output config

# Output activated?
# generate output that doesnt change in the model (collection of ROC values, ROC curve, variables, heat map)
# generally only needed to apply once per model
output.model = T

# generate output for different selected thresholds (threshold selected + info, AUC)
output.thresholds = F

############ PARAMETERS FILENAMING #################

# dataset (Basal, BCR, IL7, Pervanadate, TSLP)
dataset = c("Basal", "BCR", "IL7", "Pervanadate", "TSLP")

# subgroup (training, validation, total)
subgroup = c("Training", "Validation", "total")

# subset (full, func, func3, func6)
subset = c("full", "func", "func3", "func6")

# subject (quadrant, df, HIST, ROC, heatmap)
subject = c("quadrant", "df","ROC", "AUC","iAUC","CompV",
            "CompT","CompTotal","Variables","CompAll","KM","Thresholds")

# feature information used (absRange, mean, var, zRange, freq.green)
featureinfo = c("absRange", "mean", "var", "zRange", "freq.green")

# coef
cof = "cof0.2"

# initial & date
working.station = "FL"
initdate = paste0(working.station,substring(str_replace_all(Sys.Date(),"-",""),3))

##### command line arguments #####

if (length(args) == 2) {
  dataset[dataset.selected] = args[1]
  subset[subset.selected] = args[2]
  print(dataset[dataset.selected])
  print(subset[subset.selected]) 
}

############ PARAMETERS PATHING #################

universal_file_naming = paste(dataset[dataset.selected], "_", 
                            subset[subject.selected], "_", 
                            subject[subject.selected], "_", 
                            featureinfo[featureinfo.selected], "_", 
                            cof[cof.selected], 
                            comment, "_",
                            initdate, sep="")
                            

variant_universal_filenaming <- function(subject.variation=1, 
                                d = dataset[dataset.selected],
                                sub = subset[subset.selected],
                                feature = featureinfo[featureinfo.selected],
                                c = cof[cof.selected], 
                                com = comment,
                                init = initdate) {

        subject.new = c("quadrant", "df","ROC", "AUC","iAUC","CompV",
            "CompT","CompTotal","Variables","Heat","CompAll","KM","Thresholds")

        universal_file_naming = paste(d, "_", 
        sub, "_", 
        subject.new[subject.variation], "_", 
        feature, "_", 
        c, 
        com, "_",
        init, sep="")

        return(universal_file_naming)
}

results_path = "./Results/"

############ INPUT PATHS #################

### Path to Patient data
patient_data_path <- "./Data/patient_cohort.xlsx"

### Path to Subgroups
subgroup_selection = paste(subset[subset.selected],subject[subject.selected],featureinfo[featureinfo.selected],cof[cof.selected], sep="_")
training_data_path <- paste("./Rds/", dataset[dataset.selected],"/", dataset[dataset.selected],"_",subgroup[1],"_",subgroup_selection,".rds", sep="")
validation_data_path <- paste("./Rds/", dataset[dataset.selected],"/", dataset[dataset.selected],"_",subgroup[2],"_",subgroup_selection,".rds", sep="")
total_data_path <- paste("./Rds/",dataset[dataset.selected],"_",subgroup[3],"_",subgroup_selection, comment, ".rds", sep="")

### Path to model storage (if model doesnt exist a new one will be created)
model_rds_path = paste0("./Rds/Cox_models/Coxmodel_", alpha.value, "_", paste(dataset[dataset.selected], "_", 
                            subset[subject.selected], "_", 
                            subject[subject.selected], "_", 
                            featureinfo[featureinfo.selected], "_", 
                            cof[cof.selected], 
                            comment,
                            sep="") , ".rds")

###########################################

############ OUTPUT PATHS #################
# roc, auc and iauc curves
roc_pdf_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(3),".pdf", sep="")
auc_pdf_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(4),".pdf", sep="")
iauc_pdf_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(5),".pdf", sep="")

# comparisons of validation, training and total xlsx
comp_validation_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(6), ".xlsx", sep="")
comp_training_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(7), ".xlsx", sep="")
comp_total_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(8), ".xlsx", sep="")

# variables xlsx
variables_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(9), ".xlsx", sep="")
# heatmap pdf
heat_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(10), ".pdf", sep="")

# comparisons of SENS,FP-rate and associated thresholds
comp_all_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(11), ".xlsx", sep="")

# kaplan meier plot
kaplan_pdf_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(12), ".pdf", sep="")
#kaplan.ddpr_pdf_path <- paste(results_path, specifier,"/Kaplan-Plot_DDPR_", Sys.Date(),"_",specifier,".pdf", sep="")

# thresholds txt
thresholds_path <- paste(results_path, dataset[dataset.selected],"/",variant_universal_filenaming(13), ".txt", sep="")

###########################################

#loading patient data from excel
patient_data <- read_excel(patient_data_path)
#View(patient_data)

#printfunktion
#a usefull print function
printf <- function(...) invisible(print(sprintf(...)))

#set columstype from patien data str to numeric values
#colums with numeric value : Age at Diagnosis, WBC Count, Date of Diagnosis, Time to Relapse(Day), CCR (Day) 

patient_data$`Time to Relapse (Days)`<- as.numeric(patient_data$`Time to Relapse (Days)`)
patient_data$`CCR (Days)`<-as.numeric(patient_data$`CCR (Days)`)
patient_data$`Time to Relapse (Days)`[is.na(patient_data$`Time to Relapse (Days)`)]<-0
patient_data$`CCR (Days)`[is.na(patient_data$`CCR (Days)`)]<-0
patient_data$`Age at Diagnosis`<-as.numeric(patient_data$`Age at Diagnosis`)
patient_data$`WBC Count`<-as.numeric(patient_data$`WBC Count`)
#patien_data$`Date of Diagnosis`<-as.numeric(patien_data$`Date of Diagnosis`)

#new Collum "Survival Time (Day)
patient_data$`Survival Time (Day)` <- patient_data$`Time to Relapse (Days)`+patient_data$`CCR (Days)`
#fehlermeldung unknown column <- "Survival Time (Day) vorher erstellen?

#yes2 == yes
patient_data$`Relapse Status`[patient_data$`Relapse Status`=="Yes2"]<-"Yes"

#Set DDPR Status to nuermic binary values
patient_data$`DDPR Risk`[patient_data$`DDPR Risk`=="Low"]<-0
patient_data$`DDPR Risk`[patient_data$`DDPR Risk`=="High"]<-1

#reduce patien data to necessary colums
cohort=patient_data[,c(1,11,8,16,15,17)]

#divide cohort into training and validation set
training.set=cohort[which(cohort$Cohort=="Training"),1:6]
validation.set=cohort[which(cohort$Cohort=="Validation"),1:6]
#bind training and validation set to totalset
total.set=bind_rows(training.set,validation.set)



#safe triplot rds file as df. for training and validation set
#absRange without condition
printf("Loading RDS")
df.training <- readRDS(training_data_path)
df.validation <- readRDS(validation_data_path)

# remove 4 patients that cause no modelvar bug

#newpatients = readRDS("newpatients.rds")
#patients.idx = which(row.names(df.training) == newpatients )

if (patients.removed) {
   newpatients = readRDS("newpatients.rds")
   df.training <- df.training[!row.names(df.training) %in% newpatients,]
}


printf("RDS loaded")
printf("Dimensions of training and validation set: ")
print(dim(df.training))
print(dim(df.validation))
print(rownames(df.training))
print(rownames(df.validation))


#absRange with Condi 1.1 
#df.training<-readRDS("~/Good/RDS/TriPlotData/Basal_Training_func_quadrant_absRange_condi1.1_cof0.2.rds")
#df.validation<-readRDS("~/Good/RDS/TriPlotData/Basal_Validation_func_quadrant_absRange_condi1.1_cof0.2.rds")

#safe rownames ("patien ID")

rownames.validation = row.names(df.validation)
rownames.training = row.names(df.training)

#### Check if df.total already exists as RDS (binding is bottleneck)
#### RDS with df.total will be created once if it doesnt exist and then read every run
if (file.exists(total_data_path)){
  df.total <- readRDS(total_data_path)
  printf("df.total does exist...using created RDS")
} else {
  printf("df.total doesnt exist...creating RDS")
  df.total<-bind_rows(df.training,df.validation)
  saveRDS(df.total, file = total_data_path)

}

print("Dimensions of df.total:")
print(dim(df.total))
rownames(df.total)=c(rownames.training,rownames.validation)

printf("Num of patients:")
print(nrow(df.total))
#View(df.total)
#total 60 patien in cohort
#no condition total 54  patient
#condition 1.1 total 41 patien

typeColNum=1

#missing patien from cohort becaus of condition 1.1
missin.patien = c(row.names(df.total),total.set$`Patient ID`)
missin.patien = names(which(table(missin.patien) == 1))

#relapsstatus in training set as numeric
condition = training.set$`Relapse Status`[which(training.set$`Patient ID`%in% row.names(df.training))]
condition[condition %in% c("Yes")]=1
condition[condition %in% c("No")]=0
condition=as.numeric(condition)

#survival time in training.set as numeric
survival=training.set$`Survival Time (Day)`[which(training.set$`Patient ID`%in%row.names(df.training))]
survival=as.numeric(survival)

#DDPR STATUS as numeric
DDPR=training.set$`DDPR Risk`[which(training.set$`Patient ID`%in%row.names(df.training))]
DDPR=as.numeric(DDPR)

#add condition/survivaltime/age to df.training
df.training=bind_cols(as.data.frame(DDPR),df.training)

names(df.training)[typeColNum] = "DDPR Status"
df.training=bind_cols(as.data.frame(condition),df.training)
names(df.training)[typeColNum]="Relaps Status"
df.training=bind_cols(as.data.frame(survival),df.training)
names(df.training)[typeColNum] = "Survivaltime (Day)"

#do same for validation set
condition = validation.set$`Relapse Status`[which(validation.set$`Patient ID`%in% row.names(df.validation))]
condition[condition %in% c("Yes")]=1
condition[condition %in% c("No")]=0
condition=as.numeric(condition)
survival = validation.set$`Survival Time (Day)`[which(validation.set$`Patient ID`%in%row.names(df.validation))]
survival = as.numeric(survival)
DDPR=validation.set$`DDPR Risk`[which(validation.set$`Patient ID`%in%row.names(df.validation))]
DDPR=as.numeric(DDPR)
df.validation=bind_cols(as.data.frame(DDPR),df.validation)
names(df.validation)[typeColNum] = "DDPR Status"
df.validation=bind_cols(as.data.frame(condition),df.validation)
names(df.validation)[typeColNum]="Relaps Status"
df.validation=bind_cols(as.data.frame(survival),df.validation)
names(df.validation)[typeColNum] = "Survivaltime (Day)"

df.training.strat = df.training
rownames(df.training.strat) = rownames.training

# glmnet need input matrix for model 
# safe df as matrix
###
#names.validation <- names(df.validation)
#names.training <- names(df.training)
###

df.training=as.matrix(df.training)
df.validation=as.matrix(df.validation)

### convert NaN/+-Inf to 0 in df.training
if (any(is.nan(df.training))) df.training[is.nan(df.training) | is.infinite(df.training)] <- -0.01
if (any(is.infinite(df.training))) df.training[is.infinite(df.training)] <- -0.01
### convert NAs to sample group mean
if (any(is.na(df.training))) {
  for ( i in 2:ncol(df.training)) {
    NA.idx = which(is.na(df.training[,i]))
    
    for (j in NA.idx) {
      tmp = df.training[which(df.training[j,1]==df.training[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.training[j,i] = tmp
    }
  }
}
### convert NaN/+-Inf to 0 in df.validation

if (any(is.nan(df.validation))) df.validation[is.nan(df.validation) | is.infinite(df.validation)] <- -0.01
if (any(is.infinite(df.validation))) df.validation[is.infinite(df.validation)] <- -0.01
### convert NAs to sample group mean
if (any(is.na(df.validation))) {
  for ( i in 2:ncol(df.validation)) {
    NA.idx = which(is.na(df.validation[,i]))
    
    for (j in NA.idx) {
      tmp = df.validation[which(df.validation[j,1]==df.validation[,1]),i]
      tmp = tmp[-which(is.na(tmp))]
      tmp = round(sample(seq(mean(tmp)-sd(tmp),mean(tmp)+sd(tmp),by=0.01),1),2)
      
      df.validation[j,i] = tmp
    }
  }
}

#create survival objekt(survival: Surv()) for cox model
sur_obj_validation = Surv(df.validation[,1],df.validation[,2])
# creating alternative surv object with 4 patients
#df.validation.patiets = df.training[patients.idx,]
#sur_obj_validation_patients = Surv(df.validation.patients[,1],df.validation.patients[,2])


#set seed for reproduction 
seed.vec = sample(10^2)
it.total = 0
cluster.size =3
all.activ.Index = all.coef.value = vector()
variabl.count = vector()
min.cvm = 100

###### START OF MODEL GENERATION #######
#### only runs if model does not already exist at model_rds_path

if (!file.exists(model_rds_path)) {
printf("Cox Model does not exist.....generating new Model")
   
cl <- makeCluster(cluster.size)

registerDoParallel(cl)

timeStart = Sys.time()
ptm <- proc.time()
printf("###Start %s.###",timeStart)

while (it.total < 100) {
  it.total = it.total + 1 


set.seed(seed.vec[it.total])

#set folds for cross validation manual because of imbalance data
#set folds that 1 fold contains at least 1 relaps
# result: 4 folds with 8 patients (at least 2 relaps)


#fold.id for condition 1.1
#create fold values
#null.id = one.id= vector()
#fold.id = rep(0,nrow(df.training))
#one.id = rep(sample(1:4),3)
#length(null.id)
#null.id = rep(sample(1:4),5)
#length(one.id)
#it.set =  it.null = it.one = 1

#fold.id for no condidion
#create fold values
null.id = one.id= vector()
fold.id = rep(0,nrow(df.training))
one.id = c(rep(sample(1:4),3),3)
length(null.id)
null.id = rep(sample(1:4),8)
length(one.id)
it.set =  it.null = it.one = 1

#create folds

while (it.set < nrow(df.training)+1){
  if (df.training[it.set,2] == 0){
      fold.id[it.set] = null.id[it.null]
      it.null = it.null+1}
  else{
    fold.id[it.set] = one.id[it.one]
    it.one = it.one +1
  }
  it.set = it.set +1
}

tmp <- df.training[,-c(1:3)]
#tmp <- tmp[,sample(length(tmp))]
#create Model with cross validation
cv.fit <- cv.glmnet(tmp,Surv(df.training[,1],df.training[,2]),
                    family = "cox",
                    alpha=alpha.value,
                    foldid = fold.id,
                    parallel= TRUE)
#plot(cv.fit)


#collect all active (!=0) coef from fited model
# lambda.min = model with min cross validation error

Coefficients<-coef(cv.fit,s=cv.fit$lambda.min)
Active.Index<-which(Coefficients!=0)
coef.value <- coef(cv.fit, s=cv.fit$lambda.min)
all.coef.value <- c(all.coef.value,coef.value)
variabl.count = c(variabl.count,length(Active.Index))
print(length(Active.Index))
tmp = min(cv.fit$cvm)

#collect best model (model with min cross validation error) 
#safe best model as op.fit
#safe active coef 

if(min.cvm > tmp){
  min.cvm = tmp
  op.fit = cv.fit
  op.variabl.count = length(Active.Index)
  op.index = Active.Index

  }
if ( it.total %% 10 == 0 ) {
  printf("At Work IT:%s ",it.total)
  print(Sys.time()-timeStart)
  print(proc.time()-ptm)
}

}
printf("End Process")
print(Sys.time()-timeStart)

stopCluster(cl)


print(table(variabl.count))
printf("Bestes Model mit:%s", op.variabl.count)
printf("Cross validation error: %s", min.cvm)


###### END OF MODEL GENERATION #######

#prediction for training and validation set based on op.fit
#prediction for type cox model = relativ risk (RR)


####### saving output to txt to check on df.training and df.validation
#sink("df.training_validation-output.txt")
#print("training df")
#print("cols")
#print(ncol(df.training[,-c(1:3)]))
#print("rows")
#print(nrow(df.training[,-c(1:3)]))
#print("cols 1-3")
#print(df.training[,c(1:3)])
#print("#####################################")
#printf("size validation df")
#print("cols")
#print(ncol(df.validation[,-c(1:3)]))
#print("rows")
#print(nrow(df.validation[,-c(1:3)]))
#print("cols 1-3")
#print(df.validation[,c(1:3)])
#print("#####################################")
#print("op.fit:")
#print(op.fit)
#print("#####################################")
#print("lambda.min for cv.fit:")
#print(cv.fit$lambda.min)
#sink()
#######

##### SAVING GENERATED MODEL op.fit as RDS ##### 

saveRDS(op.fit, model_rds_path)
printf("Model has been saved at: %s", model_rds_path)
}
#### USING GENERATED MODEL op.fit #####
printf("Model has been loaded from: %s", model_rds_path)

# loading model and extracting active coefficients
op.fit = readRDS(model_rds_path)
jpeg("lambda_plot_full")
plot(op.fit)


Coefficients<-coef(op.fit,s=op.fit$lambda.min)
op.index<-which(Coefficients!=0)
p.training = predict(op.fit,newx = df.training[,-c(1:3)],s="lambda.min",type="response")
p.validation = predict(op.fit,newx = df.validation[,-c(1:3)],s="lambda.min",type="response")

#calculat threshold for cutoff
# "low Risk" < treshold > "high risk"
#scale RR between 0->1

scale.prediction = (p.training-min(0))/(max(p.training)-min(0))
scale.prediction[which(scale.prediction<0)] = 0 


#performance test for set of thresholds

threshold = seq(0,1,0.01)
predictions.roc = data.frame()

#safe all prediction cutoffs as df  

for (it in 1:length(threshold)){
  newline = findInterval(scale.prediction,threshold[it])
  predictions.roc = bind_cols(as.data.frame(newline),predictions.roc)
  names(predictions.roc)[1]= threshold[it]
}

#mirror df to set index right 
predictions.roc= predictions.roc[,ncol(predictions.roc):1]

#callculate error
#true prositiv = prediction 1 & relapsstatus 1
#calculate ssensitivity and false positiv rate 

all.spec = all.sens = all.FP.rate = AUC = vector()
for (i in 1:length(threshold)){
  TP = length(which(predictions.roc[,i]==1 & df.training[,2]==1))
  TN = length(which(predictions.roc[,i]==0 & df.training[,2]==0))
  FP = length(which(predictions.roc[,i]==1 & df.training[,2]==0))
  ALLP = length(which(df.training[,2]==1))
  ALLN = length(which(df.training[,2]==0))
  SENS = TP/ALLP
  SPEC = TN/ALLN
  all.spec = c(all.spec, SPEC)
  all.sens = c(all.sens,SENS)
  FP.rate = FP/ALLN
  all.FP.rate = c(all.FP.rate,FP.rate)
  AUC = c(AUC,auc(roc(df.training[,2],predictions.roc[,i])))
}

#safe  ROC plot as pdf

if (output.model == TRUE){
  pdf(roc_pdf_path)
  par(mfrow=c(2,2),pty = "s")
#plot Roc Kurve
  plot(all.FP.rate,all.sens,type = "l",ylab = "Sensitivity",
  xlab = "False positiv Rate", ylim = c(0,1), xlim = c(0,1),main = "ROC Kurv")
  dev.off()
}


###########################################################################
# #treshold by log-rang test. find "optimal" p-value for log rang
# p.value = c()
# for(i in 1:length(names(predictions.roc))){
#   df.new = data.frame(df.training[,1],df.training[,2],predictions.roc[,i])
#   names(df.new)[1] = "TIME"
#   names(df.new)[2] = "Status"
#   names(df.new)[3] = "PREDICTION"
#   km.type=survfit(Surv(df.new$TIME,df.new$Status) ~ df.new$PREDICTION,
#                   data = df.new,
#                   type="kaplan-meier")
#   tmp = surv_pvalue(km.type, method = "1")$pval
#   p.value = c(p.value,tmp)
# }
# p.value[which(is.na(p.value)== TRUE)] = 1
# op.thresh = threshold[which(p.value == min(p.value))[1]]
##########################################################################

#finde best threshold
#fp = 0 & sens <=90

# saving all FP-rates, SENS, SPEC and associated thresholds
if (output.model == TRUE){
  comp.all = data.frame(data.frame(threshold*max(p.training)),data.frame(all.sens),data.frame(all.spec), data.frame(all.FP.rate))
  names(comp.all)[1] = "threshold"
  names(comp.all)[2] = "SENS"
  names(comp.all)[3] = "SPEC"
  names(comp.all)[4] = "FP-rate"
  write.xlsx(comp.all,comp_all_path)
}
# sink since xlsx cannot be opened via libre
#sink(file="/home/felix/AG_Baumgrass/Results/Good_BCR/Comparison_all.txt")
#print(comp.all)
#sink()

temp = which(all.sens[which(all.FP.rate <= FP.value)]>=SENS.value)#[1]
op.thresh = which(all.sens == all.sens[which(all.FP.rate <= FP.value)][1])
op.thresholds = op.thresh
printf("Positions that have been selected as possible thresholds: %s", op.thresholds)

## selected threshold via position selected.thresh
op.thresh = threshold[op.thresh][selected.thresh]

## saved all found thresholds in op.thresholds
op.thresholds = threshold[op.thresholds]

printf("Varied threshold Values: ")
print(op.thresholds)
printf("RR over %s are interpret as Status 1(Hight Risk) and under as Status 0(Low Risk)", op.thresholds*max(p.training))
printf("Threshold that has been selected %s at position %s", op.thresholds[selected.thresh], selected.thresh)

#pediction RR Values aprox Relapsstaus with calculates Threshold
predicted.validation=findInterval(p.validation/max(p.training),op.thresh)
#print(predicted.validation)
predicted.training=findInterval(p.training/max(p.training),op.thresh)
#print(predicted.training)

####################################################################

#find error in Prediction.validation for DDPR
fehler.ddpr = df.validation[,2]-df.validation[,3]
fehler.ddpr[fehler.ddpr == 1]="FN"
fehler.ddpr[fehler.ddpr == -1]="FP"
fehler.ddpr[fehler.ddpr==0]=""

#find error in Prediction model
fehler.model = df.validation[,2]-predicted.validation
fehler.model[fehler.model==1]="FN"
fehler.model[fehler.model==-1]="FP"
fehler.model[fehler.model==0]=""

#compare Real Status / DDPR Status / Model Status 
vergleich.validation = data.frame(rownames.validation,df.validation[,2],df.validation[,3],fehler.ddpr,predicted.validation,fehler.model)
names(vergleich.validation)[1] = "Patien ID"
names(vergleich.validation)[2] = "Real Status"
names(vergleich.validation)[3] = "DDPR Status"
names(vergleich.validation)[4] = "Fehler DDPR"
names(vergleich.validation)[5] = "Predicted Status Model"
names(vergleich.validation)[6] = "Fehler Model"
if (output.thresholds == TRUE){
  write.xlsx(vergleich.validation, comp_validation_path)
}

#find error in Prediction.training for DDPR
fehler.ddpr = df.training[,2]-df.training[,3]
fehler.ddpr[fehler.ddpr == 1]="FN"
fehler.ddpr[fehler.ddpr == -1]="FP"
fehler.ddpr[fehler.ddpr==0]=""

#find error in prediction model
fehler.model = df.training[,2]-predicted.training
fehler.model[fehler.model==1]="FN"
fehler.model[fehler.model==-1]="FP"
fehler.model[fehler.model==0]=""

#compare Real Status / DDPR Status / Model Status
vergleich.training = data.frame(rownames.training,df.training[,2],df.training[,3],fehler.ddpr,predicted.training,fehler.model)
names(vergleich.training)[1] = "Patien ID"
names(vergleich.training)[2] = "Real Status"
names(vergleich.training)[3] = "DDPR Status"
names(vergleich.training)[4] = "Fehler DDPR"
names(vergleich.training)[5] = "Predicted Status Model"
names(vergleich.training)[6] = "Fehler Model"
if (output.thresholds == TRUE){
  write.xlsx(vergleich.training, comp_training_path)
}

#bind training and validation error
vergleich.total = bind_rows(vergleich.training,vergleich.validation)
if (output.thresholds == TRUE){
  write.xlsx(vergleich.total, comp_total_path)
}

## saving threshold data in txt
if (output.thresholds == TRUE){
  sink(file = thresholds_path)
  printf("Selected Parameters for SENS >= %s and alpha = %s", FP.value, SENS.value, alpha.value)
  printf("Possible thresholds for this configuration: %s", op.thresholds)
  printf("Specifically selected Threshold for this configuration: %s at position: %s", op.thresholds[selected.thresh]*max(p.training), selected.thresh)
  printf("RR for this threshold: %s", op.thresh*max(p.training))
  printf("AUC Value: %s", auc(roc(vergleich.total[,2],vergleich.total[,5])))
  sink()
}
##############################################################################
#AUC
if (output.thresholds == TRUE){
  pdf(auc_pdf_path)
  plot(roc(vergleich.total[,2],vergleich.total[,5]), )
  #aucval = auc(roc(vergleich.total[,2],vergleich.total[,5]))
  #text(3,2,paste0("AUC: ",aucval))
  

  
}
print(auc(roc(vergleich.total[,2],vergleich.total[,5])))


##############################################################################
#get active coeffizienten !!
Coef.names<-colnames(df.training[,-c(1:3)])

#get names from coef.names for all active index in op.fit
all.activ.names <- colnames(df.training[,-c(1:3)])[op.index]
all.activ.names.split <- strsplit(all.activ.names[1:length(all.activ.names)],".",fixed = TRUE)

x = vector()
y = vector()
z = vector()
modus = vector()
quat = vector()
i=0
while (i < length(all.activ.names)){
  i=i+1
  x= c(x,all.activ.names.split[[i]][1])
  y= c(y,all.activ.names.split[[i]][2])
  z= c(z,all.activ.names.split[[i]][3])
  modus= c(modus,all.activ.names.split[[i]][4])
  quat = c(quat,all.activ.names.split[[i]][5])
}

dev.off()

########################################################################
# create iAUC

if (output.thresholds == TRUE){
pdf(iauc_pdf_path)
iAUC = AUC.uno(Surv(df.training[,1],df.training[,2]),Surv(df.validation[,1],df.validation[,2]),p.validation,times=seq(10,1000,10), savesensspec=TRUE)
names(iAUC)
iAUC$iauc
plot(iAUC)
dev.off()
}

########################################################################

#safe coef names and values for op.fit as excel
df.training = as.data.frame(df.training)
df.validation = as.data.frame(df.validation)

#create heatmap for op.fit coef
df.head.training = df.training[which(colnames(df.training) %in% all.activ.names)]
df.head.validation = df.validation[which(colnames(df.validation) %in% all.activ.names)]
df.head = bind_rows(df.head.training,df.head.validation)
row.names(df.head) = c(rownames.training,rownames.validation)

#Relapsstatus as Vector Red/blue for Heatmap
relaps = as.numeric(c(df.training$`Relaps Status`,df.validation$`Relaps Status`))
relaps[relaps == 1] = "red"
relaps[relaps == 0] = "blue"

#calculate mean and var for training and validation 
all.relaps = cohort$`Patient ID`[which(cohort$`Relapse Status`== "Yes")]
all.relapsfree = cohort$`Patient ID`[which(cohort$`Relapse Status`== "No")]
all.mean.range.relaps = all.var.range.relaps = all.mean.range.relapsfree = all.var.range.relapsfree = vector()

it = 0
for (i in 1:length(names(df.head))){
  it = it +1
  all.mean.range.relaps = c(all.mean.range.relaps,mean(df.head[which(row.names(df.head) %in% all.relaps),it]))
  all.var.range.relaps = c(all.var.range.relaps,var(df.head[which(row.names(df.head) %in% all.relaps),it]))
  all.mean.range.relapsfree = c(all.mean.range.relapsfree,mean(df.head[which(row.names(df.head) %in% all.relapsfree),it]))
  all.var.range.relapsfree = c(all.var.range.relapsfree,var(df.head[which(row.names(df.head) %in% all.relapsfree),it]))
}

#safe all real absRange values 
#transponse df.head 
new.head = t(df.head)
colnames(new.head) = row.names(df.head)
row.names(new.head) = c(1:dim(new.head)[1])

#safe df as excel
all.right.index <- op.index
result.data = data.frame(all.right.index,x,y,z,modus,quat,
                         all.mean.range.relaps,all.var.range.relaps,
                         all.mean.range.relapsfree,all.var.range.relapsfree,
                         new.head)
names(result.data)[10] = "Var(AbsRange) Relapsfree"
names(result.data)[9] = "Mean(AbsRange) Relapsfree"
names(result.data)[8] = "Var(AbsRange) Relaps"
names(result.data)[7] = "Mean(AbsRange) Relaps"
names(result.data)[6] = "Quadrant"
names(result.data)[5] = "Modus"
names(result.data)[4] = "z Variable"
names(result.data)[3] = "Y Variable"
names(result.data)[2] = "x Variable"
names(result.data)[1] = "Variable Index"
if (output.model == TRUE){

write.xlsx(result.data, variables_path) 
}


##########################################################################

#create heatmap for all active coef.
df.training = as.data.frame(df.training)
df.validation = as.data.frame(df.validation)
df.head.training = df.training[which(colnames(df.training) %in% names(df.training[,-c(1:3)])[op.index])]
df.head.validation = df.validation[which(colnames(df.validation) %in% names(df.validation[,-c(1:3)])[op.index])]
df.head = bind_rows(df.head.training,df.head.validation)
row.names(df.head) = c(rownames.training,rownames.validation)
df.head = bind_cols(as.data.frame(relaps),df.head)
row.names(df.head) = c(rownames.training,rownames.validation)
df.head = df.head[order(df.head$relaps),]

if (output.model == TRUE){
  pdf(heat_path , width = 10)
  heatmap(t(as.matrix(df.head[,-1])),scale = "none",Colv = NA, ColSideColors = relaps[order(relaps)])
  dev.off()
}



#create kaplan-meier - survival Kurve
if (output.model == TRUE) {
  time = c(df.training[,1],df.validation[,1])
  status = c(df.training[,2],df.validation[,2])
  df.new = data.frame(time,status,vergleich.total[,5])
  #df.new.ddpr = data.frame(time,status,vergleich.total[,3])
  km.type=survfit(Surv(df.new[,1],df.new[,2]) ~ df.new[,3],data = df.new,type="kaplan-meier")
  #km.type.ddpr=survfit(Surv(df.new.ddpr[,1],df.new.ddpr[,2]) ~ df.new.ddpr[,3],data = df.new,type="kaplan-meier")
  survplot = ggsurvplot(km.type, conf.int = TRUE, legend.labs = c("low Risk","high Risk"),ggtheme = theme_minimal(),pval = TRUE,pval.method = TRUE, risk.table=TRUE)
  #survplot.ddpr = ggsurvplot(km.type.ddpr, conf.int = TRUE, legend.labs = c("low Risk","high Risk"),ggtheme = theme_minimal(),pval = FALSE,pval.method = TRUE, risk.table=TRUE)
  ggsave(file = kaplan_pdf_path, print(survplot))
  #ggsave(file = kaplan.ddpr_pdf_path, print(survplot.ddpr))
  #print("df.new")
  #print(df.new)
  #print(dim(df.new))

}

printf("AUC Value: %s", auc(roc(vergleich.total[,2],vergleich.total[,5])))

