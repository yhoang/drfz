
#wichtige pakete 
libarys = c("survival","glmnet","readxl","dplyr","doParallel","xlsx")
lapply(libarys,require,character.only = TRUE)

#printfunktion
printf <- function(...) invisible(print(sprintf(...)))

#laden der Patientendaten aus dem Goodpaper
patient_data <- read_excel("~/Good/patient_data.xlsx") # laden der daten
#View(patient_data)

#bearbeiten der Daten 
#NA's in collum time auf 0 setzen
#alle Daten sind vector(chr) daher richtigen datentyp hinterlegen
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

#Training/Validation Set
#nur n�tige klinische colum aus den patienten daten
cohort=patient_data[,c(1,2,3,11,15,17)]


training.set=cohort[which(cohort$Cohort=="Training"),1:6]
validation.set=cohort[which(cohort$Cohort=="Validation"),1:6]

total.set=bind_rows(training.set,validation.set)


#einlesen der RDS Files training/validation
df.training<-readRDS("~/Good/RDS/TriPlotData/Basal_Training_func_quadrant_freq_green_condi1.1_cof0.2.rds")
df.validation<-readRDS("~/Good/RDS/TriPlotData/Basal_Validation_func_quadrant_freq_green_condi1.1_cof0.2.rds")
#View(df.training)
#View(df.validation)

df.total<-bind_rows(df.training,df.validation)
sample.size=ncol(df.total)
rownames(df.total)=c(rownames(df.training),rownames(df.validation))

#View(df.total)
#von 60 Patienet 54 Training/Validation ,
#1 validation, 12 Training erf�llen nicht condition 1.1
#colum for prediction variable in ds
typeColNum=1

#relapsstatus in training set as numeric
condition = training.set$`Relapse Status`[which(training.set$`Patient ID`%in% row.names(df.training))]
condition[condition %in% c("Yes")]=1
condition[condition %in% c("No")]=0
condition=as.numeric(condition)

#survival time in training.set as numeric
survival=training.set$`Survival Time (Day)`[which(training.set$`Patient ID`%in%row.names(df.training))]
survival=as.numeric(survival)

#age as numeric
age=training.set$`Age at Diagnosis`[which(training.set$`Patient ID`%in%row.names(df.training))]
age=as.numeric(age)

#add condition/survivaltime/age to df.training
df.training=bind_cols(as.data.frame(age),df.training)
names(df.training)[typeColNum] = "Age at Diagnosis"
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
age=validation.set$`Age at Diagnosis`[which(validation.set$`Patient ID`%in%row.names(df.validation))]
age=as.numeric(age)
df.validation=bind_cols(as.data.frame(age),df.validation)
names(df.validation)[typeColNum] = "Age at Diagnosis"
df.validation=bind_cols(as.data.frame(condition),df.validation)
names(df.validation)[typeColNum]="Relaps Status"
df.validation=bind_cols(as.data.frame(survival),df.validation)
names(df.validation)[typeColNum] = "Survivaltime (Day)"


#model ben�tigt eingabe als matrix der variablen
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


#erstellen von surv.Obj , n�tig f�r cox model
sur_obj_training = Surv(df.training[,1],df.training[,2])

#initialisierung der alpha suche im training
alphalist <- seq(0,1,by=0.005)
cluster.size =10
alpha.collect = lambda1se.collect = vector()
seed.vec = sample(10^2)
it.total = 0

cl <- makeCluster(cluster.size)
registerDoParallel(cl)

timeSTART = Sys.time()
ptm <- proc.time()
printf("##### Start at %s.",timeSTART)
#printf("run %s::%s::cluster=%s::samplesize=%s",project.name,stat.info,cluster.size,sample.size)
while (it.total < 100) {
  it.total = it.total + 1  
  if ( it.total %% 10 == 0 ) {
    #printf("find_alpha %s run #%s..",stat.info,it.total)
    print(Sys.time()-timeSTART)
    print(proc.time() - ptm)
  }
  
  set.seed(seed.vec[it.total])
  ### 10-fold with 5-6 samples in one fold
  set.foldid = sample(rep(seq((1/4)*nrow(df.training)),length=nrow(df.training)))
  
  ### no need to set nfolds if foldid is provided since observations are already distributed in folds
  elasticnet <- lapply(alphalist, function(a) {
    cv.glmnet(x = df.training[,-c(1:2)], sur_obj_training,
              alpha=a, family="cox",
              foldid = set.foldid,
              parallel = TRUE
    )
  })
  
  min.err = min.err.lambda.1se = vector()
  for (i in 1:length(alphalist)) {
    #print(min(elasticnet[[i]]$cvm))
    #printf("min(CVmean)=%s lambda.1se=%s",min(elasticnet[[i]]$cvm),elasticnet[[i]]$lambda.1se)
    min.err = c(min.err,min(elasticnet[[i]]$cvm))
    #min.err.lambda.1se = c(min.err.lambda.1se,elasticnet[[i]]$lambda.1se)
  }
  min.err.idx = which(min.err==min(min.err))
  alpha.best = alphalist[min.err.idx]
  #lambda.best = elasticnet[[min.err.idx]]$lambda.1se
  # plot(min.err, ylab= "min(deviance)",
  #      main=sprintf("a=%s,seed=%s,lambda.1se=%s",alpha.best,seed.vec[it.total],lambda.best))
  
  #lambda1se.collect = c(lambda1se.collect, lambda.best)
  alpha.collect = c(alpha.collect, alpha.best)
}
stopCluster(cl)
print(Sys.time())

alpha.tab = table(alpha.collect)
alpha.max = alpha.tab[which(alpha.tab==max(alpha.tab))]


#printf("Best alpha after %s iterations: a=%s", it.total,alpha.max)
print(alpha.tab)