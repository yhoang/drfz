#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

rm(list = ls())
folder.path = "/scratch/drfz/Good2018/"
setwd(folder.path)

### load libraries
# xlsx            :   Library for Excel reading and creating
# RSQLite         :   interact with database
libraries = c("xlsx","RSQLite","dplyr")
lapply(libraries,require, character.only = TRUE)
### load functions
source("PRI_funs.R")

### data base
db.path=file.path("","data","databases")
db.name="RB_20191002_Good2018.sqlite3"

fcs$connectDb(file.path(db.path,db.name))
print(dbListTables(fcs$conn))


### table file
cohort_full=readxl::read_excel("patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort = cohort_full[,c(1,5,8,11,16,15)]

train.set = cohort[which(cohort$Cohort=="Training"),1:5]
valid.set = cohort[which(cohort$Cohort=="Validation"),1:5]
test.set = cohort[which(cohort$Cohort=="NA"),1:5]

project.name = "Basal"
project.idx = which(dbListTables(fcs$conn)==project.name)
fileID = fcs$getDFtable(paste0(project.name,"_fileIdentity"))
stainID = fcs$getDFtable(paste0(project.name,"_markerIdentity"))

mincells = 5

### initiate
data.absRange = data.frame( matrix(
  nrow = 44,
  ncol = 37*3
)
)
label.name = vector()

### go through training set
for ( i in 1:nrow(train.set)) {
  pat.id = train.set[i,1]
  file.name = paste0(pat.id,"_",project.name,".fcs")
  file.idx = fileID$file_ID[which(fileID$filename==file.name)]
  
  ### get data with cofactor = 1
  temp.data = fcs$getData(project.name,file.idx)
  ### get short names
  colnames(temp.data) = stainID$shortname[which(stainID$file_ID==file.idx)]
  
  ### get label names once
  if (i==1) {
    ### remove feature A and B
    col.vec = colnames(temp.data)[-match(c("CD34","TdT"),colnames(temp.data))]
    ### 
    for (name in col.vec) {
      label.name = c(label.name,
                     paste0("min.",name),
                     paste0("max.",name),
                     paste0("absRange.",name)
                     )
    }
  }
  
  ### get ranges
  temp.features = vector()
  for ( j in 1:length(col.vec)) {
    ### get triploT table
    fcs$bintriploT_table(data=temp.data,
                         stainTable = stainID,
                         feat.X = "CD34",
                         feat.Y = "TdT",
                         feat.Z1 = col.vec[j],
                         calc = "MFI",
                         binsize = 0.2,
                         mincells = mincells)
    
    ### fcs$my.calc created
    # fcs$my.calc
    
    ### get absolute bin range of all bins with mincells
    range.C = range(fcs$my.calc$x[which(fcs$my.calc$ncells>=mincells)])
    minVal.C = round(range.C[1],3)
    maxVal.C = round(range.C[2],3)
    absRange = round(diff(range.C),3)
    
    temp.features = c(temp.features,minVal.C,maxVal.C,absRange)
  }
  
  data.absRange[i,] = temp.features
}
### label columns
names(data.absRange) = label.name
### bind training set with bin information
bind = bind_cols(train.set,data.absRange)
### write
write.csv(bind, file=sprintf("%s_%s_%s_absRange_%s.csv",project.name,feat.X,feat.Y,calc.meth))