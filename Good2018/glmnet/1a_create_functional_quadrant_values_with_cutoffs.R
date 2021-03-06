#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high ristk

rm(list = ls())



############### initiate
# cofactor            1, 0.2, 0.1 [so far at 0.2]
# mincells            2,5,10,20 [so far at 5]
# stat.info           absRange, variance, freq_green, mean
# cluster.size        1-12
# trimming            TRUE, FALSE [so far FALSE]
# project.name        Basal, BCR, IL7, Pervanadate, TSLP
# project.name.long   Basal, BCR-Crosslink, IL-7, Pervanadate, TSLP 
# subgroup            Validation, Training
# coverage            full, func [so far func(tional)]
cofactor = 0.2
mincells = 5
stat.info = "absRange"
# stat.info = "variance"
# stat.info = "freq_green"
# stat.info = "mean"
cluster.size = 6
project.name = "TSLP"
project.name.long = "TSLP"

subgroup = "Training"
# subgroup = "Validation"
coverage = "func"
### set paths
folder.path = "/scratch/drfz/projects/Good2018/glmnet"
setwd(folder.path)
################


### load libraries
# xlsx            :   Library for Excel reading and creating
# RSQLite         :   interact with database
# reshape2        :   function melt()
# dplyr           :   faster binding of columns bind_cols() and rows bind_rows()
# foreach         :   allows for each loops
# doParallel      :   use several clusters for parallele calculations
libraries = c("RSQLite","dplyr","reshape2","foreach","doParallel")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/projects/Good2018/glmnet/PRI_funs.R")


### metatable file
# cohort_full=read.csv(sprintf("/scratch/drfz/projects/Good2018/%s_%s_patient_cohort.csv",project.name,subgroup))
# cohort = cohort_full[,c(1,5,8,11,16,15)]
sub.set = read.csv(sprintf("/scratch/drfz/projects/Good2018/%s_%s_patient_cohort.csv",project.name,subgroup))
sub.set$Relapse.Status= factor(sub.set$Relapse.Status)
sub.set$Patient.ID = factor(sub.set$Patient.ID)
col.vec.func = as.vector(unlist(read.table(file="/scratch/drfz/projects/Good2018/columns_functional.txt")))

# load data base -----------------------------------------
db.path=file.path("","data","databases")
db.name="RB_20191002_Good2018.sqlite3"

fcs$connectDb(file.path(db.path,db.name))
### select project
project.idx = which(dbListTables(fcs$conn)==project.name)
fileID = fcs$getDFtable(paste0(project.name,"_fileIdentity"))
stainID = fcs$getDFtable(paste0(project.name,"_markerIdentity"))
dbDisconnect(fcs$conn)

### and according file name
# read from dataframe RDS-file
data.table.name = sprintf("%s/Rdata/%s_%s_%s_df_cof%s",folder.path,project.name,subgroup,coverage,cofactor)
# save in quadrant table RDS-file
quad.table.name = sprintf("%s/Rdata/%s_%s_%s_quadrant_%s_cof%s",folder.path,project.name,subgroup,coverage,stat.info,cofactor)



### load matrix temp.data.all
### load temp.data.all
temp.data.all = readRDS( file = paste0(data.table.name,".rds"))
print(dim(temp.data.all))
#[1] 470928     40


## ----triplots quadrants, cache=TRUE-----------------------------------------
col.func.idx = which(names(temp.data.all) %in% col.vec.func)
temp.data.all = temp.data.all[,c(1,col.func.idx)]
len.var = ncol(temp.data.all)-1
colvec = colnames(temp.data.all)[2:ncol(temp.data.all)]
len.col = length(colvec)

### initate
it = 0
quad.df = data.frame( matrix(
  ,nrow = nrow(sub.set)
  ,ncol = 9792  # m=18, with ncol = m * (m-1) * (m-2) * 0.5 * 4
  # ,ncol = 2448  # m=18, with ncol = m * (m-1) * (m-2) * 0.5
  ,byrow = TRUE
), stringsAsFactors = FALSE)

cl <- makeCluster(cluster.size)
registerDoParallel(cl)

ptm <- proc.time()
quad.sample_id = vector()
for ( i in 1:nrow(sub.set)) {
# for ( i in 1:1 ) {
  quad.sample_id = c(quad.sample_id,as.character(sub.set$Patient.ID[i]))
  
  file.name = paste0(sub.set$Patient.ID[i],"_",project.name.long,".fcs")
  file.idx = fileID$file_ID[which(fileID$filename==file.name)]
  if (length(file.idx)==0) {
    file.name = paste0(sub.set$Patient.ID[i],"_",tolower(project.name.long),".fcs") # file names vary..
    file.idx = fileID$file_ID[which(fileID$filename==file.name)]
  }
  cutoffs = stainID$file_savedCutoffs[which(stainID$file_ID==file.idx)]
  names(cutoffs) = stainID$shortname[which(stainID$file_ID==file.idx)]
  cutoffs = cutoffs[col.func.idx-1]
  
  quad.file = vector()
  for ( v1 in 1:(len.col-1) ) {
  # for ( v1 in 1:1 ) {
    it = it +1
    
    quad.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
    # quad.oper <- foreach ( v2=11:14,.combine=cbind )  %dopar% {
      quadrant.vec = vector()
      
      for ( v3 in 1:len.col ) {
      # for ( v3 in 13:13 ) {
        if ( all(v3 != c(v1,v2)) ) {
          sampl.data = temp.data.all[which(temp.data.all$file_id==as.vector(sub.set$Patient.ID[i]))
                                     ,c(colvec[v1],colvec[v2],colvec[v3])]
          ### NEW::ONLY rows where used if v1 or v2 are >0
          # sampl.data = sampl.data[which(sampl.data[,1]>=0 & sampl.data[,2]>=0),]
          
          # calculate triplot quadrants ---------------------------------------------
          quad.results = fcs$calc_triplot_quadrant(temp.data = sampl.data, 
                                                   calc.meth = stat.info, 
                                                   min.cells = mincells,
                                                   prod.cutoff = cutoffs[c(v1,v2)])
          quadrant.vec = c(quadrant.vec, quad.results)
        }
      }
      
      return(quadrant.vec)
    }
    
    quad.file = c(quad.file,as.vector(quad.oper))
    
    # create label vector in the last step ------------------------------------------
    if ( i==nrow(sub.set) ) {
      
      label.file = vector()
      for ( v1 in 1:(len.col-1) ) {
      # for ( v1 in 1:1 ) {

        label.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
        # label.oper <- foreach ( v2=11:14,.combine=cbind )  %dopar% {
          label.vec = vector()
          
          for ( v3 in 1:len.col ) {
          # for ( v3 in 13:13 ) {
            if ( all(v3 != c(v1,v2)) ) {
              
              if (stat.info=="zRange") {
                label.vec = c(label.vec, paste0(colvec[v1],".",colvec[v2],".",colvec[v3]))
              } else {
                label.vec = c(label.vec,
                              paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q1"),
                              paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q2"),
                              paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q3"),
                              paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q4")
                )
              }
            }
          }
          
          return(label.vec)
        }
        label.file = c(label.file,as.vector(label.oper))
      }
    }
    
    printf("%s::%s::%s::%s::v1=%s[%s/%s] ready (it=%s)", i, sub.set$Patient.ID[i], subgroup,stat.info, colvec[v1],v1,len.col, it)
    print(proc.time() - ptm)
  } 
  quad.df[i,] = quad.file
  
  printf("File %s ready (it=%s)", i, it)
  print(proc.time() - ptm)
}
print("Total time:")
print(proc.time() - ptm)

colnames(quad.df) = label.file
rownames(quad.df) = quad.sample_id


stopCluster(cl)

saveRDS(quad.df, file=paste0(quad.table.name,".rds"))
printf("%s saved!",quad.table.name)




