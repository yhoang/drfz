#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019
# Goood2018

####
# SKIP STEP 1_, GO TO STEP 2a_ AFTERWARDS
####

### Basal representatives
# UPN12/15 low risk
# UPN10/22 high ristk

rm(list = ls())



############### initiate
cofactor = 0.2
mincells = 5
stat.info = "absRange"
stat.info = "freq_green"
cluster.size = 6
subgroup = "Training"
# subgroup = "Validation"
coverage = "func"
condi = 1.1
### set paths
folder.path = "/scratch/drfz/Good2018/glmnet"
setwd(folder.path)
################


### load libraries
# xlsx            :   Library for Excel reading and creating
# RSQLite         :   interact with database
# reshape2        :   function melt()
# dplyr           :   faster binding of columns bind_cols() and rows bind_rows()
# foreach         :   allows for each loops
# doParallel      :   use several clusters for parallele calculations
libraries = c("xlsx","RSQLite","dplyr","reshape2","foreach","doParallel")
lapply(libraries,require, character.only = TRUE)
### load functions
source("/scratch/drfz/Good2018/PRI_funs.R")


### metatable file
cohort_full=readxl::read_excel("/scratch/drfz/Good2018/patient_cohort.xlsx",1,col_names=TRUE,progress=T)
cohort = cohort_full[,c(1,5,8,11,16,15)]
sub.set = cohort[which(cohort$Cohort==subgroup),1:5]
sub.set$`Relapse Status`= factor(sub.set$`Relapse Status`)
sub.set$`Patient ID` = factor(sub.set$`Patient ID`)

## load data base -----------------------------------------
db.path=file.path("","data","databases")
db.name="RB_20191002_Good2018.sqlite3"

fcs$connectDb(file.path(db.path,db.name))
### select project
project.name = "Basal"
project.idx = which(dbListTables(fcs$conn)==project.name)
fileID = fcs$getDFtable(paste0(project.name,"_fileIdentity"))
stainID = fcs$getDFtable(paste0(project.name,"_markerIdentity"))
dbDisconnect(fcs$conn)

### and according file name
data.table.name = sprintf("%s/Rdata/%s_%s_full_df_cof%s",folder.path,project.name,subgroup,cofactor)
quad.table.name = sprintf("%s/Rdata/%s_%s_%s_quadrant_%s_condi%s_cof%s",folder.path,project.name,subgroup,coverage,stat.info,condi,cofactor)
sample.cond.name = sprintf("%s/%s_%s_%s_condi%s_cof%s.txt",folder.path,project.name,subgroup,coverage,condi,cofactor)




### load matrix temp.data.all
### load temp.data.all
temp.data.all = readRDS( file = paste0(data.table.name,".rds"))
print(dim(temp.data.all))
#[1] 3697393      40    for Training data
#[1] 458229     19      for Validation data

## --- marker coverage full or functional ---------------------
if (coverage=="full") {
  columns_name = "sorted"
} else columns_name = "functional"
col.vec.func = as.vector(unlist(read.table(file=sprintf("/scratch/drfz/Good2018/columns_%s.txt", columns_name))))

## --- select functional markers -----------------------------------------
col.func.idx = which(names(temp.data.all) %in% col.vec.func)
temp.data.all = temp.data.all[,c(1,col.func.idx)]
len.var = ncol(temp.data.all)-1
colvec = colnames(temp.data.all)[2:ncol(temp.data.all)]
len.col = length(colvec)

## ---- select samples which follows condition 1.1 ----------------------------
# %(CD34+CD38)>0.2
# %(TdT)>0.2
if (FALSE) {
  # frequencies = c(0.02,0.04,0.06,0.08,0.1)
  # frequencies = c(0.5,1,2,10)
  frequencies = c(0.02)
  ##### take freqs=c(0.2,0.2)
  freq.df = data.frame(matrix(NA,
                              nrow=length(frequencies),
                              ncol=length(frequencies)
  ), stringsAsFactors = FALSE)
  rownames(freq.df) = paste0(rep("CD34.CD38=",length(frequencies)),frequencies)
  colnames(freq.df) = paste0(rep("TdT=",length(frequencies)),frequencies)
  
  for (a in 1:length(frequencies)) {
      for (b in 1:length(frequencies)) {
        cond.samples = vector()
        for ( i in 1:nrow(sub.set)) {
          sampl.data = temp.data.all[which(temp.data.all$file_id==sub.set$`Patient ID`[i]),]
          
          file.name = paste0(sub.set$`Patient ID`[i],"_",project.name,".fcs")
          file.idx = fileID$file_ID[which(fileID$filename==file.name)]
          if (length(file.idx)==0) {
            file.name = paste0(sub.set$`Patient ID`[i],"_",tolower(project.name),".fcs") # file names vary..
            file.idx = fileID$file_ID[which(fileID$filename==file.name)]
          }
          cutoffs = stainID$file_savedCutoffs[which(stainID$file_ID==file.idx)]
          names(cutoffs) = stainID$shortname[which(stainID$file_ID==file.idx)]
          cutoffs = cutoffs[col.func.idx-1]
          
          if (fcs$condition_approval(df.sample = sampl.data,
                                     cutoffs = cutoffs,
                                     condition=condi,
                                     freqs = c(frequencies[a],frequencies[b]))) cond.samples = c(cond.samples,as.character(sub.set$`Patient ID`[i]))
        }
        freq.df[a,b] = length(cond.samples)
        printf("%s/%s with %%(CD34/CD38)=%s %%(TdT)=%s",length(cond.samples),nrow(sub.set),frequencies[a],frequencies[b])
        #a rows, b columns
        
      }
  }
  ##### for all frequency combinations
  # write.xlsx(freq.df,file=sprintf("frequency_samples_condi%s.xlsx",condi),sheetName = sprintf("CD34.CD38 and TdT %s",subgroup))
  ##### take freqs=c(0.2,0.2)
  write.table(cond.samples,sample.cond.name,row.names = F,col.names = F)
} else {
  cond.samples = as.vector(unlist(read.table(sample.cond.name)))
}
printf("%s out of %s %s samples to conduct",length(cond.samples),nrow(sub.set),subgroup)

# now:
# fixed X/Y with TdT, pSTAT5, CD24
marker.axes = c("TdT","pSTAT5","CD24")
marker.axes.idx = which(colvec %in% marker.axes)
############################################################################




## ----triplots quadrants, cache=TRUE-----------------------------------------
### initate
it = 0
quad.df = data.frame( matrix(
  ,nrow = length(cond.samples)
  ,ncol = 3072  # m=18, with ncol = (m-2) * (m-6) * (m-10) * 0.5 * 4
  ,byrow = TRUE
), stringsAsFactors = FALSE)


#cluster.size = 10
cl <- makeCluster(cluster.size)
registerDoParallel(cl)

ptm <- proc.time()
for ( i in 1:length(cond.samples)) {
  file.name = paste0(cond.samples[i],"_",project.name,".fcs")
  file.idx = fileID$file_ID[which(fileID$filename==file.name)]
  if (length(file.idx)==0) {
    file.name = paste0(cond.samples[i],"_",tolower(project.name),".fcs") # file names vary..
    file.idx = fileID$file_ID[which(fileID$filename==file.name)]
  }
  cutoffs = stainID$file_savedCutoffs[which(stainID$file_ID==file.idx)]
  names(cutoffs) = stainID$shortname[which(stainID$file_ID==file.idx)]
  cutoffs = cutoffs[col.func.idx-1]
  
  quad.file = vector()
  for ( v1 in 1:(len.col-1) ) {
    it = it +1
    
    quad.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
      quadrant.vec = vector()
      
      for ( v3 in 1:len.col ) {
        ### add condition=1 with TdT/pSTAT5/CD24 in the axes X/Y
        if ( !any(c(v1,v2) %in% marker.axes.idx) ) next
        
        if ( all(v3 != c(v1,v2)) ) {
          sampl.data = temp.data.all[which(temp.data.all$file_id==cond.samples[i])
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
    if ( i==length(cond.samples) ) {
      
      label.file = vector()
      for ( v1 in 1:(len.col-1) ) {
        
        label.oper <- foreach ( v2=(v1+1):len.col,.combine=cbind )  %dopar% {
          label.vec = vector()
          
          for ( v3 in 1:len.col ) {
            
            if ( !any(c(v1,v2) %in% marker.axes.idx) ) next
            
            if ( all(v3 != c(v1,v2)) ) {
              label.vec = c(label.vec,
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q1"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q2"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q3"),
                            paste0(colvec[v1],".",colvec[v2],".",colvec[v3],".",stat.info,".","Q4")
              )
            }
          }
          
          label.vec
        }
        label.file = c(label.file,as.vector(label.oper))
      }
    }
    
    printf("%s::%s::quadrants::%s::%s::v1=%s[%s/%s] ready (it=%s)", i, cond.samples[i], stat.info, subgroup,colvec[v1],v1,len.col, it)
    print(proc.time() - ptm)
  } 
  quad.df[i,] = quad.file
  
  printf("File %s ready (it=%s)", i, it)
  print(proc.time() - ptm)
}
print("Total time:")
print(proc.time() - ptm)

colnames(quad.df) = label.file
rownames(quad.df) = cond.samples


stopCluster(cl)

saveRDS(quad.df, file=paste0(quad.table.name,".rds"))
printf("%s saved!",quad.table.name)
dim(quad.df)




