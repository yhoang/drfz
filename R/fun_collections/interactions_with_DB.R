#!/usr/bin/R
# author: Yen Hoang
# DRFZ 2019

### uncomment if not created yet
# fcs = new.env()

### connect to database PRI-base
fcs$connectDb<- function(fname) {
  this=fcs
  this$conn <- dbConnect(SQLite(), dbname = fname)
  print(paste("Database opened:",fname))
}


### disconnect to database PRI-base
fcs$disconnectDb <- function () {
  this=fcs
  print(paste("Database closed:",this$conn))
  dbDisconnect(this$conn)
}


######## LOAD VARIABLES FROM DATABASE ##############
##### @param
# project.idx   project index
# file.idx      file index
##### returns data table
fcs$loaddata <- function(
  project.idx,
  file.idx)  
{
  this=fcs
  
  ### list all tables in database
  this$table.list=dbListTables(this$conn)
  
  ### get every index which has metadata
  metadata.idx=grep("markerIdentity|colnameIndex|fileIdentity|fileIndex|UserHistory|Classification|equipmentInfo|fileComments|SPILL",this$table.list)
  
  this$metadata.list = vector()
  df.num = 0
  for (i in metadata.idx){
    df.num = df.num + 1 
    this$metadata.list[df.num] = this$table.list[i]
  }
  
  ### now only the table names are listed which have our data (fluorescence intensity for each feature)
  this$project.list = this$table.list[-metadata.idx]
  ### if there are several projects in the database, choose one
  this$current.project = this$project.list[project.idx]
  
  this$current.filetable = this$getDFtable(paste0(this$current.project,"_fileIdentity"))
  this$file.list = this$current.filetable[,2]
  
  this$current.file = this$file.list[file.idx]
  this$current.staintable = this$getDFtable(paste0(this$current.project,"_markerIdentity"))
  this$current.vars = this$getVariables(index=file.idx)
}

######## GET SAMPLE DATA FROM FILE INDEX IN DATABASE #############
##### @param
# table       table name
# fileidx     file index number
# columns     columns vector (optional)
# stain       staining vector (optional)
# cofactor    cofactor for asinh transformation (default = 1)
fcs$getData <- function (table,fileidx, columns=NA, stain=NA, cofactor = NA) {
  this=fcs
  
  data=dbGetQuery(this$conn,paste("SELECT * FROM ",table," WHERE file_ID == '",fileidx,"'",sep=""))
  
  # and ignore columns with NAs
  col.NA=NA
  for (i in 1:ncol(data)) {
    if (any(is.na(data[,i])))  col.NA=c(col.NA,i)
  }
  if (!is.na(col.NA)) data=data[,-col.NA]
  
  # change column names
  column.names=colnames(data)
  if (!is.na(columns)) {
    data=data[,columns+1]
  }
  
  if ( !is.na(stain) & !is.na(columns) ) { colnames(data) = stain[columns]}
  if ( is.na(stain) & !is.na(columns) ) { colnames(data) = this$selected.vars[columns] }
  if ( !is.na(stain) & is.na(columns) ) { colnames(data) = c("file_ID",stain) }
  if ( is.na(stain) & is.na(columns) ) { colnames(data) = c("file_ID",this$selected.vars) }
  
  # set asinh cofactor to 1 if not set
  if ( is.na(cofactor) ) cofactor = 1
  
  if ( !is.na(columns) ) {
    data=asinh(data/cofactor)
  } else {
    data=asinh(data[,(2:dim(data)[2]/cofactor)])              
  }
  
  printf("w: do getData(%s) from table='%s' with fileidx=%s and asinh cofactor=%s",nrow(data),table,fileidx,cofactor)
  
  data
}

######## GET FULL TABLE FROM DATABASE #############################
##### @param
# table     table name
fcs$getDFtable <- function (table) {
  this=fcs
  
  table.df=dbGetQuery(this$conn,paste("SELECT * FROM ", table))
  
  return(table.df)
}


######## GET FEATURE SHORT NAMES ###################################
##### @param
# staintable    manual stain table (optional)
# index         file index (optional)    
fcs$getVariables <- function (staintable=NA,index=NA) {
    this=fcs

    if (is.na(staintable)) staintable=this$current.staintable
    
    if (is.na(index)) {
        if (exists("current.file",envir=fcs)) {
            file=tclvalue(tkget(this$tkchoosefile))
        } else {
            file=this$file.list[1]
        }
        index=this$current.filetable[which(this$current.filetable[,2]==file),1]
    }
    printf("do getVariables from table=%s: file idx=%s",this$current.project,index)

    vars=staintable[which(staintable[,1]==index),4]
    
    #printf("vars=%s",paste(vars,collapse=" "))
    print(vars)
    cat("Done getVariables\n\n")

    vars
}

